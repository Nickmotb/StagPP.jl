# TOex.jl — excess oxygen (Oₑₓ) shared between solid, melt and C-O-H fluid at fixed P, T; entry
# point `partition_Oₑₓ`. Sinks, method, sign and units: TOex_overview.md 5.2, TOex_math.md; numbers: TOex_bench.md.

# ---- Section 1 — bulk composition & mass bookkeeping ----

# Mantle end-members ("H"arzburgitic, "B"asaltic) mixed by fraction p of B.
function endmember_composition(p::Real)
    XH     = @SVector [0.4347, 0.4597, 0.0835, 0.0090, 0.0100, 0.0001, 0.0030, 0.0] # mass fraction
    XB     = @SVector [0.5097, 0.0988, 0.0718, 0.1268, 0.1698, 0.0225, 0.0007, 0.0] # mass fraction
    Xox    = @SVector ["SiO2", "MgO", "FeO", "CaO", "Al2O3", "Na2O", "Cr2O3", "O"]
    SymXox = @SVector [:SiO2, :MgO, :FeO, :CaO, :Al2O3, :Na2O, :Cr2O3, :O]
    mmXox  = get_Xoxmm(Xox)
    X      = p*XB + (1-p)*XH   # bulk solid-residue composition, fed to the MAGEMin table
    return (; X, XB, Xox, SymXox, mmXox)
end

# Mass and mole bookkeeping done once before the solve. `TOₑₓ` is the budget in grams.
function mass_budget(Mt, ϕ, TC, TOex, H, XB, mmXox)
    Mf   = ϕ*Mt          # melt tracer mass
    Ms   = Mt*(1-ϕ)      # solid tracer mass
    Mc   = TC*Mt         # carbon tracer mass
    Mh   = H*Mt          # free-water tracer mass
    TOₑₓ = TOex*Mt       # total excess-oxygen budget

    molfracXB   = (XB./mmXox) ./ sum(XB./mmXox)   # mole fractions of the B end-member (sum == 1)
    nB_per_mass = sum(XB./mmXox)                  # moles of B oxides per unit mass mol/g (NOT normalized)
    molMf       = Mf*nB_per_mass                  # moles of the melt's own oxide assemblage

    molCav      = Mc/mm.C                    # moles of available carbon
    # Carbon-limited ceiling on melt XCO₂; 0/0 with neither melt nor carbon, so guarded.
    maxXCO₂_raw = (molMf + molCav) > 0 ? molCav/(molMf + molCav) : 0.0

    # ABSOLUTE GRAMS OF OXYGEN, not a share of the budget — see the header note on units.
    Φ = 2*mm.O*molMf        # converts a melt XCO₂ into grams of excess oxygen

    return (; Mf, Ms, Mc, Mh, TOₑₓ, molfracXB, nB_per_mass, molMf, maxXCO₂_raw, Φ)
end

# ---- Section 2 — solid branch: MAGEMin (Stixrude & Lithgow-Bertelloni 2024) ----

# Isotonic regression by pool-adjacent-violators: the closest non-decreasing sequence in least
# squares. Used to repair a solid oxybarometer that MAGEMin returns non-monotonic.
function _isotonic(y)
    v = Float64[]; w = Int[]
    for yi in y
        push!(v, yi); push!(w, 1)
        while length(v) > 1 && v[end-1] > v[end]     # pool the violating pair
            nv = (w[end-1]*v[end-1] + w[end]*v[end])/(w[end-1] + w[end])
            nw = w[end-1] + w[end]
            pop!(v); pop!(w); v[end] = nv; w[end] = nw
        end
    end
    out = similar(y, Float64); k = 1
    for (vi, wi) in zip(v, w), _ in 1:wi
        out[k] = vi; k += 1
    end
    # Pooling leaves ties, and a flat segment is as bad as a dip — the loading that reads a given
    # fO2 stops being unique. Break them with a negligible strictly-increasing ramp.
    n = length(out)
    if n > 1
        span = max(out[end] - out[1], 1.0)
        for i in 1:n
            out[i] += 1e-6*span*(i - 1)/(n - 1)
        end
    end
    return out
end

# Tabulates solid fO2 against the TOₑₓ share the solid holds, on a SIGNED axis (+ Fe³⁺, - Fe⁰).
# Built in sb24's Fe + O basis: an FeO + O basis double-counts the iron's own oxygen.
function solid_fO2_interpolant(P, T, X, Xox, Ms, TOₑₓ, s_stoich; nr=25, Rmax=0.2, nrm=15, R0=2e-3, data=nothing, verbose=false,
                               mag_cache=nothing)
    # If no solid phase at all (ϕ = 1) return a placeholder and skip MAGEMin call
    Ms <= 0 && return (_ -> 0.0), 0.0, 0.0, true, 0.0, NaN, (_ -> 0.0)

    Xoxv   = Vector{String}(Xox)
    iO     = findfirst(==("O"), Xoxv)
    Xred   = oxidized_Xox(Xoxv, "Fe_O")     # the labels MAGEMin is actually given
    Xdummy = (X=zeros(length(Xoxv)), Xox=copy(Xoxv), mm=get_Xoxmm(Xoxv))
    # wt_in: composition arrives as MASS fractions but R must be imposed in mole space. wt_out: sys_in="wt".
    oxb(r) = oxidize_bulk(X, r, Xdummy; wt_in=true, wt_out=true, frac=true, FeFormat="Fe_O", vector=true)

    # Signed r: Fe³⁺/ΣFe above the FeO reference, -2× the metal fraction below; -2 is singular.
    # asinh-spaced as `solve_sH2O_fO2`: uniform halves put every node above r ≈ Rmax/nr, so raising nr made it WORSE.
    rlist = Rasinhspace(-1.96, Rmax, nr + nrm; R0=R0)   # ascending: reduced → oxidised
    ntot  = length(rlist)
    Xlist = [oxb(r) for r in rlist]
    # Excess O vs the FeO reference, mass fraction at that r: iron holds nFe(1 + r/2) O against nFe,
    # so the excess is just the r/2 part of the O column — no need to re-derive nFe.
    sOwt  = [Xlist[i][iO]*(0.5rlist[i])/(1 + 0.5rlist[i]) for i in 1:ntot]

    # mag_cache (a Dict) keeps the one multi_point result per (P, T, bulk, grid): Ms and TOₑₓ only rescale it
    key = (P, T, Tuple(X), nr, Rmax, nrm, R0)
    out = isnothing(mag_cache) ? nothing : get(mag_cache, key, nothing)
    if isnothing(out)
        flag = isnothing(data)
        flag && (data = Initialize_MAGEMin("sb24", verbose=false))
        out  = multi_point_minimization(10P*ones(ntot), k2c(T)*ones(ntot), data, X=Xlist, Xoxides=Xred, name_solvus=true, sys_in="wt", progressbar=false)
        flag && Finalize_MAGEMin(data)
        isnothing(mag_cache) || (mag_cache[key] = out)
    end

    # Knot vector in GRAMS of excess oxygen: dividing by a signed or near-zero TOₑₓ reversed the
    # knots for a reduced cell or inflated them past convergence.
    grid = sOwt .* Ms
    fO2  = [out[i].fO2 for i in eachindex(out)]

    # Native iron off the assemblage — the only honest source: disproportionation moves no oxygen
    # (2 Fe³⁺ consume the 1 O that 1 Fe⁰ frees), so sOₑₓ cannot see metal made that way.
    isnativeFe(ph) = lowercase(ph) in ("feg", "fea", "fee")
    metal = [sum(o.ph_frac_wt[i] for i in eachindex(o.ph) if isnativeFe(o.ph[i]); init=0.0) for o in out]

    # QFM = fO2 - dQFM is free and (P,T)-only, so identical at every scan point; median guards a bad
    # one. NaN ⇒ `reference_buffers` rebuilds it from Gibbs.
    qfmv = filter(isfinite, [out[i].fO2 - out[i].dQFM for i in eachindex(out)])
    sort!(qfmv)
    qfm  = isempty(qfmv) ? NaN : qfmv[(length(qfmv) + 1) ÷ 2]

    # MAGEMin can fail to converge at some sampled points; keep only what it actually solved.
    good = isfinite.(fO2)
    if !all(good)
        verbose  && (@warn "solid_fO2_interpolant: MAGEMin returned $(count(!, good))/$ntot non-finite fO2 value(s) (Fe³⁺/ΣFe ∈ [$(rlist[1]), $(rlist[end])]) — dropping those grid points.")
        grid, fO2, metal = grid[good], fO2[good], metal[good]
    end

    # Must be monotone or the inversion is ill-posed and the root can vanish; MAGEMin dips it where
    # the assemblage changes. Isotonic regression is the least-squares fit, ordered points untouched.
    monotonic = issorted(fO2)
    if !monotonic
        if verbose
            @warn "solid_fO2_interpolant: MAGEMin fO2 is non-monotonic in the Fe loading at P=$P GPa, T=$T K " *
                "(largest dip $(round(maximum(fO2[i]-minimum(fO2[i:end]) for i in eachindex(fO2)), digits=3)) log units) — " *
                "likely a phase-assemblage change (native Fe coming in, or an Fe³⁺ host). " *
                "Repaired by isotonic regression to keep the solve well-posed."
        end
        fO2 = _isotonic(fO2)
    end

    itp = extrapolate(Interpolations.interpolate((grid,), fO2, Gridded(Linear())), Flat())
    # Metal is monotonically DEcreasing in the loading, so it needs no repair; Flat() is right at
    # both ends (saturated below the table, zero above it).
    mitp = extrapolate(Interpolations.interpolate((grid,), metal, Gridded(Linear())), Flat())

    # ABOVE: continue log-linear in R/(1-R), so it diverges at the stoichiometric limit rather than
    # at the sampling edge. BELOW: flat is correct — Fe⁰ + ½O₂ = FeO buffers at a FINITE fO2.
    g1, f1 = grid[1], fO2[1]
    gN, fN = grid[end], fO2[end]
    slope  = 1/0.1917
    RN     = gN/s_stoich
    upper  = RN < 1 - 1e-9                       # is there any room left above the table at all?
    lnN    = upper ? log10(RN/(1 - RN)) : 0.0

    function sfO2(s)
        s < g1  && return f1
        s <= gN && return itp(s)
        upper   || return fN
        R = min(s/s_stoich, 1 - 1e-12)
        return fN + slope*(log10(R/(1 - R)) - lnN)
    end

    # Full metal saturation (r = -2): read the Fe column, not O which is zero there — the oxygen
    # given back is one per iron, i.e. the iron mass fraction reweighted from mm.Fe to mm.O.
    iFe = findfirst(==("Fe"), Xred)
    return sfO2, grid[1], grid[end], monotonic, -oxb(-2.0)[iFe]*(mm.O/mm.Fe)*Ms, qfm,
           s -> max(0.0, mitp(clamp(s, grid[1], grid[end])))
end

# ---- Section 3 — melt branch: Hirschmann (2022) oxybarometer ----

# log₁₀ fO2 of a melt holding a fraction mOₑₓ of the TOₑₓ budget as Fe³⁺.
function melt_fO2(mOₑₓ, T, TOₑₓ, Mf, molfracXB, nB_per_mass, IDV)
    a=0.1917; b=-1.961; c=4158.1; ΔCₚ=33.25; T₀=1673.15                             # Hirschmann (2022)
    y1=-520.46; y3=494.39; y4=1838.34; y5=2888.48; y8=-1245.09; y9=-1156.86         # regression terms

    cat = (1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 0.0)   # cations per mole of oxide, Xox order

    nOxideTotal = Mf*nB_per_mass           # moles of the melt's own unoxidized oxides
    nOₑₓ_added  = mOₑₓ/mm.O                # moles of excess O added to the melt (mOₑₓ in grams)
    nFeO        = molfracXB[3]*nOxideTotal # moles of Fe, however it is oxidized

    nOxideTotal <= 0 && return oftype(mOₑₓ, 14.0)
    (2nOₑₓ_added >= nFeO) && return oftype(mOₑₓ, 14.0)  # more O than there is Fe to host as Fe₂O₃

    # Single-cation basis: N counts cations, so Fe counts once per Fe whether it sits as FeO or
    # FeO₁.₅ and N does not move as the melt oxidizes. The excess O is NOT in the denominator.
    N      = nOxideTotal*sum(cat[i]*molfracXB[i] for i in eachindex(molfracXB))
    XSiO2  =  molfracXB[1]*nOxideTotal/N
    XMgO   =  molfracXB[2]*nOxideTotal/N
    XCaO   =  molfracXB[4]*nOxideTotal/N
    XAlO15 = 2molfracXB[5]*nOxideTotal/N
    XNaO05 = 2molfracXB[6]*nOxideTotal/N

    ratio = 2nOₑₓ_added/(nFeO - 2nOₑₓ_added)   # 2 FeO + O -> 2 FeO₁.₅, so Fe³⁺/Fe²⁺ = 2nO/(nFeO-2nO)

    # IDV arrives as kJ/mol; the RT ln10 division is done here and nowhere else.
    _T, _ln10 = 1/T, 1/log(10)
    return (log10(ratio) - b - c*_T + (ΔCₚ*_ln10/R)*(1 - T₀*_T - log(T/T₀)) + (IDV/(1e-3R))*_T*_ln10
              - _T*(y1*XSiO2 + y3*XMgO + y4*XCaO + y5*XNaO05 + y8*XSiO2*XAlO15 + y9*XSiO2*XMgO)) / a
end

# ---- Section 4 — melt CO₂ branch: EMOG(D) buffer (Stagno & Frost, 2010-style) ----

# XCO₂ of a melt at graphite/diamond saturation for a given fO2 (inverse of the EMOG(D)-style melt buffer).
invert_carbon_buffer_fO2(f, P, T) = 10.0^(f - 5.44 + 21380/T - 0.078*(1e4P - 1)/T)

# ---- Section 4b — C-O-H fluid: Holland & Powell (1998) Kᵢ(T), Zhang & Duan (2009) fugacities ----

# Holland & Powell (1998): G0/H [J/mol] at 298.15 K, S [J/mol/K], Cp = c1 + c2·T + c3·T⁻² + c5·T⁻⁰·⁵.
# Carbon polymorphs add Murnaghan volume terms (V0, b1, b5, b6=K₀, b7=dK/dT, b8=K′) for their μ°(P).
const HP98 = (
    H2O = (G0=-228542.3,   H=-241810.0, S=188.8,  c1=40.1,  c2=8.656e-3,  c3=487500.0,  c5=-251.2),
    CO2 = (G0=-394332.9,   H=-393510.0, S=213.7,  c1=87.8,  c2=-2.644e-3, c3=706400.0,  c5=-998.9),
    CO  = (G0=-137163.7,   H=-110530.0, S=197.67, c1=45.7,  c2=-9.7e-5,   c3=662700.0,  c5=-414.7),
    CH4 = (G0=-50695.63,   H=-74810.0,  S=186.26, c1=150.1, c2=2.062e-3,  c3=3427700.0, c5=-2650.4),
    H2  = (G0=-0.0045494,  H=0.0,       S=130.7,  c1=23.3,  c2=4.627e-3,  c3=0.0,       c5=76.3),
    O2  = (G0=0.0,         H=0.0,       S=205.2,  c1=48.3,  c2=-6.91e-4,  c3=499200.0,  c5=-420.7),
    gph = (G0=-32.79654,   H=0.0,       S=5.85,   c1=51.0,  c2=-4.428e-3, c3=488600.0,  c5=-805.5,
           V0=0.5298, b1=0.484e-4, b5=-0.484e-3, b6=390000.0, b7=-58.5, b8=4.0),
    diam = (G0=3125.636,  H=2100.0,    S=2.3,    c1=24.3,  c2=6.272e-3,  c3=-377400.0, c5=-273.4,
           V0=0.3417, b1=0.165e-4, b5=-0.165e-3, b6=5800000.0, b7=-870.0, b8=4.0),
)

# Gibbs energy thermal contribution up to temperature T
function hp98_gibbs(T, d)
    Tr = 298.15
    ΔH = d.c1*(T-Tr) + 0.5d.c2*(T^2-Tr^2) - d.c3*(1/T-1/Tr) + 2d.c5*(√T-√Tr)
    ΔS = d.c1*log(T/Tr) + d.c2*(T-Tr) - 0.5d.c3*(1/T^2-1/Tr^2) - 2d.c5*(1/√T-1/√Tr)
    return d.G0 - (T-Tr)*d.S + ΔH - T*ΔS
end

# ΔG of compressing a solid from 1 bar to Pbar (Holland & Powell modified-Tait/Murnaghan volume EOS).
function murnaghan_dG(T, Pbar, d)
    Tr = 298.15
    V0T = d.V0*exp(d.b1*(T-Tr) + 2d.b5*(√T-√Tr))
    K   = d.b6 + d.b7*(T-Tr)
    Kp  = d.b8
    return V0T*K/(Kp-1) * ((1 + Kp*(Pbar-1)/K)^(1-1/Kp) - 1)
end

# Equilibrium constants for  C+O₂=CO₂,  C+½O₂=CO,  H₂+½O₂=H₂O,  C+2H₂=CH₄.
function reaction_Ks(T, Pbar)
    Ggph = hp98_gibbs(T, HP98.gph)  + murnaghan_dG(T, Pbar, HP98.gph)
    Gdia = hp98_gibbs(T, HP98.diam) + murnaghan_dG(T, Pbar, HP98.diam)
    μC   = min(Ggph, Gdia)
    μO2  = hp98_gibbs(T, HP98.O2);  μCO2 = hp98_gibbs(T, HP98.CO2)
    μCO  = hp98_gibbs(T, HP98.CO);  μH2O = hp98_gibbs(T, HP98.H2O)
    μH2  = hp98_gibbs(T, HP98.H2);  μCH4 = hp98_gibbs(T, HP98.CH4)
    K1 = exp(-(μCO2 - μC - μO2)/(R*T))
    K2 = exp(-(μCO  - μC - 0.5μO2)/(R*T))
    K3 = exp(-(μH2O - μH2 - 0.5μO2)/(R*T))
    K4 = exp(-(μCH4 - μC - 2μH2)/(R*T))
    return K1, K2, K3, K4
end

# --- Zhang & Duan (2009) equation of state. Species order: CH4, H2O, CO2, H2, CO, O2, C2H6. ---
const ZD_c = (2.95177298930e-2, -6.33756452413e3, -2.75265428882e5, 1.29128089283e-3, -1.45797416153e2,
              7.65938947237e4, 2.58661493537e-6, 0.52126532146e0, -1.39839523753e2, -2.36335007175e-8,
              5.35026383543e-3, -0.27110649951e0, 2.50387836486e4, 0.73226726041e0, 1.54833359970e-2)
const ZD_eps = (154.0, 510.0, 235.0, 31.2, 105.6, 124.5, 246.1)   # Lennard-Jones ε/k_B [K]
const ZD_sig = (3.691, 2.88, 3.79, 2.93, 3.66, 3.36, 4.35)        # Lennard-Jones σ [Å]
const ZD_Rgas = R/100   # L·bar/(mol·K); same constant as R, since 1 J = 0.01 L·bar

# Binary interaction parameters (Zhang & Duan 2009, Table 7); unity elsewhere.
function zd_default_K1K2()
    K1 = ones(7,7); K2 = ones(7,7)
    K1[2,3] = K1[3,2] = 0.85; K1[1,2] = K1[2,1] = 0.8   # H2O-CO2, CH4-H2O
    K2[2,3] = K2[3,2] = 1.02
    return K1, K2
end
const ZD_K1, ZD_K2 = zd_default_K1K2()

# Residual of the reduced-variable EOS; its root in V is the reduced molar volume.
function zd_EOS(V, Tm, targetPm)
    V <= 0 && return Inf
    T2, T3 = Tm^2, Tm^3
    V2, V4 = V^2, V^4
    c = ZD_c
    Z = 1.0
    Z += (c[1] + c[2]/T2 + c[3]/T3) / V
    Z += (c[4] + c[5]/T2 + c[6]/T3) / V2
    Z += (c[7] + c[8]/T2 + c[9]/T3) / V4
    Z += (c[10] + c[11]/T2 + c[12]/T3) / (V4*V)
    Z += c[13]/T3/V2 * (c[14] + c[15]/V2) * exp(-c[15]/V2)
    return ZD_Rgas*Tm*Z/V - targetPm
end

# Root acceptance is RELATIVE to the residual's own scale Pm (F has units of Pm, ~1e5 here):
# an absolute 1e-6 test rejected converged roots on round-off alone.
const ZD_FTOL = 1e-10

# The EOS has several real roots; brackets are tried widest-first for the fluid-like one. Inside a bracket:
# Newton in density d = 1/V with the analytic dP/dd, halving the bracket when a step leaves it (TOex_numerics.md 2.2).
function zd_findroot(Tm, Pm)
    T2, T3 = Tm^2, Tm^3; c = ZD_c
    g1 = c[1] + c[2]/T2 + c[3]/T3; g2 = c[4] + c[5]/T2 + c[6]/T3; g3 = c[7] + c[8]/T2 + c[9]/T3; g4 = c[10] + c[11]/T2 + c[12]/T3; g5 = c[13]/T3
    function FdF(d)
        d2 = d*d; e = exp(-c[15]*d2); q = c[14] + c[15]*d2
        F  = ZD_Rgas*Tm*d*(1 + g1*d + g2*d2 + g3*d2*d2 + g4*d2*d2*d + g5*d2*q*e) - Pm
        dF = ZD_Rgas*Tm*(1 + 2g1*d + 3g2*d2 + 5g3*d2*d2 + 6g4*d2*d2*d + g5*e*d2*(3q + 2c[15]*d2*(1 - q)))
        return F, dF
    end
    for (lo, hi) in ((1.0,1000.0), (0.01,10.0), (0.001,1.0))
        dl, dh = 1/hi, 1/lo; Fl = FdF(dl)[1]
        Fl*FdF(dh)[1] > 0 && continue
        d = sqrt(dl*dh)
        for _ in 1:200
            F, dF = FdF(d)
            sign(F) == sign(Fl) ? (dl = d; Fl = F) : (dh = d)
            dn = d - F/dF
            (isfinite(dn) && dl < dn < dh) || (dn = 0.5*(dl + dh))
            conv = abs(dn - d) <= 1e-13*dn; d = dn
            conv && break
        end
        abs(zd_EOS(1/d, Tm, Pm)) <= ZD_FTOL*max(abs(Pm), 1.0) && return 1/d
    end
    return NaN
end

# Molar volume [cc/mol] of a fluid of composition vx (contains Xi of COH species).
function zd_volume(vx, Pbar, TK)
    ε = 0.0; σ = 0.0
    for i in 1:7, j in 1:7
        ε += vx[i]*vx[j]*ZD_K1[i,j]*√(ZD_eps[i]*ZD_eps[j])
        σ += vx[i]*vx[j]*ZD_K2[i,j]*(ZD_sig[i]+ZD_sig[j])/2
    end
    (ε <= 0 || σ <= 0) && return NaN
    Pm = 3.06258806161*σ^3*Pbar/ε   # ε(CH4)/σ(CH4)³ at full precision, not the paper's rounded 3.0636
    Tm = 154.0*TK/ε
    Vm = zd_findroot(Tm, Pm)
    isnan(Vm) && return NaN
    return 1000.0*Vm*(σ/3.691)^3
end

# Returns (volume [cc/mol], fugacities [bar]) for all 7 species at composition vx.
function zd_fugacity(vx, Pbar, TK)
    vol = zd_volume(vx, Pbar, TK) / 1000.0   # L/mol
    isnan(vol) && return NaN, fill(NaN,7)
    Z = Pbar*vol/ZD_Rgas/TK
    εi, σi = zeros(7), zeros(7)
    ε, σ = 0.0, 0.0
    for i in 1:7, j in 1:7
        εi[i] += ZD_K1[i,j]*vx[j]*√(ZD_eps[i]*ZD_eps[j])
        σi[i] += ZD_K2[i,j]*vx[j]*(ZD_sig[i]+ZD_sig[j])/2
        ε += vx[i]*vx[j]*ZD_K1[i,j]*√(ZD_eps[i]*ZD_eps[j])
        σ += vx[i]*vx[j]*ZD_K2[i,j]*(ZD_sig[i]+ZD_sig[j])/2
    end
    Tm = 154.0*TK/ε
    Vm = vol*(3.691/σ)^3
    c = ZD_c
    S1  = (c[1] + c[2]/Tm^2 + c[3]/Tm^3)/Vm
    S1 += (c[4] + c[5]/Tm^2 + c[6]/Tm^3)/(2Vm^2)
    S1 += (c[7] + c[8]/Tm^2 + c[9]/Tm^3)/(4Vm^4)
    S1 += (c[10] + c[11]/Tm^2 + c[12]/Tm^3)/(5Vm^5)
    S1 += c[13]/(2c[15]*Tm^3) * (c[14] + 1 - (c[14]+1+c[15]/Vm^2)*exp(-c[15]/Vm^2))
    S2  = (2c[2]/Tm^2 + 3c[3]/Tm^3)/Vm
    S2 += (2c[5]/Tm^2 + 3c[6]/Tm^3)/(2Vm^2)
    S2 += (2c[8]/Tm^2 + 3c[9]/Tm^3)/(4Vm^4)
    S2 += (2c[11]/Tm^2 + 3c[12]/Tm^3)/(5Vm^5)
    S2 += 3c[13]/(2c[15]*Tm^3) * (c[14] + 1 - (c[14]+1-c[15]/Vm^2)*exp(-c[15]/Vm^2))
    fug = zeros(7)
    for i in 1:7
        lnγ = Z - 1 - log(Z) + S1 - 2S2*(1 - εi[i]/ε) + 6(1-Z)*(1 - σi[i]/σ)
        fug[i] = Pbar*exp(lnγ)*vx[i]
    end
    return vol*1000.0, fug
end

# ---- Section 4b2 — low-temperature C-O-H fluid: hybrid MRK + Pitzer & Sterner EoS, blended into ZD09 ----
# ZD09 has no liquid-water root below ~600 K (it returns ~3 cc/mol and φ ~ 1e-238). Below FL_T_LO the fluid is
# the hybrid of Connolly (Perple_X): MRK for all species and their mixing (Holloway 1977, de Santis et al. 1974),
# with the pure-H2O and pure-CO2 parts replaced by Pitzer & Sterner (1994). ln φ is blended into ZD09 up to FL_T_HI.
const FL_T_LO, FL_T_HI = 673.0, 773.0
const RK_R = 83.1441                 # cm³·bar/(mol·K), the value the MRK constants were fitted with

# MRK a [bar·cm⁶·K^½/mol²] and b [cm³/mol] in ZD09's species order (CH4, H2O, CO2, H2, CO, O2, C2H6).
# H2O, CO2 and H2 carry the temperature-dependent dispersion terms of Holloway (1977) / Connolly's H2 fit.
function mrk_ab(T)
    aw = T > 300.0 ? 1.452535403e7 + T*(306893.3587 + T*(-307.9995871 + T*(0.09226256008 - 2.930106337e-6T))) :
                     127354240.0                     # 298-300 K: the liquid-water pair (b = 16) of the same scheme
    bw = T > 300.0 ? 14.6 : 16.0
    ah = 3.91950132949994654e7 + T*(-8.81231157499978144e4 + T*(89.0185987380923081 - 2.86881183333320412e-2T))
    a = (32.154811e6, aw, 92935540.0 + T*(-82130.73 + 21.29T), ah, 16.98e6, 17.4026e6, 98774720.4)
    b = (29.681, bw, 29.7, 12.81508162, 27.38, 22.07, 45.139)
    return map(x -> max(x, 1.0), a), b
end
# H2O–CO2 cross term with the association correction of de Santis et al. (1974); a_mix gets x_w·x_c·mrk_awc
mrk_awc(T) = exp(-11.218 + (6032.0 + (-2782000.0 + 4.708e8/T)/T)/T)*6912.824964*T^2*√T + 79267647.0

# Real roots of V³ + c1·V² + c2·V + c3 = 0
function cubic_roots(c1, c2, c3)
    q = (c1^2 - 3c2)/9; r = (2c1^3 - 9c1*c2 + 27c3)/54
    if r^2 < q^3
        θ = acos(clamp(r/sqrt(q^3), -1.0, 1.0)); s = -2sqrt(q)
        return sort([s*cos(θ/3) - c1/3, s*cos((θ + 2π)/3) - c1/3, s*cos((θ - 2π)/3) - c1/3])
    end
    A = -sign(r)*cbrt(abs(r) + sqrt(r^2 - q^3)); B = A == 0 ? 0.0 : q/A
    return [A + B - c1/3]
end

# MRK molar volume: of three roots, the one of lower Gibbs energy (liquid vs vapour), else the positive one
function mrk_volume(a, b, Pbar, T)
    rt, st = RK_R*T, √T
    c1 = -rt/Pbar
    v = cubic_roots(c1, c1*b + a/(st*Pbar) - b^2, -a*b/(st*Pbar))
    v = filter(x -> x > b, v)
    isempty(v) && return NaN
    length(v) == 1 && return v[1]
    vmin, vmax = first(v), last(v)
    pdv = Pbar*(vmax - vmin) - rt*log((vmax - b)/(vmin - b)) - a/(b*st)*log((vmax + b)/(vmin + b)*vmin/vmax)
    return pdv > 0 ? vmin : vmax
end

# MRK ln φᵢ of every species in the mixture vx (also at infinite dilution) and the mixture volume
function mrk_lnφ(vx, Pbar, T)
    a, b = mrk_ab(T); awc = mrk_awc(T)
    aij(i, j) = (i, j) in ((2, 3), (3, 2)) ? awc/2 : √(a[i]*a[j])
    am = 0.0; bm = 0.0
    for i in 1:7
        bm += vx[i]*b[i]
        for j in 1:7
            am += vx[i]*vx[j]*aij(i, j)
        end
    end
    v = mrk_volume(am, bm, Pbar, T)
    isfinite(v) || return NaN, fill(NaN, 7)
    rt, st = RK_R*T, √T
    d2 = log((v + bm)/v)/(rt*st*bm)
    d3 = am*d2/bm - am/((v + bm)*rt*st*bm) + 1/(v - bm)
    d6 = log(rt/(v - bm)) - log(Pbar)
    lnφ = [b[i]*d3 - 2sum(vx[j]*aij(i, j) for j in 1:7)*d2 + d6 for i in 1:7]
    return v, lnφ
end

# Pitzer & Sterner (1994) pure H2O (k = 1) or CO2 (k = 2): P(V) and ln f. Every root of P(V) = P on a log
# grid is refined and the one of lowest Gibbs energy (lowest ln f) is kept, so liquid water cannot fall
# onto the vapour-like branch. Returns (V [cm³/mol], ln φ).
function ps94_coeffs(T, k)
    t2 = T^2
    k == 1 && return (0.24657688e6/T + 0.51359951e2,
                      0.58638965/T - 0.28646939e-2 + 0.31375577e-4T,
                      -0.62783840e1/T + 0.14791599e-1 + T*(0.35779579e-3 + 0.15432925e-7T),
                      -0.42719875 - 0.16325155e-4T,
                      0.56654978e4/T - 0.16580167e2 + 0.76560762e-1T,
                      0.10917883,
                      ((0.38878656e13/t2 - 0.13494878e9)/T + 0.30916564e6)/T + 0.75591105e1,
                      -0.65537898e5/T + 0.18810675e3,
                      ((-0.14182435e14/t2 + 0.18165390e9)/T - 0.19769068e6)/T - 0.23530318e2,
                      0.92093375e5/T + 0.12246777e3)
    return (0.18261340e7/T + 0.79224365e2,
            0.66560660e-4 + 0.57152798e-5T + 0.30222363e-9t2,
            0.59957845e-2 + 0.71669631e-4T + 0.62416103e-8t2,
            -0.13270279e1/T - 0.15210731 + 0.53654244e-3T - 0.71115142e-7t2,
            0.12456776/T + 0.49045367e1 + 0.98220560e-2T + 0.55962121e-5t2,
            0.75522299,
            ((-0.39344644e12/t2 + 0.90918237e8)/T + 0.42776716e6)/T - 0.22347856e2,
            0.40282608e3/T + 0.11971627e3,
            (0.22995650e8/T - 0.78971817e5)/T - 0.63376456e2,
            0.95029765e5/T + 0.18038071e2)
end
function ps94_pure(Pbar, T, k)
    c1, c2, c3, c4, c5, c6, c7, c8, c9, c0 = ps94_coeffs(T, k)
    rt = 83.14*T                                       # the paper's R in cm³·bar/(mol·K)
    a1(v) = c2 + (c3 + (c4 + (c5 + c6/v)/v)/v)/v
    function pv(v)                                     # P(V)/RT
        e1, e2 = c7*exp(-c8/v), c9*exp(-c0/v)
        return (1 + (c1 + e1 + e2)/v - (c3 + (2c4 + (3c5 + 4c6/v)/v)/v)/v/a1(v)^2)/v
    end
    lnf(v) = c1/v + 1/a1(v) - 1/c2 - (c7*exp(-c8/v) - c7)/c8 - (c9*exp(-c0/v) - c9)/c0 + log(rt/v) + Pbar*v/rt - 1
    F(v) = pv(v) - Pbar/rt
    lv = range(log(5.0), log(max(1e6, 50rt/Pbar)), length=400)
    best, vb = Inf, NaN
    Fl = F(exp(lv[1]))
    for k2 in 2:length(lv)
        Fh = F(exp(lv[k2]))
        if isfinite(Fl) && isfinite(Fh) && Fl*Fh <= 0
            lo, hi, flo = lv[k2-1], lv[k2], Fl        # bisection in ln V: P(V) is not monotone across branches
            for _ in 1:80
                m = 0.5(lo + hi); fm = F(exp(m))
                (fm*flo > 0) ? (lo = m; flo = fm) : (hi = m)
            end
            v = exp(0.5(lo + hi)); g = lnf(v)
            (isfinite(g) && g < best) && (best = g; vb = v)
        end
        Fl = Fh
    end
    return vb, best - log(Pbar)
end

# Hybrid fugacities: MRK mixture, with H2O and CO2 shifted by ln φ(PS94, pure) - ln φ(MRK, pure).
# The volume takes the same per-species shift, weighted by mole fraction.
function hyb_fugacity(vx, Pbar, TK)
    v, lnφ = mrk_lnφ(vx, Pbar, TK)
    isfinite(v) || return NaN, fill(NaN, 7)
    a, b = mrk_ab(TK)
    for (i, k) in ((2, 1), (3, 2))
        vps, lps = ps94_pure(Pbar, TK, k)
        isfinite(lps) || continue
        e = zeros(7); e[i] = 1.0
        vm, lm = mrk_lnφ(e, Pbar, TK)
        lnφ[i] += lps - lm[i]
        v += vx[i]*(vps - vm)
    end
    fug = [vx[i] > 0 ? vx[i]*Pbar*exp(lnφ[i]) : 0.0 for i in 1:7]
    return v, fug
end

# The fluid EoS used by every TOex speciation: hybrid below FL_T_LO, ZD09 above FL_T_HI, ln φ and ln V
# blended with a smoothstep weight in between (both are finite there; either alone if the other fails).
function fluid_fugacity(vx, Pbar, TK)
    TK >= FL_T_HI && return zd_fugacity(vx, Pbar, TK)
    vh, fh = hyb_fugacity(vx, Pbar, TK)
    TK <= FL_T_LO && return vh, fh
    vz, fz = zd_fugacity(vx, Pbar, TK)
    okz = isfinite(vz) && all(i -> vx[i] <= 0 || (isfinite(fz[i]) && fz[i] > 0), 1:7)
    okh = isfinite(vh) && all(i -> vx[i] <= 0 || (isfinite(fh[i]) && fh[i] > 0), 1:7)
    okz || return vh, fh
    okh || return vz, fz
    s = (TK - FL_T_LO)/(FL_T_HI - FL_T_LO); w = s^2*(3 - 2s)
    fug = [vx[i] > 0 ? exp((1 - w)*log(fh[i]) + w*log(fz[i])) : 0.0 for i in 1:7]
    return exp((1 - w)*log(vh) + w*log(vz)), fug
end

# ---- Section 3b — molecular H₂ dissolved in the melt (Hirschmann et al. 2012) ----

# ln K = ln(X_H2/f_H2) at T₀ = 1400 °C, P₀ = 1 bar. Their T range is too narrow to calibrate ΔH, so
# it defaults to 0 and lnK carries no explicit T dependence — only the ΔV(P-P₀)/RT term.
const H2_MELT = (basalt     = (lnK0 = -11.403, ΔV = 10.6, IP = 45.6),   # Table 2 / Fig. 4 regression
                 andesite   = (lnK0 = -10.591, ΔV = 11.3, IP = 48.3),
                 peridotite = (lnK0 = -12.5,   ΔV = 10.6, IP = 42.1))   # extrapolated on ionic porosity
const H2_T0, H2_P0 = 1673.15, 1.0   # K, bar

# Henry's-law constant for H₂(fluid) ⇌ H₂(melt). ΔV in cm³/mol against P in bar carries the 0.1.
function h2_melt_K(T, Pbar; comp=:basalt, ΔH=0.0)
    d = H2_MELT[comp]
    return exp(d.lnK0 + (ΔH/R)*(1/H2_T0 - 1/T) - 0.1*d.ΔV*(Pbar - H2_P0)/(R*T))
end

# Sun & Yao (2026) Eq. (7) with their Table 1: total dissolved H₂O as molecular + hydroxyl, both
# returned in wt% H₂O. P and PH2O in bar, T in K, X oxide mole fractions on a volatile-free basis.
const SY_H2O = (0.807805, -0.700336, 2.949559e-6, -74.546597, 7146.275452,
                0.449387, -0.360194, -0.308259,   2.984703,  -1.693710)

# Takes the PARTIAL PRESSURE X_H2O·P as calibrated, never a ZD09 fugacity. The H₂O half is fitted
# only to 2 GPa and its a₃P²/T term diverges above, so P is clamped there at fixed X_H2O (math S.10).
const SY_PMAX = 2.0e4    # bar

# The one melt-water law: the solver passes its own a_H2O, everyone else 1 (pure-water saturation).
# :katz is the legacy 12·P^0.6 + P; with :sunyao activity enters via P_H2O, so `n` survives for the solid only.
function melt_h2o_wt(P, T, aH2O; law=:sunyao, clamp_P=true, n_katz=1.0)
    a = clamp(aH2O, 0.0, 1.0)
    if law === :katz
        # Katz (2003) with the Rozel (2018) cutoff. No speciation, no T dependence, so the
        # activity has to be applied as the a_H2O^n post-multiply rather than through P_H2O.
        S = min(12*P^0.6 + P, 12*10.0^0.6 + 10.0)*a^n_katz
        return (S, 0.0)          # reported wholly as "molecular"; the split is not resolved
    end
    Pbar = 1e4*P
    return h2o_melt_wt(a*Pbar, Pbar, T, MELT_XNA2O, MELT_XAL2O3; clamp_P=clamp_P)
end

# The melt is held at a fixed basaltic XB, so its Na2O and Al2O3 are constants.
const MELT_XNA2O, MELT_XAL2O3 = let (; XB, mmXox) = endmember_composition(1.0)
    m = (XB./mmXox)./sum(XB./mmXox); (m[6], m[5])
end

function h2o_melt_wt(PH2O, Pbar, T, XNa2O, XAl2O3; clamp_P=true)
    (PH2O <= 0 || !isfinite(PH2O)) && return (0.0, 0.0)
    Pc, Ph = Pbar, PH2O
    if clamp_P && Pbar > SY_PMAX
        Pc = SY_PMAX;  Ph = PH2O*(SY_PMAX/Pbar)     # hold P at the ceiling, keep the mole fraction
    end
    a = SY_H2O
    lnSm = a[1]*log(Ph) + a[2]*log(T) + a[3]*Pc^2/T + (a[4]*sqrt(Pc) + a[5])*XNa2O/T
    lnSo = a[6]*log(Ph) + a[7]*log(T) + (a[8] + a[9]*XAl2O3 + a[10]*XNa2O)*Pc/T
    return (exp(lnSm), exp(lnSo))          # wt% as molecular H₂O, wt% as OH⁻ counted as H₂O
end

# ---- Section 4b3 — EMOG(D) buffer: graphite/diamond against carbonate ----

# EMOG(D), the graphite/diamond ⇄ carbonate boundary: 2 MgCO₃ + Mg₂Si₂O₆ = 2 Mg₂SiO₄ + 2C + 2O₂,
# so log₁₀fO2 = -ΔG/(2·R·T·ln10) at unit activities. Free carbon below it, carbonate above.
const EMOG_MAG = (G0=-1027504.0, S=65.1,  c1=186.4, c2=-0.3772e-2, c3=0.0,       c5=-1886.2,
                  V0=2.803, b1=0.648e-4, b5=-0.648e-3, b6=1460000.0, b7=-219.0, b8=4.0)
const EMOG_EN  = (G0=-2915760.0, S=132.5, c1=356.2, c2=-0.299e-2,  c3=-596900.0, c5=-3185.3,
                  V0=6.262, b1=0.505e-4, b5=-0.505e-3, b6=1070000.0, b7=-160.5, b8=4.0)
const EMOG_FO  = (G0=-2053138.0, S=95.1,  c1=233.3, c2=0.1494e-2,  c3=-603800.0, c5=-1869.7,
                  V0=4.366, b1=0.613e-4, b5=-0.613e-3, b6=1250000.0, b7=-187.5, b8=4.0)
function logfO2_EMOG(T, Pbar)
    G(d) = hp98_gibbs(T, d) + murnaghan_dG(T, Pbar, d)
    μC = min(G(HP98.gph), G(HP98.diam))   # whichever carbon polymorph is stable here
    ΔG = 2G(EMOG_FO) + 2μC + 2hp98_gibbs(T, HP98.O2) - 2G(EMOG_MAG) - G(EMOG_EN)
    return -ΔG/(2*R*T*log(10))
end

# ---- Section 4c — reference buffers: FMQ and IW, both out of sb24 ----

# Both buffers out of sb24 so an offset shares its dataset with the solved fO2 (an empirical curve put
# metal-bearing rock above IW). IW built here because MAGEMin's `iw` phase segfaults; caveats: TOex_math.md S.13.

const T0_SLB = 300.0        # SLB reference temperature; F₀ and the thermal terms are relative to it

# sb24 endmembers at full precision (stx24_data.json): F0 J/mol | n | V0 J/bar/mol | K0 bar | Kp | θ0 K | γ0 | q0 | fe.
# `wu` is Fe₄O₄ (SLB B1 phases sit on four oxygens); `mod` = Landau (SD,TC,VD) or Hillert-Jarl magnetic (TC,β,pf,dE,dS).
const SLB_EM = (
    wu   = (F0=-962050.35, n=8.0, V0=4.9024, K0=1607000.0, Kp=4.0, θ0=427.00102, γ0=1.45033, q0=1.5487, fe=0,
            mod=(SD=53.5254, TC=191.0, VD=0.0)),
    wuls = (F0=-609335.18, n=8.0, V0=4.33997, K0=1997000.0, Kp=4.0, θ0=492.99392, γ0=1.45033, q0=1.5487, fe=0,
            mod=nothing),
    fea  = (F0=11.82, n=1.0, V0=0.7093, K0=1634211.4, Kp=6.01732, θ0=398.00486, γ0=1.66729, q0=0.90658, fe=1,
            mod=(TC=1043.01, β=2.1199287186834, pf=0.4, dE=8932.7399, dS=9.46028)),
    feg  = (F0=4675.83, n=1.0, V0=0.698835, K0=1650000.0, Kp=5.97, θ0=285.19315, γ0=1.84924, q0=1.02176, fe=3,
            mod=nothing),
    fee  = (F0=3428.48, n=1.0, V0=0.67651, K0=1650000.0, Kp=4.97, θ0=379.74281, γ0=2.13092, q0=1.02253, fe=2,
            mod=nothing),
    qtz  = (F0=-858043.33, n=3.0, V0=2.242105, K0=614259.6, Kp=19.78014, θ0=883.46813, γ0=-0.03957, q0=1.0, fe=0,
            mod=(SD=5.76, TC=847.0, VD=0.135936)),
    fa   = (F0=-1371708.34, n=7.0, V0=4.629, K0=1364855.8, Kp=4.88157, θ0=618.96116, γ0=1.08388, q0=2.88055, fe=0,
            mod=(SD=26.7627, TC=65.0, VD=0.0)),
    mag  = (F0=-1015909.99, n=7.0, V0=4.4528, K0=1838766.4, Kp=5.22819, θ0=529.46966, γ0=1.35821, q0=0.41872, fe=0,
            mod=(TC=845.5, β=178.9816942215617, pf=0.4, dE=33048.0797, dS=43.1758)),
)

# SLB property modifier for one endmember, J/mol. Landau is taken relative to the fully ordered
# (Q = 1) state; the magnetic branch is Hillert & Jarl (1978) as SLB adopt it (A5-A9, Table C3).
function slb_modifier(d, Pbar, T)
    m = d.mod
    m === nothing && return 0.0
    G = 0.0

    if haskey(m, :SD)
        Tc = m.TC + m.VD*Pbar/m.SD
        G += -m.SD*((T - Tc) + m.TC/3)
        # below Tc the order parameter is finite; Q² is capped at 4 to stop runaway stabilisation
        if T < Tc
            Q2 = min(sqrt((Tc - T)/m.TC), 4.0)
            G += m.SD*((T - Tc)*Q2 + m.TC*Q2^3/3)
        end
    end

    if haskey(m, :β)
        D = 518/1125 + (11692/15975)*(1/m.pf - 1)
        τ = T/m.TC
        f = τ < 1 ? 1 - (79/(140*m.pf*τ) + (474/497)*(1/m.pf - 1)*(τ^3/6 + τ^9/135 + τ^15/600))/D :
                       -(τ^-5/10 + τ^-15/315 + τ^-25/1500)/D
        G += R*T*log(m.β + 1)*f
    end

    haskey(m, :dE) && (G += m.dE - T*m.dS)
    return G
end

# Metal electronic term ℱ_el = -½·β_el0·(V/V₀)^γ_el·T² (SLB 2024 A2-A3); fe = 1 α, 2 ε, 3 γ.
# Worth ~0.45 log units on IW at 1700 K. The only values the paper omits — from its repository.
const SLB_FE_EL = ((0.00388, 1.47960), (0.00411, 1.69270), (0.00375, 1.47960))

# Debye free-energy integral, HeFESTo/Perple_X form: plg(x) = -π⁴/45 + Σ e^{-ix}(x²/i² + 2x/i³ + 2/i⁴),
# so ℱ_th = 9nRT·plg(θ/T)/(θ/T)³ (the zero-point term cancels against the T₀ reference).
function slb_plg(t)
    p0 = exp(-t); p1 = 1.0; p2 = t*t; p3 = 2t
    s  = -2.1646464674222763831            # -π⁴/45
    for i in 1:100000
        p1  *= p0
        dinc = (p2 + (p3 + 2.0/i)/i)*p1/(i*i)
        s   += dinc
        abs(dinc/(1 + abs(s))) < 1e-15 && break
    end
    return s
end

# Gibbs energy of one sb24 endmember at (Pbar, T) in J/mol — SLB Mie-Grüneisen-Debye on a BM3 cold
# curve, plus metal-only electronic and magnetic terms. V is found by Newton on P_EOS(V,T) = Pbar.
function slb_gibbs(d, Pbar, T)
    nr9    = 9d.n*R
    c1     = 9d.V0*d.K0                    # = 9K₀V₀, the finite-strain energy scale
    c2     = d.Kp/2 - 2
    c3     = 3c1*c2
    aii    = 6d.γ0
    aiikk2 = 0.5*(-12d.γ0 + 36d.γ0^2 - 18d.q0*d.γ0)
    aii2   = 0.5aii
    nr9t0  = nr9*T0_SLB;  nr9t = nr9*T
    t1     = d.θ0/T;      t2   = T/T0_SLB
    delt2  = T*T - T0_SLB*T0_SLB
    b1, b2 = d.fe > 0 ? SLB_FE_EL[d.fe] : (0.0, 0.0)

    # Starting volume: the same expression with θ frozen at θ₀, which inverts in closed form.
    tht   = t1;  tht0 = tht*t2
    dfth  = nr9t *d.γ0/d.V0*(3slb_plg(tht )/tht^3  - log(1 - exp(-tht )))
    dfth0 = nr9t0*d.γ0/d.V0*(3slb_plg(tht0)/tht0^3 - log(1 - exp(-tht0)))
    root  = d.K0*((2 + 2d.Kp)*(Pbar + dfth - dfth0) + d.K0)
    v = root > 0 ? ((2 + d.Kp) - sqrt(root)/d.K0)*d.V0/(1 + d.Kp) : d.V0
    (v < d.V0/10 || v > 10d.V0) && (v = d.V0)

    for _ in 1:2048
        v23  = (d.V0/v)^(2/3)
        f    = 0.5v23 - 0.5                # Eulerian finite strain
        df   = -v23/(3v)
        d2f  = (5/9)*v23/(v*v)
        dfc  = (c3*f + c1)*f*df
        d2fc = (2c3*f + c1)*df*df + (c3*f + c1)*f*d2f
        fel = dfel = d2fel = 0.0
        if d.fe > 0
            fel   = -b1/2*(v/d.V0)^b2*delt2
            dfel  = fel*b2/v
            d2fel = dfel*(b2 - 1)/v
        end
        z = 1 + (aii + aiikk2*f)*f         # (θ/θ₀)²
        (z < 0 || v/d.V0 > 100 || v/d.V0 < 0.01) && break
        root = sqrt(z)
        tht  = t1*root;  tht0 = tht*t2
        a2f  = aii2 + aiikk2*f
        da   = a2f/root
        dtht  = t1*da*df
        d2tht = t1*((aiikk2/root - a2f*a2f/z^1.5)*df*df + da*d2f)
        dtht0 = dtht*t2;  d2tht0 = d2tht*t2
        fpoly  = 3slb_plg(tht )/tht^3
        fpoly0 = 3slb_plg(tht0)/tht0^3
        etht   = exp(-tht);   letht  = log(1 - etht)
        etht0  = exp(-tht0);  letht0 = log(1 - etht0)
        dfth   = (letht - fpoly)*nr9t*dtht/tht
        d2fth  = ((4dtht*dtht/tht - d2tht)*(fpoly - letht) + dtht*dtht*etht/(1 - etht))*nr9t/tht
        dfth0  = (letht0 - fpoly0)*nr9t0*dtht0/tht0
        d2fth0 = ((4dtht0*dtht0/tht0 - d2tht0)*(fpoly0 - letht0) +
                  dtht0*dtht0*etht0/(1 - etht0))*nr9t0/tht0
        f1  = -dfc - dfth + dfth0 - dfel - Pbar     # P_EOS - P, as a function of V
        df1 = -d2fc - d2fth + d2fth0 - d2fel
        dv  = f1/df1
        (v - dv < 0) && (dv = v/2)
        v  -= dv
        abs(dv/(1 + v)) < 1e-12 && break
    end

    f    = 0.5*(d.V0/v)^(2/3) - 0.5
    z    = 1 + (aii + aiikk2*f)*f
    root = sqrt(z);  tht = t1*root;  tht0 = tht*t2
    fel  = d.fe > 0 ? -b1/2*(v/d.V0)^b2*delt2 : 0.0
    G    = d.F0 + c1*f*f*(0.5 + c2*f) +
           nr9*(T/tht^3*slb_plg(tht) - T0_SLB/tht0^3*slb_plg(tht0)) + fel + Pbar*v

    # Ordering/magnetic modifiers; zero for the endmembers that carry none
    return G + slb_modifier(d, Pbar, T)
end

# O₂ standard state: sb24 keeps O₂ on a HP Cp polynomial fitted to JANAF, same five constants as
# `HP98.O2`. The T coefficient sign is MAGEMin's (a negated slot); matching it keeps ΔFMQ/ΔIW exact.
const SLB_O2 = (S0=205.3614, c1=48.3, c2=6.91e-4, c3=499200.0, c5=-420.7)
function slb_gibbs_O2(T)
    Tr = 298.15; d = SLB_O2
    cp = d.c1*(T - Tr) + d.c2*(T^2 - Tr^2)/2 - d.c3*(1/T - 1/Tr) + 2d.c5*(√T - √Tr) -
         T*(d.c1*log(T/Tr) + d.c2*(T - Tr) - d.c3/2*(T^-2 - Tr^-2) - 2d.c5*(T^-0.5 - Tr^-0.5))
    return -(T - Tr)*d.S0 + cp
end

# ΔG of each buffer reaction per mole of O₂: oxidised side minus reduced side, O₂ itself excluded.
slb_ΔG_FMQ(Pbar, T) = 3slb_gibbs(SLB_EM.qtz, Pbar, T) + 2slb_gibbs(SLB_EM.mag, Pbar, T) -
                      3slb_gibbs(SLB_EM.fa, Pbar, T)                  # 3fa + O₂ = 3qtz + 2mag
# Lowest-energy iron polymorph (α/γ/ε) and wüstite spin state, so the buffer stays valid across the
# whole mantle range rather than only where one of them happens to be the stable one.
slb_ΔG_IW(Pbar, T)  = (min(slb_gibbs(SLB_EM.wu, Pbar, T), slb_gibbs(SLB_EM.wuls, Pbar, T)) -
                       4minimum(slb_gibbs(SLB_EM[n], Pbar, T) for n in (:fea, :feg, :fee)))/2
                                                                      # 4Fe + 2O₂ = Fe₄O₄

# μ(O₂) = G⁰(O₂) + RT ln fO2, so log₁₀fO2 on the buffer is (ΔG - G⁰(O₂))/(RT ln10).
logfO2_FMQ(T, Pbar) = (slb_ΔG_FMQ(Pbar, T) - slb_gibbs_O2(T))/(R*T*log(10))
logfO2_IW(T, Pbar)  = (slb_ΔG_IW(Pbar, T)  - slb_gibbs_O2(T))/(R*T*log(10))

# The pair `partition_Oₑₓ` consumes. `qfm` is MAGEMin's own fO2 - dQFM off the solid scan; NaN (no
# solid, so no minimisation) falls back to the Gibbs form above, which agrees to <0.01 log units.
function reference_buffers(P, T; qfm=NaN)
    Pbar = 1e4P
    return (isfinite(qfm) ? qfm : logfO2_FMQ(T, Pbar)), logfO2_IW(T, Pbar)
end

# Frost (1991) Rev. Mineral. 25 Table 1 (P in bar) — the retired empirical curves, kept for the
# benchmark scripts. NOT used: subtracting them from an sb24 fO2 is a category error, see above.
logfO2_FMQ_frost(T, Pbar) = -25096.3/T + 8.735 + 0.110*(Pbar - 1)/T
logfO2_IW_frost(T, Pbar)  = -27489.0/T + 6.702 + 0.055*(Pbar - 1)/T

# ---- Section 5 — smooth saturation and finiteness helpers ----

@inline function softmax2(a, b, β)
    m = max(a, b)
    return m + log(exp(β*(a-m)) + exp(β*(b-m)))/β
end
@inline softmin2(a, b, β) = -softmax2(-a, -b, β)

# Smooth ceiling: identity well below `cap`, rolling over without a kink. Caps from ABOVE only —
# rounding the bottom would put a phantom floor in XCO₂ that Φ amplifies into a spurious demand.
@inline softcap(x, cap; sharpness=30.0) = max(zero(x), softmin2(x, cap, sharpness/cap))

# Throws on Inf as well as NaN: an Inf compares silently and would otherwise poison the solve.
@inline function nancheck(label, val, inputs...)
    isfinite(val) || error("$label evaluated to $val (non-finite) — inputs: $inputs")
    return val
end

# ---- Section 6 — the residual (1 unknown: the shared log₁₀ fO2) ----

# Inverts an increasing oxybarometer for the loading that reads `target`, saturating at the bounds.
function _invert_demand(g, target, lo, hi; iters=80)
    hi <= lo && return lo             # zero-capacity sink (phase absent): nothing to invert
    g(lo) >= target && return lo      # floor already too oxidized for this fO2 → holds nothing
    g(hi) <= target && return hi      # cap still too reduced → this sink is saturated
    for _ in 1:iters
        mid = 0.5*(lo + hi)
        g(mid) < target ? (lo = mid) : (hi = mid)
    end
    return 0.5*(lo + hi)
end

# ---- Section 6b — regime diagnosis (reporting only) ----

# Labels the converged state. Never an input to the solve, so it cannot disagree with what was solved.
function classify_regime(st, ctx, res)
    atcap(x, cap)  = cap > 0 && x >= 0.999*cap
    atfloor(x)     = x <= 1.000001*ctx.lowclip
    # Zero capacity means ABSENT (ϕ = 0 or 1), not empty: no oxybarometer, so its fO2 is meaningless.
    solid_absent = ctx.maxsOₑₓ <= 0;  melt_absent = ctx.maxmOₑₓ <= 0
    # The solid's lower end is total METAL saturation, not "holds ~no Fe³⁺" — a very different rock
    # from one sitting at the FeO reference. The melt keeps the near-zero floor; it has no metal branch.
    solid_sat   = atcap(st.sOₑₓ, ctx.maxsOₑₓ)
    solid_metal = !solid_absent && ctx.minsOₑₓ > 0 && st.sOₑₓ <= -0.999*ctx.minsOₑₓ
    melt_sat    = atcap(st.mOₑₓ, ctx.maxmOₑₓ);  melt_floored  = !melt_absent  && atfloor(st.mOₑₓ)
    graphite   = st.cs.graphite       # free graphite/diamond coexists
    carbonate  = st.cs.carbonate      # solid carbonate coexists (above EMOG(D))

    regime = if !res.converged        ; :capacity_exceeded
             elseif solid_absent          ; :melt_only
             elseif melt_absent           ; :solid_only
             elseif solid_sat && melt_sat ; :oxidant_exhausted
             elseif solid_sat             ; :solid_saturated
             elseif melt_sat              ; :melt_saturated
             elseif solid_metal           ; :solid_metal_saturated
             elseif melt_floored          ; :melt_depleted
             elseif carbonate && graphite ; :carbonate_graphite_buffered
             elseif carbonate             ; :carbonate_saturated
             elseif !graphite             ; :carbon_undersaturated
             else                         ; :buffered
             end
    # fO2 equality may legitimately fail only for a sink pinned at one of its two ends, or absent.
    pinned = solid_sat || solid_metal || melt_sat || melt_floored || solid_absent || melt_absent
    return (; regime, solid_sat, melt_sat, solid_metal, melt_floored, solid_absent, melt_absent,
              pinned, graphite, carbonate,
              aC = st.cs.aC, fE = st.cs.fE, fO2_equal = abs(st.fₛ - st.fₗ) <= 1e-6,
              has_fluid = ctx.Mh > 0, fluid_reducing = st.fₑₓ < 0)
end

# ---- Section 7 — bracketed bisection solver ----

# Last roots of the unified solve (ln aC, ln fH2, ln n), carried across ideality passes and fO2 trials of one cell
mutable struct RootSeed; x::Float64; y::Float64; u::Float64; end
RootSeed() = RootSeed(NaN, NaN, NaN)

# Illinois false position with a halving fallback: stays bracketed, which the tabulated solid sink
# (step-function slope) needs. Converged on the RESIDUAL alone: a collapsed bracket with |F| large is no root.
function bisect_solve(F::Function, lo::Real, hi::Real; maxiter=200, atol=1e-10, xtol=1e-14,
                      newton=true, debug=false)
    flo, fhi = F(lo), F(hi)
    # No sign change: the root is outside the bracket, so return the closer endpoint, unconverged.
    if sign(flo) == sign(fhi)
        best = abs(flo) <= abs(fhi) ? lo : hi
        return (converged=false, x=best, resnorm=min(abs(flo), abs(fhi)), iters=0)
    end
    x, fx, side = lo, flo, 0
    it = 0
    while it < maxiter
        it += 1
        if newton && (flo - fhi) != 0
            x = (flo*hi - fhi*lo)/(flo - fhi)
            (!isfinite(x) || x <= min(lo, hi) || x >= max(lo, hi)) && (x = 0.5*(lo + hi))
        else
            x = 0.5*(lo + hi)
        end
        fx = F(x)
        debug && println("  it=$it: x=$x F=$fx  bracket=[$lo, $hi] width=$(hi-lo)")
        (abs(fx) <= atol || (hi - lo) <= xtol*max(1.0, abs(x))) && break
        if sign(fx) == sign(flo)
            lo, flo = x, fx;  side == -1 && (fhi *= 0.5);  side = -1
        else
            hi, fhi = x, fx;  side == +1 && (flo *= 0.5);  side = +1
        end
    end
    return (converged = abs(fx) <= atol, x=x, resnorm=abs(fx), iters=it)
end

# ---- Section 9 — solid water sink: StagH2Oa capacity C(P, T, p, a) in wt%, from scripts/sH2O/sfO2.jl ----

const _SH2OA = Dict{String, Any}()

# Reads StagH2Oa.dat once per path. Blocks per activity, sector, composition (basalt, harzburgite, mixes)
function load_StagH2Oa(path=get(ENV, "STAGH2OA", joinpath(@__DIR__, "..", "scripts", "sH2O", "out", "aw16", "StagH2Oa.dat")))
    get!(_SH2OA, abspath(path)) do
        ln = readlines(path)
        h  = split(ln[1]); nP, nT, nA, nEM = parse.(Int, h[1:4]); xm = parse.(Float64, h[5:end])
        sec = [Tuple(parse.(Float64, split(ln[k]))) for k in 3:5]
        aw  = parse.(Float64, split(ln[6]))
        v   = parse.(Float64, split(join(ln[7:end], " ")))
        length(v) == nP*nT*3*nEM*nA || error("StagH2Oa: $(length(v)) values, expected $(nP*nT*3*nEM*nA)")
        S = Array{Float64}(undef, nP, nT, 3, nEM, nA); k = 0
        for ka in 1:nA, q in 1:3, c in 1:nEM, i in 1:nP, j in 1:nT
            S[i, j, q, c, ka] = v[k += 1]
        end
        pf = vcat(1.0, 0.0, xm); o = sortperm(pf)              # compositions by basalt fraction
        (; sectors = sec, pX = pf[o], aw, S = S[:, :, :, o, :], nP, nT)
    end
end

# C(a1 … an) at (P, T, p): bilinear in the P-T sector, linear in the basalt fraction, ends clamped
function solid_h2o_curve(tb, P, T, p)
    q = P < tb.sectors[2][1] ? 1 : P < tb.sectors[3][1] ? 2 : 3
    P0, P1, T0, T1 = tb.sectors[q]
    x = 1 + clamp((P - P0)/(P1 - P0), 0, 1)*(tb.nP - 1);  i = min(floor(Int, x), tb.nP - 1);  dx = x - i
    y = 1 + clamp((T - T0)/(T1 - T0), 0, 1)*(tb.nT - 1);  j = min(floor(Int, y), tb.nT - 1);  dy = y - j
    pc = clamp(p, tb.pX[1], tb.pX[end]); c = min(searchsortedlast(tb.pX, pc), length(tb.pX) - 1)
    dc = (pc - tb.pX[c])/(tb.pX[c+1] - tb.pX[c])
    C = zeros(length(tb.aw))
    for (ii, wi) in ((i, 1 - dx), (i + 1, dx)), (jj, wj) in ((j, 1 - dy), (j + 1, dy)), (cc, wc) in ((c, 1 - dc), (c + 1, dc))
        w = wi*wj*wc; w > 0 && (C .+= w .* @view tb.S[ii, jj, q, cc, :])
    end
    return C
end

# C(a) in wt%: log-log between the activity nodes (descending), C(1) at a >= 1, the last slope below
function solid_h2o_wt(a, C)
    a <= 0 && return 0.0
    la, lC = C.la, C.lC; l = log(a); n = length(la)
    l >= la[1] && return exp(lC[1])
    k = l <= la[n] ? n - 1 : findfirst(t -> la[t+1] <= l, 1:n-1)
    return exp(lC[k] + (lC[k+1] - lC[k])*(l - la[k])/(la[k+1] - la[k]))
end

# Coldest Sun & Yao (2026) experiment (660 °C): colder melts take the law's value there, as P takes it at SY_PMAX
const SY_TMIN = 933.15

# Fugacity coefficient of pure water at (P, T), same fluid EoS: the reference of the table activity
function pure_water_φ(P, T)
    vw = zeros(7); vw[2] = 1.0; Pbar = 1e4*P
    _, fug = fluid_fugacity(vw, Pbar, T); φ = fug[2]/Pbar
    return isfinite(φ) && φ > 0 ? φ : 1.0
end

# ---- Section 11 — residual and top-level driver (solid + melt + fluid, unified) ----

# A trial fO2 → the full state: each sink's demand, its own fO2 reading, the total. Used for both
# solving and reporting, so the reported state is exactly the one the residual saw.
function _coh_derive(f, ctx)
    (; sfO2, T, P, IDV, TOₑₓ, Mf, Mc, molfracXB, nB_per_mass, Φ, Mh, nOxideTotal, maxmOₑₓ, maxsOₑₓ,
       minsOₑₓ, lowclip, melt_H2, melt_H2O, h2_comp, XH2max, XCO2max, h2o_law, fluid_Pmax, sol) = ctx

    # The solid's bracket straddles zero (below IW it hands oxygen back); the melt has no metal branch.
    sOₑₓ = _invert_demand(s -> sfO2(s), f, -minsOₑₓ + lowclip, maxsOₑₓ - lowclip)
    mOₑₓ = _invert_demand(m -> melt_fO2(m, T, TOₑₓ, Mf, molfracXB, nB_per_mass, IDV), f,
                          lowclip, maxmOₑₓ - lowclip)

    # One joint solve: aC enforces carbon mass balance exactly, so melt XCO₂ needs no ceiling of its own.
    cs      = coh_state_u(f, T, P, Mh, Mc, nOxideTotal, TOₑₓ;
                        melt_H2=melt_H2, melt_H2O=melt_H2O, h2_comp=h2_comp, XH2max=XH2max,
                        XCO2max=XCO2max, Mf=Mf, molfracXB=molfracXB, h2o_law=h2o_law,
                        fluid_Pmax=fluid_Pmax, φ0=ctx.φseed[], sol=sol, seed=ctx.rseed)
    cs.φ_seed === nothing || (ctx.φseed[] = cs.φ_seed)
    XCO₂    = cs.XCO2_melt
    carbonO = Φ*XCO₂/(1 - XCO₂)

    fₑₓ   = nancheck("coh_state.fₑₓ", cs.fₑₓ, f, Mh)
    carbO = nancheck("coh_state.carbO", cs.carbO, f, Mc)   # solid carbonate: 2 O per carbon
    fₛ    = nancheck("sfO2(sOₑₓ)", sfO2(sOₑₓ), sOₑₓ)
    fₗ    = nancheck("melt_fO2", melt_fO2(mOₑₓ, T, TOₑₓ, Mf, molfracXB, nB_per_mass, IDV), mOₑₓ)

    # Melt H₂ is the sixth sink and monotone like the others: fO2 up ⇒ f_H2 down ⇒ less dissolved H₂,
    # so its deficit shrinks toward zero and the unique root survives.
    mH2ex = nancheck("coh_state.meltH2ex", cs.meltH2ex, f, Mh)
    return (; f, sOₑₓ, mOₑₓ, XCO₂, carbonO, carbO, fₑₓ, mH2ex, fₛ, fₗ, cs,
              total = sOₑₓ + mOₑₓ + carbonO + carbO + fₑₓ + mH2ex)
end

# The one equation, in grams: total sink demand equals the budget. The budget appears only here,
# so its sign and nearness to zero never reach the sink brackets (math S.3).
unified_residual(f, ctx) = _coh_derive(f, ctx).total - ctx.TOₑₓ


# Partitions the Oₑₓ budget at (P,T); `out` returns a NamedTuple. `p`: 0 harzburgitic, 1 basaltic. `H` is the
# TOTAL water (bound + free) shared by solid, melt and fluid at one activity; solid C(P,T,p,a) from `sH2Oa`.
function partition_Oₑₓ(P::Real, T::Real, p::Real, ϕ::Real, TOex::Real, TC::Real, H::Real;
                           nr=100, Rmax=0.4, nrm=15, R0=2e-3, maxiter=500, tol=1e-8,
                           tol_floor=1e-12, verbose=true,
                           data=nothing, out=false, debug=false, mag_cache=nothing,
                           melt_H2=true, melt_H2O=true, h2_comp=:basalt, XH2max=0.3,
                           XCO2max=0.5, h2o_law=:sunyao, fluid_Pmax=50.0, sH2Oa=nothing, T_quench=298.15,
                           closure_tol=1e-6)
    # Floor at the HP98 reference T, below which the gas Cp data are not defined (the fluid EoS has its own cold branch)
    T = max(T, T_quench)

    (; X, XB, Xox, SymXox, mmXox) = endmember_composition(p)
    Mt = 1e3   # reference cell mass (1 kg in grams, matching mm.* in g/mol)
    (; Mf, Ms, Mc, Mh, TOₑₓ, molfracXB, nB_per_mass, maxXCO₂_raw, Φ) = mass_budget(Mt, ϕ, TC, TOex, H, XB, mmXox)

    # Solid water sink: the cell's C(a) curve, and pure water's ZD09 φ (1 above fluid_Pmax, where the
    # fluid is ideal pure water) to turn the fluid's f_H2O into the table's activity
    tbl  = sH2Oa === nothing ? load_StagH2Oa() : sH2Oa
    Csat = solid_h2o_curve(tbl, P, T, p)
    sol  = (; Ms, C = (la = log.(tbl.aw), lC = log.(max.(Csat, 1e-12))), φw0 = P <= fluid_Pmax ? pure_water_φ(P, T) : 1.0)

    IDV = nancheck("solve_∫ΔVdP", solve_∫ΔVdP(P, T), P, T)   # kJ/mol, 1 bar to P

    # Ceilings at Fe³⁺/ΣFe = 1: moles of Fe, half an O each, over the budget. NOT clipped at 100% of
    # TOₑₓ — that is the residual's job, and doing it twice would break fO2 equality.
    maxsOₑₓ = 0.5*(Ms*X[3]/mmXox[3])*mm.O

    nOxideTotal = Mf*nB_per_mass
    maxmOₑₓ = 0.5*molfracXB[3]*nOxideTotal*mm.O
    maxXCO₂ = clamp(maxXCO₂_raw, 0.0, 1.0)   # all carbon dissolved in the melt; reported, not imposed
    lowclip = 1e-10

    # Bounds are stoichiometric, not tabulation limits — `sfO2` extrapolates past both ends.
    # `tbl_lo`/`tbl_hi` are kept only so the report can flag an answer that relies on extrapolation.
    sfO2, tbl_lo, tbl_hi, tbl_mono, s_floor, qfm_out, smetal = solid_fO2_interpolant(P, T, X, Xox, Ms, TOₑₓ, maxsOₑₓ;
                                                            nr=nr, Rmax=Rmax, nrm=nrm, R0=R0, data=data, verbose=verbose,
                                                            mag_cache=mag_cache)

    # The solid's FLOOR: every Fe²⁺ reduced to Fe⁰, handing all its iron-bound oxygen back. Taken
    # from the interpolant, not as 2·maxsOₑₓ. The melt has no counterpart — its floor is `lowclip`.
    minsOₑₓ = -s_floor

    # Both melt-H terms default ON. They are live even with no carbon (a pure H-O fluid is nearly
    # all H₂ at low fO2), so switch them OFF to reproduce a pre-2026-09-04 benchmark.
    ctx = (; sfO2, T, P, IDV, TOₑₓ, Mf, Mc, molfracXB, nB_per_mass, Φ, Mh, nOxideTotal, maxmOₑₓ,
             maxsOₑₓ, minsOₑₓ, lowclip, melt_H2, melt_H2O, h2_comp, XH2max, XCO2max, h2o_law,
             fluid_Pmax, sol, φseed = Ref{Union{Nothing,Vector{Float64}}}(nothing),   # φ seed of the ideality loop, reset per cell as in Fortran
             rseed = RootSeed())                           # and the inner roots' seeds, likewise

    # Bracket far wider than physical; demands saturate outside their range, so it only straddles the root.
    # The residual is in grams, so `tol` is relative: scaled by the budget, floored by sink capacity (S.12).
    atol_g = tol*max(abs(TOₑₓ), 1e-6*max(maxsOₑₓ + s_floor, maxmOₑₓ, floatmin(Float64)))
    # Absolute floor on the tolerance: at a zero budget the relative rule collapses to
    # 1e-14·capacity, below the inner solves' own gram-level noise. tol_floor = 0 disables it.
    atol_g = max(atol_g, tol_floor*(Mc + Mh + max(maxsOₑₓ + s_floor, maxmOₑₓ)))
    res  = bisect_solve(f -> unified_residual(f, ctx), -120.0, 40.0; maxiter=maxiter, atol=atol_g, debug=debug)
    st   = _coh_derive(res.x, ctx)
    # Converged only if the H and C balances close as well, and no hydrogen floor fell on a melt
    conv = res.converged && st.cs.res_ideality <= closure_tol && !st.cs.h_floor_melt
    sOₑₓ, mOₑₓ, XCO₂ = st.sOₑₓ, st.mOₑₓ, st.XCO₂

    # Oxidation state at the loading actually solved, so it cannot drift from the answer.
    # Metal comes from MAGEMin, not from sOₑₓ, so ferric and metal can now both be non-zero at once.
    nFe_solid = Ms > 0 ? Ms*X[3]/mmXox[3] : 0.0          # moles Fe in the solid residue
    nFe_melt  = molfracXB[3]*nOxideTotal                 # moles Fe in the melt
    nO_solid  = st.sOₑₓ/mm.O
    FeR_melt  = nFe_melt  > 0 ? 2*(st.mOₑₓ/mm.O)/nFe_melt  : 0.0
    # Native iron from MAGEMin's assemblage, so it includes disproportionated metal.
    # `FeM_net` is the blind coordinate's own answer, a strict LOWER BOUND on FeM_solid.
    metal_wt   = Ms > 0 ? smetal(st.sOₑₓ) : 0.0                        # wt fraction of the solid
    nFe0_solid = metal_wt*Ms/mm.Fe                                     # moles of native iron
    FeM_solid  = nFe_solid > 0 ? nFe0_solid/nFe_solid : 0.0            # Fe⁰/ΣFe, total
    FeM_net    = nFe_solid > 0 ? max(0.0, -nO_solid/nFe_solid) : 0.0   # Fe⁰/ΣFe, net-oxygen only
    FeR_net   = nFe_solid > 0 ? max(0.0, 2nO_solid/nFe_solid) : 0.0
    # The two routes, split by the oxygen budget since MAGEMin reports only a total: dissociation
    # frees exactly 1 O per Fe⁰ and that O IS -sOₑₓ, so the rest was made at constant oxygen.
    FeM_diss = min(FeM_net, FeM_solid)                                 # paid for by O leaving the solid
    FeM_disp = max(0.0, FeM_solid - FeM_diss)                          # made at constant oxygen
    # IRON BALANCE. 3 Fe²⁺ → 2 Fe³⁺ + Fe⁰, so the disproportionated metal cannot exceed a third
    # of the Fe²⁺ left after the net oxygen and the dissociated metal have taken theirs.
    R_net      = nFe_solid > 0 ? 2nO_solid/nFe_solid : 0.0
    disp_cap   = max(0.0, (1 - max(R_net, 0.0) - FeM_diss)/3)
    fe_clipped = FeM_disp - min(FeM_disp, disp_cap) > 1e-9
    FeM_disp   = min(FeM_disp, disp_cap)
    FeM_solid  = FeM_diss + FeM_disp
    # With metal present the net oxygen alone no longer gives Fe³⁺: balance is nFe³⁺ = 2(nO + nFe⁰),
    # which returns the 2:1 ratio at nO = 0 and zero when metal is exactly paid for by released O.
    FeR_solid  = min(1.0, max(0.0, R_net + 2FeM_solid))

    eqfO2     = st.f       # the shared fO2 all unpinned sinks equilibrate to
    cohsink   = st.fₑₓ
    # A FRACTION of the budget, and the budget can legitimately be zero.
    shortfall = TOₑₓ != 0 ? (TOₑₓ - st.total)/TOₑₓ : 0.0

    cs = st.cs
    # Both out of sb24, the dataset that produced `eqfO2`: FMQ off MAGEMin's own dQFM, IW from its
    # wüstite/iron endmembers. Section 4c explains why an empirical curve cannot go here.
    fmq, iw = reference_buffers(P, T; qfm=qfm_out)
    if verbose
        XH2O, XCO2f, XCH4, XCOf, XH2 = cs.XH2O, cs.XCO2, cs.XCH4, cs.XCO, cs.XH2
        rg = classify_regime(st, ctx, res)
        # Layout helpers. No boxes: two columns are set by padding alone, so a long line runs on.
        L, C = 4, 48                                  # left margin, second-column offset
        pad(s, w) = s * " "^max(0, w - length(s))
        # rstrip everywhere: a column of trailing blanks is invisible until someone diffs it.
        say(s)    = println(rstrip(" "^(L+2) * s))
        head(t)   = (println(); println(" "^L * t); println(" "^L * "─"^length(t)))
        two(a, b) = println(rstrip(" "^(L+2) * pad(a, C) * b))
        function twohead(a, b)
            println()
            println(" "^L * pad(a, C+2) * b)
            println(" "^L * pad("─"^length(a), C+2) * "─"^length(b))
        end

        # ══ header ═════════════════════════════════════════════════════════════════════
        println()
        println(" "^L * "═"^86)
        @printf("%sP %.4g GPa    T %.5g K    ϕ %.4g %%    mix %.2f        TOex %.4g wt%%    TC %.4g wt%%    H %.4g  wt%%\n",
                " "^(L+1), P, T, 1e2ϕ, p, 1e2TOex, 1e2TC, 1e2H)
        println(" "^L * "═"^86)
        println()
        @printf("%s%s in %d iterations       |F| %.3g %% of the budget\n",
                " "^(L+1), conv ? "CONVERGED" : "DID NOT CONVERGE", res.iters,
                TOₑₓ != 0 ? 1e2*abs(st.total-TOₑₓ)/abs(TOₑₓ) : 0.0)
        @printf("%sregime   %s\n", " "^(L+1), rg.regime)
        if !conv
            println()
            println(" "^(L+1) * "⚠  CAPACITY EXCEEDED — even at every cap the assemblage cannot hold the budget.")
            shortfall >  1e-6 && @printf("%s   %.4g %% of TOₑₓ has no sink to go into.\n", " "^(L+1),  1e2shortfall)
            shortfall < -1e-6 && @printf("%s   the sinks demand %.4g %% more than TOₑₓ holds.\n", " "^(L+1), -1e2shortfall)
        end
        (T < 673 || T > 2573 || P > 10.0) && println(" "^(L+1) *
            "⚠  (T,P) outside the Zhang & Duan (2009) COH EOS box — fluid speciation extrapolated.")

        # ══ everything about fO2, in one place ═════════════════════════════════════════
        head("OXYGEN FUGACITY")
        @printf("%slog fO₂  %10.4f            ΔFMQ  %+8.3f            ΔIW  %+8.3f\n",
                " "^(L+2), eqfO2, eqfO2-fmq, eqfO2-iw)
        println()
        spin = rg.solid_sat   ? "  pinned at its ferric cap" :
               rg.solid_metal ? "  pinned at the metal floor — all Fe as Fe⁰" : ""
        mpin = rg.melt_sat    ? "  pinned at its ferric cap" :
               rg.melt_floored ? "  pinned at its floor" : ""
        sread = rg.solid_absent ? @sprintf("%10s", "— no solid") : @sprintf("%10.4f%s", st.fₛ, spin)
        mread = rg.melt_absent  ? @sprintf("%10s", "— no melt")  : @sprintf("%10.4f%s", st.fₗ, mpin)
        @printf("%ssolid reads %s\n", " "^(L+2), sread)
        @printf("%smelt  reads %s\n", " "^(L+2), mread)
        println()
        if rg.solid_absent || rg.melt_absent
            say("one phase is absent, so the remaining sinks carry the whole budget.")
        elseif rg.fO2_equal
            say("✓ solid and melt agree — the equilibrium the whole model rests on.")
        elseif rg.pinned
            @printf("%sthey differ by %.4g, which is legitimate: a pinned sink has left the\n", " "^(L+2), abs(st.fₛ-st.fₗ))
            say("equilibrium and no longer defines an fO₂.")
        else
            @printf("%s⚠ they differ by %.4g with NEITHER pinned. Coexisting phases must share one\n", " "^(L+2), abs(st.fₛ-st.fₗ))
            say("fO₂ — treat this result as suspect.")
        end

        # ══ the two IRON sinks only; the melt's volatiles get their own section so the columns stay parallel
        unphys(x) = x > 1.05 ? "  ⚠ >1" : ""
        shr(x)    = TOₑₓ != 0 ? @sprintf("%+.2f %%", 1e2x/TOₑₓ) : "n/a"
        cap(x, c) = c > 0 ? @sprintf("%.2f %% of its cap", 1e2x/c) : "n/a"
        SL, ML = String[], String[]
        if rg.solid_absent
            push!(SL, "absent (ϕ = 1)")
        else
            push!(SL, @sprintf("Fe³⁺/ΣFe   %.6f%s", FeR_solid, unphys(FeR_solid)))
            push!(SL, @sprintf("Fe⁰/ΣFe    %.6f   (%.4g wt%%)", FeM_solid, 1e2metal_wt))
            FeM_solid > 1e-12 && push!(SL, @sprintf("   dissociation      %.6f", FeM_diss))
            FeM_solid > 1e-12 && push!(SL, @sprintf("   disproportionation %.6f", FeM_disp))
            push!(SL, "")
            push!(SL, @sprintf("holds      %s", shr(sOₑₓ)))
            push!(SL, @sprintf("           %s", cap(sOₑₓ, maxsOₑₓ)))
            push!(SL, @sprintf("can hold   %s … %s", shr(-minsOₑₓ), shr(maxsOₑₓ)))
            push!(SL, "")
            inside = tbl_lo ≤ sOₑₓ ≤ tbl_hi
            push!(SL, @sprintf("table      %s … %s", shr(tbl_lo), shr(tbl_hi)))
            push!(SL, inside ? "           ✓ solution inside" :
                      sOₑₓ > tbl_hi ? "           ⚠ above — extrapolated" :
                                      "           ⚠ below — held flat")
            tbl_mono || push!(SL, "           ⚠ non-monotonic, repaired")
        end
        if rg.melt_absent
            push!(ML, "absent (ϕ = 0)")
        else
            push!(ML, @sprintf("Fe³⁺/ΣFe   %.6f%s", FeR_melt, unphys(FeR_melt)))
            push!(ML, "")
            push!(ML, "")
            push!(ML, @sprintf("holds      %s", shr(mOₑₓ)))
            push!(ML, @sprintf("           %s", cap(mOₑₓ, maxmOₑₓ)))
            push!(ML, @sprintf("can hold   0 … %s", shr(maxmOₑₓ)))
            push!(ML, "")
            push!(ML, "no metal branch — a melt has no Fe⁰")
        end
        twohead("SOLID  iron", "MELT  iron")
        for i in 1:max(length(SL), length(ML))
            two(i ≤ length(SL) ? SL[i] : "", i ≤ length(ML) ? ML[i] : "")
        end

        # ══ volatiles: what the melt dissolved, then what the fluid is ═════════════════
        if !rg.melt_absent
            head("DISSOLVED IN THE MELT")
            wtCO2 = 1e2*cs.nC_melt*(mm.C + 2mm.O)/max(Mf, eps())
            wtH2  = 1e2*cs.nH2_melt*2.016/max(Mf, eps())
            pH(n) = cs.molH_tot > 0 ? 1e2*n/cs.molH_tot : 0.0
            wtot = cs.SH2Om + cs.SOH
            # The molecular : hydroxyl split rides on the H₂O line as a ratio. Two numbers that
            # only mean anything together do not deserve a row of their own.
            sp = (melt_H2O && wtot > 0) ?
                 @sprintf("   mol %.0f / OH⁻ %.0f", 1e2cs.SH2Om/wtot, 1e2cs.SOH/wtot) : ""
            @printf("%s%-5s %8.4g wt%%  %8.4g mol%%   %s%s\n", " "^(L+2), "H₂O", wtot,
                    1e2cs.XH2O_melt,
                    cs.molH_tot > 0 ? @sprintf("%5.1f %% of cell H", pH(cs.nH_H2O_melt)) : "", sp)
            @printf("%s%-5s %8.4g wt%%  %8.4g mol%%   %s\n", " "^(L+2), "H₂", wtH2,
                    1e2cs.XH2_melt,
                    cs.molH_tot > 0 ? @sprintf("%5.1f %% of cell H", pH(2*cs.nH2_melt)) : "")
            @printf("%s%-5s %8.4g wt%%  %8.4g mol%%   %s\n", " "^(L+2), "CO₂", wtCO2, 1e2XCO₂,
                    maxXCO₂ > 0 ? @sprintf("%5.2f %% of its cap", 1e2XCO₂/maxXCO₂) : "")
            if melt_H2 || melt_H2O
                println()
                @printf("%sP_H₂O %.4g bar        f_H₂ %.4g bar\n", " "^(L+2), cs.PH2O_melt, cs.fH2)
            end
        end

        head("C–O–H FLUID")
        if cs.co_fluid
            @printf("%shydrogen-free CO₂–CO fluid     XCO₂ %.5f     XCO %.3e\n", " "^(L+2), XCO2f, XCOf)
            @printf("%sholds %.2f %% of the cell carbon at aC %.4g\n", " "^(L+2),
                    1e2*cs.nC_fluid/max(Mc/mm.C, eps()), cs.aC)
            println()
            say("net OXIDISING — a SINK: this fluid takes oxygen from the budget.")
        elseif Mh <= 0.0
            say("absent — the cell holds no hydrogen, and carbon stays in the solid or the melt.")
        elseif cs.fluid_absent
            say("absent — melt and solid hold the whole hydrogen inventory.")
            @printf("%sfugacities sit at λ = %.4g of saturation; below 1 nothing exsolves.\n", " "^(L+2), cs.lam)
        else
            @printf("%sXH₂O %.5f     XCO₂ %.5f     XCO %.3e     XCH₄ %.3e     XH₂ %.3e\n",
                    " "^(L+2), XH2O, XCO2f, XCOf, XCH4, XH2)
            println()
            @printf("%shydrogen → fluid %.2f %%,  melt %.2f %%,  solid %.2f %%\n", " "^(L+2),
                    1e2*cs.molH_fl/max(cs.molH_tot,eps()),
                    1e2*(1-(cs.molH_fl + cs.nH_solid)/max(cs.molH_tot,eps())), 1e2*cs.nH_solid/max(cs.molH_tot,eps()))
            println()
            say((rg.fluid_reducing ?
                "net REDUCING — a SOURCE: this fluid hands oxygen back to the rock." :
                "net OXIDISING — a SINK: this fluid takes oxygen from the budget."))
        end
        println()
        say((rg.carbonate && rg.graphite ? "graphite/diamond and carbonate coexist, on EMOG(D)" :
                            rg.carbonate ? @sprintf("solid carbonate present   aC %.4g", rg.aC) :
                            rg.graphite  ? "graphite/diamond present   aC 1" :
                                           @sprintf("carbon-undersaturated, no solid carbon   aC %.4g", rg.aC)))

        # ══ budgets. ONLY the rows that can go either way are annotated. ═══════════════
        gp(x) = TOₑₓ != 0 ? 1e2x/TOₑₓ : 0.0
        head("OXYGEN BUDGET")
        @printf("%sthe sinks account for %+.4f %% of the entering budget  (residual %+.3g %%)\n",
                " "^(L+2), gp(st.total), gp(st.total-TOₑₓ))
        println()
        orow(lbl, g, note) = println(rstrip(@sprintf("%s%-20s %+9.3f %%   %s", " "^(L+2), lbl, gp(g), note)))
        orow("solid Fe₂O₃",     sOₑₓ, rg.solid_absent ? "" :
             sOₑₓ < 0 ? "SOURCE — below FeO, Fe⁰ hands oxygen back" : "")
        orow("melt Fe₂O₃",      mOₑₓ, "")
        orow("melt CO₂",        st.carbonO, "")
        orow("melt H₂",         st.mH2ex, "")
        orow("solid carbonate", st.carbO, "")
        orow("C–O–H fluid",     cohsink, cohsink > 0 ? "SINK — oxidising carbon to CO₂/CO" :
                                         cohsink < 0 ? "SOURCE — reducing H₂O to CH₄/H₂" : "")
        say("─"^33)
        @printf("%s%-20s %+9.3f %%\n", " "^(L+2), "", gp(st.total))
        cohsink < 0 && @printf("\n%sH₂O dissociation ADDED %.3f %%; the rock shares out %.3f %%\n",
                               " "^(L+2), -gp(cohsink), 1e2 - gp(cohsink))
        println()
        say("signs: melt and carbonate rows are always sinks; melt H₂ always a source.")

        head("CARBON BUDGET")
        if Mc > 0
            crow(lbl, n) = @printf("%s%-20s %+9.3f %%\n", " "^(L+2), lbl, 1e2*n*mm.C/Mc)
            crow("solid graphite/dia", cs.nC_gph)
            crow("solid carbonate",    cs.nC_carb)
            crow("melt CO₂",           cs.nC_melt)
            crow("C–O–H fluid",        cs.nC_fluid)
            tot = cs.nC_gph + cs.nC_carb + cs.nC_melt + cs.nC_fluid
            say("─"^33)
            @printf("%s%-20s %+9.3f %%\n", " "^(L+2), "", 1e2*tot*mm.C/Mc)
        else
            say("TC = 0 — no carbon in this tracer.")
        end

        head("HYDROGEN BUDGET")
        if Mh > 0
            hrow(lbl, n) = @printf("%s%-20s %+9.3f %%\n", " "^(L+2), lbl, 1e2*n/max(cs.molH_tot,eps()))
            hrow("C–O–H fluid",  cs.molH_fl)
            hrow("melt H₂O",     cs.nH_H2O_melt)
            hrow("melt H₂",      2*cs.nH2_melt)
            hrow("solid H₂O",    cs.nH_solid)
            say("─"^33)
            @printf("%s%-20s %+9.3f %%\n", " "^(L+2), "",
                    1e2*(cs.molH_fl + cs.nH_H2O_melt + 2cs.nH2_melt + cs.nH_solid)/max(cs.molH_tot,eps()))
            println()
            @printf("%ssolid a_H₂O %.4g   holds %.4g wt%% of %.4g at a = 1   (φ°_H₂O %.4g)\n", " "^(L+2),
                    cs.a_solid, solid_h2o_wt(cs.a_solid, sol.C), Csat[1], sol.φw0)
        else
            say("H = 0 — no hydrogen in this tracer.")
        end
        println()
    end

    # Interface in SHARES, solve in masses: oxygen as a fraction of TOₑₓ, carbon of Mc, hydrogen of molH.
    # At TOₑₓ = 0 a share is undefined, so oS returns 0 and the gram values below carry the answer.
    oS(x) = TOₑₓ != 0 ? x/TOₑₓ : 0.0
    fC(x) = Mc    > 0 ? x*mm.C/Mc : 0.0
    fH(x) = cs.molH_tot > 0 ? x/cs.molH_tot : 0.0
    # cohsink is SIGNED (negative ⇒ the fluid donates); Oex_added is what H₂O dissociation ADDED and
    # must be booked against hydrogen. rock_pool = 1 - cohsink, NOT 1 + Oex_added, which over-reports.
    if out; return (converged=conv, fO2=eqfO2, FMQ=eqfO2-fmq, IW=eqfO2-iw, aC=st.cs.aC, sOₑₓ=oS(sOₑₓ), mOₑₓ=oS(mOₑₓ), XCO₂=XCO₂,
                    carbonO=oS(st.carbonO), carbO=oS(st.carbO), cohsink=oS(cohsink),
                    fluid_O_uptake=oS(cs.O_up), fluid_O_release=oS(cs.O_rel),
                    Oex_added=max(0.0, -oS(cohsink)), rock_pool=1.0-oS(cohsink),
                    FeR_solid=FeR_solid, FeR_net=FeR_net, FeM_solid=FeM_solid, FeM_net=FeM_net,
                    FeM_diss=FeM_diss, FeM_disp=FeM_disp, metal_wt=metal_wt, fe_clipped=fe_clipped,
                    nFe_solid=nFe_solid, FeR_melt=FeR_melt,
                    fC_gph=fC(st.cs.nC_gph), fC_carb=fC(st.cs.nC_carb),
                    fC_melt=fC(st.cs.nC_melt), fC_fluid=fC(st.cs.nC_fluid),
                    fH_fluid=fH(cs.molH_fl), fH_melt=fH(cs.molH_tot - cs.molH_fl - cs.nH_solid),
                    fH_solid=fH(cs.nH_solid), aH2O_solid=cs.a_solid, C_solid=solid_h2o_wt(cs.a_solid, sol.C),
                    C_solid_sat=Csat[1], φ_H2O_pure=sol.φw0,
                    iters=res.iters, resnorm=res.resnorm, res_id=cs.res_ideality, it_id=cs.it_ideality,
                    h_floor=cs.h_floor, h_floor_melt=cs.h_floor_melt,
                    XH2O=cs.XH2O, XCO2=cs.XCO2, XCH4=cs.XCH4, XCO=cs.XCO, XH2=cs.XH2,
                    meltH2ex=oS(st.mH2ex), nH2_melt=cs.nH2_melt, XH2_melt=cs.XH2_melt, fH2=cs.fH2,
                    molH_fluid=cs.molH_fl, molH_total=cs.molH_tot,
                    XH2O_melt=cs.XH2O_melt, nH2O_melt=cs.nH2O_melt, SH2O_melt=cs.SH2Om+cs.SOH,
                    SH2Om=cs.SH2Om, SOH=cs.SOH, PH2O_melt=cs.PH2O_melt, fluid_absent=cs.fluid_absent,
                    fluid_restricted=cs.fluid_restricted, vol_fluid=cs.vol_fluid,
                    lam=cs.lam,
                    # GRAMS, the solve's own units. The only meaningful oxygen output at TOₑₓ = 0,
                    # where every share above divides by a zero budget.
                    TOₑₓ_g=TOₑₓ, sOₑₓ_g=sOₑₓ, mOₑₓ_g=mOₑₓ, carbonO_g=st.carbonO, carbO_g=st.carbO,
                    cohsink_g=cohsink, meltH2ex_g=st.mH2ex, atol_g=atol_g,
                    mC_gph=st.cs.nC_gph*mm.C, mC_carb=st.cs.nC_carb*mm.C,
                    mC_melt=st.cs.nC_melt*mm.C, mC_fluid=st.cs.nC_fluid*mm.C,
                    mH_fluid=cs.molH_fl*0.5*mm.H2O,
                    mH_melt=(cs.molH_tot - cs.molH_fl - cs.nH_solid)*0.5*mm.H2O,
                    mH_solid=cs.nH_solid*0.5*mm.H2O); end

end

# ---- Section 10 — unified C-O-H fluid solve: Newton on ln n around a 2D Newton on (ln aC, ln f_H2) ----
# The fluid amount n closes S = ΣX ≤ 1 with n·(1 - S) = 0; carbon and hydrogen close at every trial n.
#  Derivatives by ForwardDiff through the closures below; bracketed fallbacks keep every root safe.


# Everything one fO2 trial and one φ pass hold fixed, concrete-typed so the closures compile to straight code
struct UTrial
    fl_ok::Bool                 # P ≤ fluid_Pmax: a C-O-H fluid; above it pure ideal water, y = ln a_H2O
    c::NTuple{5,Float64}        # X_H2O = c1·e^y, X_CO2 = c2·e^x, X_CH4 = c3·e^(x+2y), X_CO = c4·e^x, X_H2 = c5·e^y
    caw::Float64                # a_H2O = caw·fH2 when a fluid is allowed
    nC::Float64                 # carbon budget [mol]
    nH::Float64                 # hydrogen budget [mol H]
    nOx::Float64                # oxide moles of the melt
    Mf::Float64                 # melt mass [g]
    Ms::Float64                 # solid mass [g]
    Bmelt::Float64              # melt CO2 per unit aC (inverted carbon buffer)
    XCO2cap::Float64            # melt CO2 ceiling, mole fraction
    x_cap::Float64              # ln aC at graphite/diamond or carbonate saturation (EMOG)
    wC::Float64                 # e-folding width of the solid-carbon ramp below the cap, in ln aC
    melt_H2::Bool               # H2 dissolved in the melt
    KH2::Float64                # its Henry constant: X_H2 = KH2·fH2 below the cap
    XH2max::Float64             # melt H2 ceiling, mole fraction
    melt_H2O::Bool              # H2O + OH dissolved in the melt
    katz::Bool                  # melt-water law: Katz (true) or Sun & Yao (false)
    kSm::Float64                # Sun & Yao molecular H2O: ln S_m = a1·ln a + kSm
    kSo::Float64                # Sun & Yao hydroxyl:      ln S_o = a6·ln a + kSo
    kKatz::Float64              # Katz: S = kKatz·a
    sol_on::Bool                # the solid takes water from the C(a) table
    la::Vector{Float64}         # ln a of the table's activity nodes, descending
    lC::Vector{Float64}         # ln C [wt%] of this cell's curve at those nodes
    xlo::Float64                # lowest ln aC searched
    xhi::Float64                # highest ln aC searched (the cap)
    ylo::Float64                # lowest ln fH2 searched (hydrogen floor below it)
    yhi::Float64                # highest ln fH2 searched
end

# Work counters for comparisons: 2D F evaluations, 2D fallbacks, S evaluations, n bisection steps
const U_COUNT = zeros(Int, 4)

# Fluid species mole fractions at the potentials (x, y), generic so ForwardDiff can run through it
@inline function u_species(t::UTrial, x, y)
    ex, ey = exp(x), exp(y)                         # aC and fH2 (or a_H2O above fluid_Pmax)
    z = zero(ex*ey)                                 # a zero of the common (dual) type
    t.fl_ok || return SVector(ey + z, z, z, z, z)   # above fluid_Pmax: pure water at activity e^y
    return SVector(t.c[1]*ey + z,                   # H2O
                   t.c[2]*ex + z,                   # CO2
                   t.c[3]*ex*ey*ey,                 # CH4
                   t.c[4]*ex + z,                   # CO
                   t.c[5]*ey + z)                   # H2
end

# Solid water [wt%] at activity a: log-log between the activity nodes, C(1) at a ≥ 1, last slope below
@inline function u_solid(t::UTrial, a)
    a <= 0 && return zero(a)                        # no activity, no water
    la, lC = t.la, t.lC                             # table nodes, descending in ln a
    l = log(a)
    m = length(la)
    l >= la[1] && return zero(l) + exp(lC[1])       # at or above a = 1: the saturated capacity
    k = 1                                           # segment holding l
    if l <= la[m]
        k = m - 1                                   # below the table: extend the last segment
    else
        while la[k+1] > l                           # walk down to the segment [la[k+1], la[k]]
            k += 1
        end
    end
    return exp(lC[k] + (lC[k+1] - lC[k])*(l - la[k])/(la[k+1] - la[k]))   # straight line in log-log
end

# Melt water [wt% H2O, molecular + OH] at activity a
@inline function u_meltS(t::UTrial, a)
    a <= 0 && return zero(a)                        # no water activity, no dissolved water
    t.katz && return t.kKatz*a                      # Katz: linear in a
    la = log(a)
    return exp(SY_H2O[1]*la + t.kSm) + exp(SY_H2O[6]*la + t.kSo)   # Sun & Yao: two power laws in a
end

# Carbon held at ln aC = x with fluid species sp and fluid moles n [mol C]
@inline function u_cdem(t::UTrial, x, sp, n)
    ex = exp(x)                                     # aC
    cm = zero(ex)                                   # melt CO2
    if t.nOx > 0
        X  = softcap(ex*t.Bmelt, t.XCO2cap)         # melt XCO2, smoothly capped
        cm = X*t.nOx/(1 - X)
    end
    cf = n*(sp[2] + sp[3] + sp[4])                  # fluid CO2 + CH4 + CO
    cs = t.nC*exp((x <= t.x_cap ? x - t.x_cap : zero(x))/t.wC)   # solid ramp; at the cap it keeps the slope from below
    return cm + cf + cs
end

# Hydrogen held by melt and solid at ln fH2 = y, and the activity that sets them
@inline function u_hsinks(t::UTrial, y, sp)
    ey  = exp(y)
    fH2 = t.fl_ok ? ey : zero(ey)                   # no fH2 above fluid_Pmax (pure water)
    nH2 = zero(ey)                                  # melt H2 [mol]
    if t.melt_H2 && t.nOx > 0
        X   = softcap(t.KH2*fH2, t.XH2max)          # melt XH2, smoothly capped
        nH2 = X*t.nOx/(1 - X)
    end
    aw  = min(one(ey), t.fl_ok ? t.caw*ey : sp[1] + zero(ey))              # water activity, at most 1
    nHw = (t.melt_H2O && t.Mf > 0) ? 2e-2*t.Mf*u_meltS(t, aw)/mm.H2O : zero(ey)    # melt H2O + OH [mol H]
    nHs = t.sol_on ? 2e-2*t.Ms*u_solid(t, aw)/mm.H2O : zero(ey)            # solid water [mol H]
    return (; nH2, nHw, nHs, aw, fH2)
end

# The two closures as log residuals in z = (x, y). A budget of zero freezes its potential at -Inf
#  and its residual becomes the dummy z = 0, so the system stays 2x2 and nonsingular in every case.
@inline function u_F(t::UTrial, z, n)
    x  = t.nC > 0 ? z[1] : oftype(z[1], -Inf) + 0*z[1]                      # ln aC (frozen without carbon)
    y  = t.nH > 0 ? z[2] : oftype(z[2], -Inf) + 0*z[2]                      # ln fH2 (frozen without hydrogen)
    sp = u_species(t, x, y)
    F1 = t.nC > 0 ? log(u_cdem(t, x, sp, n)/t.nC) : z[1]                   # carbon held / carbon budget
    s  = u_hsinks(t, y, sp)
    hf = n*(2sp[1] + 4sp[3] + 2sp[5])                                       # fluid H2O + CH4 + H2 [mol H]
    F2 = t.nH > 0 ? log((2s.nH2 + s.nHw + s.nHs + hf)/t.nH) : z[2]          # hydrogen held / hydrogen budget
    return SVector(F1 + zero(F2), F2 + zero(F1))      # one element type; zero() stays finite at F = -Inf
end

# Value and Jacobian of f at z from one dual pass, allocation-free for static vectors
@inline function u_valjac(f::F, z::SVector{N,Float64}) where {F,N}
    Tg = typeof(ForwardDiff.Tag(f, Float64))        # tag, so nested differentiation cannot mix
    zd = SVector(ntuple(i -> ForwardDiff.Dual{Tg}(z[i], ntuple(j -> Float64(i == j), Val(N))), Val(N)))   # seeded duals
    fd = f(zd)
    M  = length(typeof(fd))                         # number of outputs, known at compile time
    v  = SVector(ntuple(i -> ForwardDiff.value(fd[i]), Val(M)))                                         # f(z)
    J  = SMatrix{M,N}(ntuple(k -> ForwardDiff.partials(fd[(k-1)%M + 1], (k-1)÷M + 1), Val(M*N)))        # ∂f/∂z, column-major
    return v, J
end

# Bracketed Newton–bisection on an increasing g over [a, b]; gd(x) returns (g, g').
#  No sign change in [a, b]: the end nearer the root, flagged.
function u_rtsafe(gd::F, a, b, x0; atol=1e-12, maxit=100) where {F}
    ga, _ = gd(a)
    ga >= 0 && return (x = a, ok = abs(ga) <= atol)   # root at or below a
    gb, _ = gd(b)
    gb <= 0 && return (x = b, ok = abs(gb) <= atol)   # root at or above b
    x = (isfinite(x0) && a < x0 < b) ? x0 : 0.5*(a + b)   # the seed if it lies inside, else the middle
    for _ in 1:maxit
        g, dg = gd(x)
        abs(g) <= atol && return (x = x, ok = true)
        g < 0 ? (a = x) : (b = x)                   # shrink the bracket on the sign
        xn = x - g/dg                               # Newton step
        x  = (isfinite(xn) && a < xn < b) ? xn : 0.5*(a + b)   # outside the bracket: bisect instead
        b - a <= 1e-14*max(1.0, abs(x)) && return (x = x, ok = true)
    end
    return (x = x, ok = false)
end

# Fallback for the pair: carbon closed at every trial y, hydrogen on its total derivative. Also
#  decides the hydrogen floor: the solid's capacity at the lowest activity already exceeds the budget.
function u_xy_nested(t::UTrial, n, sd::RootSeed)
    yd(y) = t.nH > 0 ? y : 0.0                      # the dummy y when there is no hydrogen
    # x(y): carbon closed at y
    function xof(y)
        t.nC > 0 || return -Inf
        gx(x) = (vj = u_valjac(z -> u_F(t, z, n), SVector(x, yd(y))); (vj[1][1], vj[2][1,1]))
        return u_rtsafe(gx, t.xlo, t.xhi, sd.x).x
    end
    t.nH > 0 || return (x = xof(-Inf), y = -Inf, hfloor = false)   # carbon only
    # g(y): hydrogen at y with carbon closed; slope F2_y + F2_x·dx/dy, dx/dy = -F1_y/F1_x
    function gy(y)
        x = xof(y)
        v, J = u_valjac(z -> u_F(t, z, n), SVector(t.nC > 0 ? x : 0.0, y))
        dxdy = t.nC > 0 ? -J[1,2]/J[1,1] : 0.0
        return (v[2], J[2,2] + J[2,1]*dxdy)
    end
    r = u_rtsafe(gy, t.ylo, t.yhi, sd.y)
    hfloor = t.sol_on && r.x <= t.ylo && gy(t.ylo)[1] > 0   # demand above budget even at the lowest fH2
    return (x = xof(r.x), y = r.x, hfloor = hfloor)
end

# 2D Newton on (x, y) at fluid moles n: steps capped at 2 log units, projected on the box, Armijo
#  backtracking on |F|². Seeds come from the cell's previous solve; a failure goes to the fallback.
function u_xy(t::UTrial, n, sd::RootSeed)
    ax, ay = t.nC > 0, t.nH > 0                     # which potentials are unknowns
    x = ax ? clamp(isfinite(sd.x) ? sd.x : t.x_cap - 1.0, t.xlo, t.xhi) : 0.0      # seed or just below the cap
    y = ay ? clamp(isfinite(sd.y) ? sd.y : t.ylo + 140.0, t.ylo, t.yhi) : 0.0      # seed or the H2O-H2 reference
    ycap = min(t.fl_ok ? -log(t.caw) : 0.0, t.yhi)  # ln fH2 where a_H2O reaches 1: D is flat in y above it
    ay && n <= 0 && (y = min(y, ycap - 1e-6))       # never start on the flat stretch
    f = z -> u_F(t, z, n)
    z = SVector(x, y)
    v, J = u_valjac(f, z)
    U_COUNT[1] += 1
    # With no fluid, hydrogen still short at a = 1 has no root: the bracket end, as the fallback
    #  returns, and S > 1 then says a fluid exists
    if ay && n <= 0 && u_valjac(f, SVector(x, ycap))[1][2] < 0
        r = u_xy_nested(t, n, sd)
        ax && isfinite(r.x) && (sd.x = r.x)
        return r
    end
    for _ in 1:40
        if maximum(abs, v) <= 1e-12                 # both closures met
            ax && (sd.x = z[1])
            ay && (sd.y = z[2])
            return (x = ax ? z[1] : -Inf, y = ay ? z[2] : -Inf, hfloor = false)
        end
        δ = -(J \ v)                                # Newton step
        all(isfinite, δ) || break
        s = maximum(abs, δ)
        s > 2 && (δ *= 2/s)                         # at most 2 log units per step
        m0, λ, acc = sum(abs2, v), 1.0, false
        for _ in 1:8                                # halve until |F|² drops enough
            zn = SVector(ax ? clamp(z[1] + λ*δ[1], t.xlo, t.xhi) : 0.0,
                         ay ? clamp(z[2] + λ*δ[2], t.ylo, t.yhi) : 0.0)
            vn, Jn = u_valjac(f, zn)
            U_COUNT[1] += 1
            if sum(abs2, vn) < (1 - 1e-4λ)*m0
                z, v, J, acc = zn, vn, Jn, true
                break
            end
            λ /= 2
        end
        acc || break                                # no descent: hand over
    end
    U_COUNT[2] += 1
    r = u_xy_nested(t, n, sd)
    ax && isfinite(r.x) && (sd.x = r.x)
    ay && isfinite(r.y) && (sd.y = r.y)
    return r
end

# Fluid amount: none if the potentials at n = 0 leave S ≤ 1, else S(n) = 1 by Newton on u = ln n,
#  with dS/du through the implicit (dx, dy)/du of the two closures; bisection keeps it bracketed.
function u_solve_n(t::UTrial, sd::RootSeed)
    nmax = t.nC + t.nH/2                            # one C or two H per fluid molecule
    nmax <= 0 && return (n = 0.0, xy = (x = -Inf, y = -Inf, hfloor = false))
    xy0 = u_xy(t, 0.0, sd)                          # potentials with no fluid
    sum(u_species(t, xy0.x, xy0.y)) <= 1 && return (n = 0.0, xy = xy0)   # undersaturated: no fluid
    a, b = log(nmax) - 70.0, log(nmax)              # bracket in u; S falls with u, so G > 0 below the root
    u  = (isfinite(sd.u) && a < sd.u < b) ? sd.u : b
    xy = xy0
    for _ in 1:100
        n  = exp(u)
        xy = u_xy(t, n, sd)
        U_COUNT[3] += 1
        G  = log(sum(u_species(t, xy.x, xy.y)))     # ln S
        if abs(G) <= 1e-12 || xy.hfloor
            sd.u = u
            return (n = n, xy = xy)
        end
        if G > 0
            a = u
            if u == b                               # still oversaturated at nmax: that end
                sd.u = u
                return (n = n, xy = xy)
            end
        else
            b = u
        end
        # dG/du = -G_z·F_z⁻¹·F_u from one dual pass in (x, y, u)
        h = w -> begin
            xs = t.nC > 0 ? w[1] : oftype(w[1], -Inf)
            ys = t.nH > 0 ? w[2] : oftype(w[2], -Inf)
            F  = u_F(t, SVector(w[1], w[2]), exp(w[3]))
            SVector(F[1], F[2], log(sum(u_species(t, xs, ys))) + 0*w[3])
        end
        _, J3 = u_valjac(h, SVector(t.nC > 0 ? xy.x : 0.0, t.nH > 0 ? xy.y : 0.0, u))
        Fz = SMatrix{2,2}(J3[1,1], J3[2,1], J3[1,2], J3[2,2])   # closures in (x, y)
        Fu = SVector(J3[1,3], J3[2,3])                          # closures in u
        Gz = SVector(J3[3,1], J3[3,2])                          # ln S in (x, y)
        dG = -(Gz' * (Fz \ Fu))
        un = u - G/dG                               # Newton step in u
        if isfinite(un) && a < un < b
            u = un
        else
            u = 0.5*(a + b)                         # outside the bracket: bisect
            U_COUNT[4] += 1
        end
        b - a <= 1e-14*max(1.0, abs(u)) && break
    end
    sd.u = u
    return (n = exp(u), xy = xy)
end

# Everything UTrial holds, from the trial fO2, the current φ and the cell
function u_trial(f, T, P, Pbar, fO₂, K, φ, nC, nH, nOx, Mf, Ms, fl_ok, sol_on, la, lC, φw0, ylo, yhi,
                 melt_H2, melt_H2O, h2_comp, XH2max, XCO2max, h2o_law, wGs)
    K1, K2, K3, K4 = K                              # reaction constants of the four C-O-H equilibria
    sf  = sqrt(fO₂)
    c   = (K3*sf/(φ[2]*Pbar),                       # H2O
           K1*fO₂/(φ[3]*Pbar),                      # CO2
           K4/(φ[1]*Pbar),                          # CH4
           K2*sf/(φ[5]*Pbar),                       # CO
           1/(φ[4]*Pbar))                           # H2
    caw = K3*sf/((sol_on ? φw0 : φ[2])*Pbar)        # a_H2O = f_H2O / f_H2O(pure water)
    x_cap = log(min(1.0, 10.0^(logfO2_EMOG(T, Pbar) - f)))   # ln aC where solid carbon saturates
    # Sun & Yao at a = 1 with P held at their 2 GPa ceiling; the activity then enters as a power law
    Tm  = max(T, SY_TMIN)
    Pc  = min(Pbar, SY_PMAX)
    a   = SY_H2O
    kSm = a[1]*log(Pc) + a[2]*log(Tm) + a[3]*Pc^2/Tm + (a[4]*sqrt(Pc) + a[5])*MELT_XNA2O/Tm
    kSo = a[6]*log(Pc) + a[7]*log(Tm) + (a[8] + a[9]*MELT_XAL2O3 + a[10]*MELT_XNA2O)*Pc/Tm
    kKatz = min(12*P^0.6 + P, 12*10.0^0.6 + 10.0)   # Katz (2003) with the Rozel (2018) cutoff
    return UTrial(fl_ok, c, caw, nC, nH, nOx, Mf, Ms, invert_carbon_buffer_fO2(f, P, T),
                  min(max(XCO2max, 1e-6), 0.95), x_cap, wGs*log(10.0),
                  melt_H2, h2_melt_K(T, Pbar; comp=h2_comp), XH2max,
                  melt_H2O, h2o_law === :katz, kSm, kSo, kKatz,
                  sol_on, la, lC, x_cap - 700.0, x_cap, ylo, yhi)
end

# Fluid, carbon and hydrogen of every phase at one trial fO2: fluid amount and potentials at frozen φ,
#  then φ from the fluid EoS until it settles. Signature and output as before the Newton rewrite.
function coh_state_u(f, T, P, Mh, Mc, nOxideTotal, TOₑₓ; nit_min=2, nit_max=30, id_tol=1e-8,
                     wEMOG=0.05, wCO=0.05, XCO2max=0.5,
                     melt_H2=true, melt_H2O=true, h2_comp=:basalt, XH2max=0.3,
                     Mf=0.0, molfracXB=nothing, h2o_law=:sunyao, fluid_Pmax=50.0, φ0=nothing, sol=nothing, wGs=0.01,
                     seed=nothing)
    Pbar, fO₂ = 1e4*P, 10.0^f
    K  = reaction_Ks(T, Pbar)
    nH, nC = 2*Mh/mm.H2O, Mc/mm.C                   # budgets: mol H, mol C
    fl_ok = P <= fluid_Pmax
    φ  = (φ0 === nothing || !fl_ok) ? ones(7) : Vector{Float64}(φ0)   # φ from the previous trial, or ideal
    sd = seed === nothing ? RootSeed() : seed       # roots carried across passes and trials of this cell
    sol_on = sol !== nothing
    Ms  = sol_on ? Float64(sol.Ms) : 0.0
    la  = sol_on ? Vector{Float64}(sol.C.la) : Float64[]
    lC  = sol_on ? Vector{Float64}(sol.C.lC) : Float64[]
    φw0 = sol_on ? Float64(sol.φw0) : 1.0
    y_ref = fl_ok ? log(Pbar/(1/φ[4] + K[3]*sqrt(fO₂)/φ[2])) : 0.0   # ln fH2 of a carbon-free H2O-H2 fluid
    ylo, yhi = y_ref - 140.0, y_ref + 20.0          # hydrogen search range
    nothing_there = nC <= 0 && nH <= 0
    mk() = u_trial(f, T, P, Pbar, fO₂, K, φ, nC, nH, Float64(nOxideTotal), Float64(Mf), Ms, fl_ok, sol_on,
                   la, lC, φw0, ylo, yhi, melt_H2, melt_H2O, h2_comp, XH2max, XCO2max, h2o_law, wGs)
    t  = mk()
    n, x, y, hfloor, it, vol_fluid = 0.0, -Inf, -Inf, false, 0, NaN
    vx = zeros(7)                                   # fluid composition in ZD09 order: CH4 H2O CO2 H2 CO O2 C2H6
    for outer it in 1:(nothing_there ? 0 : (fl_ok ? max(1, nit_max) : 1))
        r = u_solve_n(t, sd)                        # fluid amount and potentials at frozen φ
        n, x, y, hfloor = r.n, r.xy.x, r.xy.y, r.xy.hfloor
        fl_ok || break                              # pure water: φ does not apply
        sp = u_species(t, x, y)
        S  = sum(sp)
        vx .= 0.0
        vx[1], vx[2], vx[3], vx[4], vx[5] = sp[3], sp[1], sp[2], sp[5], sp[4]
        S > 0 && (vx ./= S)                         # the fluid, or the incipient one
        vol_fluid, fug = fluid_fugacity(vx, Pbar, T)
        dφ = 0.0                                    # largest relative change of φ this pass
        for i in 1:7
            vx[i] > 0 || continue
            nφ = fug[i]/(vx[i]*Pbar)                # new fugacity coefficient
            (isfinite(nφ) && nφ > 0) || continue
            dφ = max(dφ, abs(nφ - φ[i])/φ[i])
            φ[i] = it > 4 ? sqrt(φ[i]*nφ) : nφ      # geometric damping after 4 passes: the loop can cycle
        end
        it >= nit_min && dφ < id_tol && break
        t = mk()                                    # the next pass at the new φ
    end
    # The state at the potentials actually solved (φ of the last solve, not the refreshed one)
    sp  = nothing_there ? SVector(0.0, 0.0, 0.0, 0.0, 0.0) : u_species(t, x, y)
    S   = sum(sp)
    hs  = nothing_there ? (nH2 = 0.0, nHw = 0.0, nHs = 0.0, aw = 0.0, fH2 = 0.0) : u_hsinks(t, y, sp)
    XH2O, XCO2, XCH4, XCO, XH2 = S > 0 ? (sp[1]/S, sp[2]/S, sp[3]/S, sp[4]/S, sp[5]/S) : (0.0, 0.0, 0.0, 0.0, 0.0)
    aC  = isfinite(x) ? exp(x) : 0.0
    fE  = logfO2_EMOG(T, Pbar)
    aC_cap   = exp(t.x_cap)
    molH_fl  = n*(2sp[1] + 4sp[3] + 2sp[5])         # fluid hydrogen [mol H]
    nC_fluid = n*(sp[2] + sp[3] + sp[4])            # fluid carbon [mol]
    XCO2_melt = nOxideTotal > 0 ? softcap(aC*t.Bmelt, t.XCO2cap) : 0.0
    nC_melt  = XCO2_melt*nOxideTotal/(1 - XCO2_melt)   # melt carbon [mol]
    nC_solid = nC > 0 ? nC*exp(min(x - t.x_cap, 0.0)/t.wC) : 0.0   # solid carbon [mol]
    x_carb   = 1/(1 + exp(-(f - fE)/wEMOG))         # carbonate share of it across EMOG
    nC_carb  = x_carb*nC_solid
    nC_gph   = nC_solid - nC_carb
    carbO    = 2*nC_carb*mm.O                       # solid carbonate: 2 O per carbon [g]
    netO  = n > 0 ? n*(2sp[2] + sp[4] - 2sp[3] - sp[5])*mm.O : 0.0   # fluid excess O [g], signed
    fₑₓ   = isfinite(netO) ? clamp(netO, -1e3, 1e3) : NaN
    O_up  = n > 0 ? n*(2sp[2] + sp[4])*mm.O : 0.0   # O the fluid takes up as CO2 and CO
    O_rel = n > 0 ? n*(2sp[3] + sp[5])*mm.O : 0.0   # O it releases as CH4 and H2
    if !fl_ok && n > 0                              # Darcy density above Pmax: the pure-water ZD09 root
        vw = zeros(7); vw[2] = 1.0
        vpw, _ = fluid_fugacity(vw, Pbar, T)
        isfinite(vpw) && (vol_fluid = vpw)
    end
    nH2_melt, nH_H2O_melt, fH2 = hs.nH2, hs.nHw, hs.fH2
    PH2O_melt = hs.aw*Pbar
    meltH2ex  = -nH2_melt*mm.O                      # melt H2 hands back oxygen [g]
    XH2_melt  = (melt_H2 && nOxideTotal > 0) ? softcap(t.KH2*fH2, XH2max) : 0.0
    SH2Om, SOH = (melt_H2O && Mf > 0) ? melt_h2o_wt(P, max(T, SY_TMIN), PH2O_melt/Pbar; law=h2o_law) : (0.0, 0.0)
    nH2O_melt  = nH_H2O_melt/2
    XH2O_melt  = nOxideTotal > 0 ? nH2O_melt/(nH2O_melt + nOxideTotal) : 0.0
    # Hydrogen floor without melt: the solid takes the closure exactly, at an activity below the table
    #  (reported as its lowest node). With melt the split is undetermined there: h_floor_melt fails the cell.
    nHs, a_solid = hs.nHs, sol_on ? hs.aw : 0.0
    h_floor_melt = hfloor && Mf > 0
    if hfloor && !h_floor_melt
        nHs     = max(nH - molH_fl - 2nH2_melt - nH_H2O_melt, 0.0)
        a_solid = exp(la[end])
    end
    rH = nH > 0 ? abs(molH_fl + 2nH2_melt + nH_H2O_melt + nHs - nH)/nH : 0.0   # hydrogen closure error
    rC = nC > 0 ? abs(nC_melt + nC_fluid + nC_solid - nC)/nC : 0.0              # carbon closure error
    ctol = 1e-6*max(nC, eps())
    co_fluid = n > 0 && molH_fl <= 1e-9*max(nH, 1e-30) + 1e-300                  # a fluid with no hydrogen
    return (; aC, aC_cap, fE, vol_fluid, it_ideality = it, res_ideality = max(rC, rH),
              graphite = nC_gph > ctol, carbonate = nC_carb > ctol, valid = true,
              XH2O, XCO2, XCH4, XCO, XH2, n_fluid = n, nC_fluid,
              XCO2_melt, nC_melt, nC_solid, nC_carb, nC_gph, carbO, fₑₓ, O_up, O_rel,
              nH2_melt, molH_fl, molH_tot = nH, fH2, XH2_melt, meltH2ex,
              nH2O_melt, nH_H2O_melt, XH2O_melt, SH2Om, SOH, PH2O_melt, lam = min(1.0, S),
              nH_solid = nHs, a_solid, S_fluid = S, h_floor = hfloor, h_floor_melt,
              fluid_absent = n <= 0, fluid_restricted = !fl_ok, co_fluid,
              φ_seed = fl_ok && !nothing_there ? copy(φ) : nothing)
end
