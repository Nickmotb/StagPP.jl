# =============================================================================================
# TOex.jl — excess oxygen (Oₑₓ) partitioned between solid, melt and C-O-H fluid at fixed P, T.
#           Entry point `partition_Oₑₓ`. Numbers, benchmarks, open issues: TOex_bench.txt.
# =============================================================================================
    # FOUR SINKS, ONE SHARED fO2. Each holds more oxygen as fO2 rises, so every demand is monotone.
    #   solid Fe₂O₃   Fe³⁺ in the residue      MAGEMin sb24 (Stixrude & Lithgow-Bertelloni 2024)
    #   melt  Fe₂O₃   Fe³⁺ in the melt         Hirschmann (2022) oxybarometer
    #   melt  CO₂     C dissolved in the melt  EMOG(D), Stagno & Frost (2010)-style
    #   solid carb.   carbonate, 2 O per C     EMOG(D) from Holland & Powell (1998)
    #   C-O-H fluid   H₂O/CO₂/CH₄/CO/H₂       HP98 Kᵢ(T) + Zhang & Duan (2009) fugacities
    # METHOD. fO2 is the one unknown, mass balance (Σ demands = 1) the one equation.
    #   Outer bisection on fO2, inner on aC; φ frozen inside the aC solve and refreshed once per
    #   outer step, so the EOS runs once per step rather than once per aC trial.
    # fO2 equality is exact BY CONSTRUCTION: an unsaturated sink IS the loading reading that fO2.
    # Monotone demands ⇒ unique root ⇒ plain bisection. `classify_regime` only labels the answer.
    # SIGN. Oₑₓ counted against a REDUCED reference (Fe as FeO, C as graphite, H as H₂O).
    #   Negative share = SOURCE: the fluid when H₂O → CH₄/H₂, the solid below IW when Fe²⁺ → Fe⁰.
    # UNITS. P GPa (bar internally); T K; TOex/TC/H mass fractions 0-1; shares are fractions of
    #   TOₑₓ; every fO2 is log₁₀.
    # LIMITS (evidence in TOex_bench.txt): ~1.9 log unit absolute fO2 offset, largely MAGEMin never
    #   applying sb24's Landau term (Section 4c); no precipitating Fe³⁺ phase, so an over-large
    #   budget reports `capacity_exceeded`; melt CO₂ capped at `XCO2max`; carbonate is pure MgCO₃;
    #   C₂H₆ excluded; ΔFMQ/ΔIW are offsets WITHIN sb24; Zhang & Duan valid 673-2573 K, ≤10 GPa.
# =============================================================================================

# ------------------------------------------------------------
# Section 1 — bulk composition & mass bookkeeping
# ------------------------------------------------------------

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

    Φ = 2*mm.O*molMf/TOₑₓ   # converts a melt XCO₂ into a fraction of the TOₑₓ budget

    return (; Mf, Ms, Mc, Mh, TOₑₓ, molfracXB, nB_per_mass, molMf, maxXCO₂_raw, Φ)
end

# ------------------------------------------------------------
# Section 2 — solid branch: MAGEMin (Stixrude & Lithgow-Bertelloni 2024)
# ------------------------------------------------------------

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
function solid_fO2_interpolant(P, T, X, Xox, Ms, TOₑₓ, s_stoich; nr=25, Rmax=0.2, nrm=15, R0=2e-3, data=nothing, verbose=false)
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

    flag = isnothing(data)
    flag && (data = Initialize_MAGEMin("sb24", verbose=false))
    out  = multi_point_minimization(10P*ones(ntot), k2c(T)*ones(ntot), data, X=Xlist, Xoxides=Xred, name_solvus=true, sys_in="wt", progressbar=false)
    flag && Finalize_MAGEMin(data)

    grid = sOwt .* (Ms/TOₑₓ)   # re-express the x-axis as a fraction of the TOₑₓ budget
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
    return sfO2, grid[1], grid[end], monotonic, -oxb(-2.0)[iFe]*(mm.O/mm.Fe)*(Ms/TOₑₓ), qfm,
           s -> max(0.0, mitp(clamp(s, grid[1], grid[end])))
end

# ------------------------------------------------------------
# Section 3 — melt branch: Hirschmann (2022) oxybarometer
# ------------------------------------------------------------

# log₁₀ fO2 of a melt holding a fraction mOₑₓ of the TOₑₓ budget as Fe³⁺.
function melt_fO2(mOₑₓ, T, TOₑₓ, Mf, molfracXB, nB_per_mass, IDV)
    a=0.1917; b=-1.961; c=4158.1; ΔCₚ=33.25; T₀=1673.15                             # Hirschmann (2022)
    y1=-520.46; y3=494.39; y4=1838.34; y5=2888.48; y8=-1245.09; y9=-1156.86         # regression terms

    nOxideTotal = Mf*nB_per_mass           # moles of the melt's own unoxidized oxides
    nOₑₓ_added  = TOₑₓ*mOₑₓ/mm.O           # moles of excess O added to the melt
    θ           = nOxideTotal + nOₑₓ_added # renormalization after adding that oxygen

    XO     = nOₑₓ_added/θ
    XSiO2  = molfracXB[1]*nOxideTotal/θ
    XMgO   = molfracXB[2]*nOxideTotal/θ
    XFeO   = molfracXB[3]*nOxideTotal/θ
    XCaO   = molfracXB[4]*nOxideTotal/θ
    XAl2O3 = molfracXB[5]*nOxideTotal/θ
    XNa2O  = molfracXB[6]*nOxideTotal/θ

    (XO >= 0.5XFeO) && return oftype(mOₑₓ, 14.0)   # more O than there is Fe to host it as Fe₂O₃

    _T, _ln10 = 1/T, 1/log(10)
    return (log10(XO/(XFeO - 2XO)) - b - c*_T + (ΔCₚ*_ln10/R)*(1 - T₀*_T - log(T/T₀)) + (IDV/(1e-3R))*_T*_ln10
              - _T*(y1*XSiO2 + y3*XMgO + y4*XCaO + y5*XNa2O + y8*XSiO2*XAl2O3 + y9*XSiO2*XMgO)) / a
end

# ------------------------------------------------------------
# Section 4 — melt CO₂ branch: EMOG(D) buffer (Stagno & Frost, 2010-style)
# ------------------------------------------------------------

# log₁₀ fO2 at which graphite/diamond coexists with a melt of the given XCO₂.
function carbon_buffer_fO2(XCO₂, P, T)
    @assert XCO₂ > 0.0 "XCO₂ must be strictly positive."
    return 5.44 - 21380/T + 0.078*(1e4P - 1)/T + log10(XCO₂)
end
# Closed-form inverse: the XCO₂ that puts the buffer at a given fO2.
invert_carbon_buffer_fO2(f, P, T) = 10.0^(f - 5.44 + 21380/T - 0.078*(1e4P - 1)/T)


# ------------------------------------------------------------
# Section 4b — C-O-H fluid: Holland & Powell (1998) Kᵢ(T), Zhang & Duan (2009) fugacities
# ------------------------------------------------------------

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

function zd_bisect(f, lo, hi; iters=200)
    flo, fhi = f(lo), f(hi)
    flo*fhi > 0 && return NaN
    for _ in 1:iters
        mid = 0.5*(lo+hi); fm = f(mid)
        sign(fm)==sign(flo) ? (lo=mid; flo=fm) : (hi=mid)
    end
    return 0.5*(lo+hi)
end

# The EOS has several real roots; try brackets widest-first to land on the physical fluid-like one.
function zd_findroot(f)
    for (lo,hi) in ((1.0,1000.0), (0.01,10.0), (0.001,1.0))
        r = zd_bisect(f, lo, hi)
        (isfinite(r) && abs(f(r)) < 1e-6) && return r
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
    Vm = zd_findroot(V -> zd_EOS(V, Tm, Pm))
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

# One speciation step at fixed φᵢ and carbon activity aC: a quadratic in XH2.
# `valid=false` ⇒ the carbon species alone would exceed the fluid; `_solve_aC` fixes that by lowering aC.
function _coh_reaction_step(fO₂, Pbar, K1, K2, K3, K4, φCH4, φH2O, φCO2, φH2, φCO, aC)
    # Rewrite species fugacity in terms of ideality coefficient, mol fraction and pressure. Isoalte then fractions:
    XCO2 = aC*K1*fO₂/(φCO2*Pbar)
    XCO  = aC*K2*√fO₂/(φCO*Pbar)
    # molar proportions within the fluid mass sum up to 1. Here we express the available room for H-bearing
    # species by subtracting XCO2 and XCO from 1.
    aq2  = aC*K4*φH2^2*Pbar/φCH4
    aq1  = K3*φH2*√fO₂/φH2O + 1
    c_q  = 1 - XCO2 - XCO

    # XH2O and XCH4 reduce to functions of XH2 alone (the pressure terms cancel), leaving
    #       aq2·XH2² + aq1·XH2 - c_q = 0

    # Equation descriminant
    disc = aq1^2 + 4aq2*c_q
    # Return sentinel composition if impossible or no space left for H-bearing species
    if disc < 0 || c_q <= 0
        return false, 0.5, 0.5, 0.0, 0.0, 0.0
    end
    # Equation real root.
    XH2 = 2c_q/(aq1 + √disc)
    # Again, return dummy if non-sensical
    (XH2 < 0) && return false, 0.5, 0.5, 0.0, 0.0, 0.0
    # Retrieve XH2O and XCH4 fractions based off root-found XH2
    XH2O = K3*φH2*√fO₂/φH2O * XH2
    XCH4 = aq2*XH2^2
    return true, XH2O, XCO2, XCH4, XCO, XH2
end

# Moles of carbon that melt + fluid together demand at a given carbon activity.
function _carbon_demand(aC, fO₂, Pbar, K1, K2, K3, K4, φ, molH, nOxideTotal, Bmelt, XCO2max)
    Xm = softcap(aC*Bmelt, XCO2max)
    nC_melt = Xm*nOxideTotal/(1 - Xm)

    # If no free water -> demand is set by carbon in the melt
    molH <= 0 && return nC_melt

    ok, XH2O, XCO2, XCH4, XCO, XH2 = _coh_reaction_step(fO₂, Pbar, K1, K2, K3, K4, φ[1], φ[2], φ[3], φ[4], φ[5], aC)
    ok || return Inf
    nH = 2XH2O + 4XCH4 + 2XH2
    nC_fluid = nH > 0 ? molH*(XCO2 + XCH4 + XCO)/nH : 0.0
    return nC_melt + nC_fluid
end

# Carbon activity satisfying carbon mass balance, bisected in log₁₀ aC and capped by whichever solid
# host is stable: pinned at the cap ⇒ that solid is present; below it ⇒ no solid carbon at all.
function _solve_aC(fO₂, Pbar, K1, K2, K3, K4, φ, molH, nC_total, nOxideTotal, Bmelt, aC_cap, XCO2max; iters=120)
    nC_total <= 0 && return 0.0                        # no carbon at all ⇒ no carbon species
    d1 = _carbon_demand(aC_cap, fO₂, Pbar, K1, K2, K3, K4, φ, molH, nOxideTotal, Bmelt, XCO2max)
    (isfinite(d1) && d1 <= nC_total) && return aC_cap   # a solid carbon phase is present
    lo, hi = -300.0, log10(aC_cap)                     # log₁₀ aC bracket
    for _ in 1:iters
        mid = 0.5*(lo + hi)
        _carbon_demand(10.0^mid, fO₂, Pbar, K1, K2, K3, K4, φ, molH, nOxideTotal, Bmelt, XCO2max) < nC_total ?
            (lo = mid) : (hi = mid)
    end
    return 10.0^(0.5*(lo + hi))
end

# Fluid, carbon activity and oxygen draw solved together at one fO2. φ is frozen while aC is
# bisected, then refreshed from the result, so the EOS runs once per outer iteration, not per trial.
function coh_state(f, T, P, Mh, Mc, nOxideTotal, TOₑₓ; outer_iters=6, wEMOG=0.25, XCO2max=0.3)
    Pbar  = 1e4*P                  # GPa → bar, matching the 1-bar standard state of the data
    fO₂   = 10.0^f
    K1, K2, K3, K4 = reaction_Ks(T, Pbar)
    molH     = 2*Mh/mm.H2O
    nC_total = Mc/mm.C
    Bmelt    = invert_carbon_buffer_fO2(f, P, T)       # melt XCO₂ at aC=1

    fE     = logfO2_EMOG(T, Pbar)   # EMOG(D): graphite/diamond below it, carbonate above

    # aC ceiling from the stable solid host: 1 for graphite/diamond, 10^(f_EMOG - f) for carbonate.
    # `_solve_aC` prices demand at this cap — under budget ⇒ solid C takes the rest, over ⇒ aC drops.
    aC_cap = min(1.0, 10.0^(fE - f))

    φ  = ones(7) # Ideal gas as starting point
    aC = aC_cap
    ok, XH2O, XCO2, XCH4, XCO, XH2 = true, 0.5, 0.5, 0.0, 0.0, 0.0
    for _ in 1:outer_iters
        # Solve for carbon mass balance -> i.e pin-pointing aC
        aC = _solve_aC(fO₂, Pbar, K1, K2, K3, K4, φ, molH, nC_total, nOxideTotal, Bmelt, aC_cap, XCO2max)
        # Resolve C-O-H fluid given P,T,fO2 and aC
        ok, XH2O, XCO2, XCH4, XCO, XH2 = _coh_reaction_step(fO₂, Pbar, K1, K2, K3, K4, φ[1], φ[2], φ[3], φ[4], φ[5], aC)
        vx = zeros(7); vx[1]=XCH4; vx[2]=XH2O; vx[3]=XCO2; vx[4]=XH2; vx[5]=XCO
        s = sum(vx); (s > 0) && (vx ./= s)
        # C-O-H species ideality coefficients update
        _, fug = zd_fugacity(vx, Pbar, T)
        for i in 1:7
            newφ = fug[i]/(vx[i]*Pbar)
            (isfinite(newφ) && newφ > 0) && (φ[i] = newφ)
        end
    end

    # Extract absolute mols of fluid through molH/(molH/molfl) based off speciation
    nH       = 2XH2O + 4XCH4 + 2XH2
    n_fluid  = (Mh > 0 && nH > 0) ? molH/nH : 0.0

    # Extract absolute amount of carbon in the fluid based off speciation
    nC_fluid = n_fluid*(XCO2 + XCH4 + XCO)

    # Extract absolute amount of carbon in the melt as CO2
    XCO2_melt = softcap(aC*Bmelt, XCO2max)
    nC_melt   = XCO2_melt*nOxideTotal/(1 - XCO2_melt)

    # Leftover carbon sits in a solid phase, but only while aC is pinned at its ceiling.
    solid_present = nC_total > 0 && aC >= aC_cap*(1 - 1e-9)
    nC_solid = solid_present ? max(0.0, nC_total - nC_melt - nC_fluid) : 0.0

    # Which solid host: graphite/diamond below EMOG(D), carbonate above; ramped over `wEMOG`, not switched.
    x_carb   = 1/(1 + exp(-(f - fE)/wEMOG))
    nC_carb  = x_carb*nC_solid                 # carbonate: MgCO₃ vs MgO+C ⇒ 2 excess O per carbon
    nC_gph   = nC_solid - nC_carb              # graphite/diamond: carries no excess oxygen
    carbO    = 2*nC_carb*mm.O/TOₑₓ

    # Net draw vs the reduced reference (H₂O itself neutral): CO₂/CO consume O (+2/+1), CH₄/H₂
    # release it (-2/-1). Negative is physical — a reduced fluid donates oxygen back to the rock.
    netO = n_fluid*(2XCO2 + XCO - 2XCH4 - XH2)*mm.O / TOₑₓ
    fₑₓ  = clamp(netO, -1e3, 1e3)   # |·|>1 already exceeds the whole budget; bound for the bracket

    # The two GROSS flows behind that net. A host model needs them: a reducing fluid GROWS the Oₑₓ
    # inventory, and that source has to be booked against the hydrogen reservoir, not re-derived.
    O_up  = n_fluid*(2XCO2 + XCO)*mm.O / TOₑₓ    # consumed oxidising C to CO₂/CO
    O_rel = n_fluid*(2XCH4 + XH2)*mm.O / TOₑₓ    # released reducing H₂O to CH₄/H₂

    # Thresholded: the logistic ramp's exponential tails would make a bare `> 0` report phantom phases.
    ctol = 1e-6*max(nC_total, eps())
    return (; aC, aC_cap, fE, graphite = nC_gph > ctol, carbonate = nC_carb > ctol, valid = ok,
              XH2O, XCO2, XCH4, XCO, XH2, n_fluid, nC_fluid,
              XCO2_melt, nC_melt, nC_solid, nC_carb, nC_gph, carbO, fₑₓ, O_up, O_rel)
end

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

# ------------------------------------------------------------
# Section 4c — reference buffers: FMQ and IW, both out of sb24
# ------------------------------------------------------------
# An offset means nothing unless buffer and fO2 share a dataset: Frost (1991) put ΔFMQ up to +4.0
#   off and metal-BEARING rock ABOVE its own IW, which is impossible. Numbers in TOex_bench.txt.
# sb24's IW sits ~1.4 log units high because MAGEMin never applies the Landau term (`cme` = 0), so
#   wüstite is short R·ln5 per Fe²⁺. It cancels in ΔIW — adding it to IW alone breaks metal.
# FMQ is free from MAGEMin (QFM = fO2 - dQFM); `logfO2_FMQ` covers the no-solid case.
# IW is built here because MAGEMin's `iw` phase calls an endmember `iron` sb24 lacks, and SEGFAULTS.
# Parameters below are SLB (2024) Table C1 as printed: ≤0.0084 log units on FMQ, ≤0.0012 on IW.

const T0_SLB = 300.0        # SLB reference temperature; F₀ and the thermal terms are relative to it

# SLB (2024) Table C1: F0 J/mol | n atoms/fu | V0 J/bar/mol | K0 bar | Kp | θ0 K | γ0 | q0 | fe.
# NOTE `wu` is Fe₄O₄, not FeO — SLB's B1 phases are written on four oxygens (`pe` is Mg₄O₄ too).
const SLB_EM = (
    wu   = (F0=-962e3,  n=8.0, V0=4.902, K0=161e4, Kp=4.0,  θ0=427.0, γ0= 1.45, q0=1.5, fe=0),
    wuls = (F0=-609e3,  n=8.0, V0=4.340, K0=200e4, Kp=4.0,  θ0=493.0, γ0= 1.45, q0=1.5, fe=0),
    fea  = (F0=    0.0, n=1.0, V0=0.709, K0=163e4, Kp=6.0,  θ0=398.0, γ0= 1.67, q0=0.9, fe=1),
    feg  = (F0=  4.7e3, n=1.0, V0=0.699, K0=165e4, Kp=6.0,  θ0=285.0, γ0= 1.85, q0=1.0, fe=3),
    fee  = (F0=  3.4e3, n=1.0, V0=0.677, K0=165e4, Kp=5.0,  θ0=380.0, γ0= 2.13, q0=1.0, fe=2),
    qtz  = (F0= -858e3, n=3.0, V0=2.242, K0= 61e4, Kp=19.8, θ0=883.0, γ0=-0.04, q0=1.0, fe=0),
    fa   = (F0=-1372e3, n=7.0, V0=4.629, K0=136e4, Kp=4.9,  θ0=619.0, γ0= 1.08, q0=2.9, fe=0),
    mag  = (F0=-1016e3, n=7.0, V0=4.453, K0=184e4, Kp=5.2,  θ0=529.0, γ0= 1.36, q0=0.4, fe=0),
)

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

    # Magnetic ordering, α-iron only: Hillert & Jarl (1978), which SLB adopt for bcc in place of
    # Landau (A5-A9, Table C3). `p` = 0.4 is the standard bcc value, defined in A7 but not given.
    if d.fe == 1
        SD, TC, pf = 9.46, 1043.01, 0.4
        D = 518/1125 + (11692/15975)*(1/pf - 1)
        τ = T/TC
        G += τ <= 1 ?
             T*SD/D*(1 - (79/τ/140/pf + (474/497)*(1/pf - 1))*(τ^3/6 + τ^9/135 + τ^15/600)) :
            -T*SD/D*(τ^-5/10 + τ^-15/315 + τ^-25/1500)
    end
    return G
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

# hp02ver FMQ (3Fe₂SiO₄ + O₂ = 2Fe₃O₄ + 3SiO₂), an independent third opinion: a near-constant +3.4
# to +3.8 above sb24's, i.e. a dataset offset. No hp98 IW twin — HP carry no solid FeO endmember.
const FMQ_FA = (G0=-1378945.0, S=151.0, c1=201.1, c2=0.01733,   c3=-1960600.0, c5=-900.9,
                 V0=4.631,  b1=0.505e-4, b5=-0.505e-3, b6=1330000.0, b7=-199.5, b8=4.0)
const FMQ_QZ = (G0=-856424.6,  S=41.5,  c1=110.7, c2=-0.005189, c3=0.0,        c5=-1128.3,
                 V0=2.2688, b1=0.65e-5,  b5=-0.65e-4,  b6=750000.0,  b7=-112.5, b8=4.0)
const FMQ_MT = (G0=-1012288.0, S=146.1, c1=262.5, c2=-0.007204, c3=-1926200.0, c5=-1655.7,
                 V0=4.452,  b1=0.696e-4, b5=-0.696e-3, b6=1850000.0, b7=-277.5, b8=4.0)
function logfO2_FMQ_hp98(T, Pbar)
    dG0 = 2hp98_gibbs(T,FMQ_MT) + 3hp98_gibbs(T,FMQ_QZ) - 3hp98_gibbs(T,FMQ_FA) - hp98_gibbs(T,HP98.O2)
    dGP = 2murnaghan_dG(T,Pbar,FMQ_MT) + 3murnaghan_dG(T,Pbar,FMQ_QZ) - 3murnaghan_dG(T,Pbar,FMQ_FA)
    return (dG0 + dGP)/(R*T*log(10))
end

# ------------------------------------------------------------
# Section 5 — smooth saturation and finiteness helpers
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# Section 6 — the residual (1 unknown: the shared log₁₀ fO2)
# ------------------------------------------------------------

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

# A trial fO2 → the full state: each sink's demand, its own fO2 reading, the total. Used for both
# solving and reporting, so the reported state is exactly the one the residual saw.
function _coh_derive(f, ctx)
    (; sfO2, T, P, IDV, TOₑₓ, Mf, Mc, molfracXB, nB_per_mass, Φ, Mh, nOxideTotal, maxmOₑₓ, maxsOₑₓ, minsOₑₓ, lowclip) = ctx

    # The solid's bracket straddles zero (below IW it hands oxygen back); the melt has no metal branch.
    sOₑₓ = _invert_demand(s -> sfO2(s), f, -minsOₑₓ + lowclip, maxsOₑₓ - lowclip)
    mOₑₓ = _invert_demand(m -> melt_fO2(m, T, TOₑₓ, Mf, molfracXB, nB_per_mass, IDV), f,
                          lowclip, maxmOₑₓ - lowclip)

    # One joint solve: aC enforces carbon mass balance exactly, so melt XCO₂ needs no ceiling of its own.
    cs      = coh_state(f, T, P, Mh, Mc, nOxideTotal, TOₑₓ)
    XCO₂    = cs.XCO2_melt
    carbonO = Φ*XCO₂/(1 - XCO₂)

    fₑₓ   = nancheck("coh_state.fₑₓ", cs.fₑₓ, f, Mh)
    carbO = nancheck("coh_state.carbO", cs.carbO, f, Mc)   # solid carbonate: 2 O per carbon
    fₛ    = nancheck("sfO2(sOₑₓ)", sfO2(sOₑₓ), sOₑₓ)
    fₗ    = nancheck("melt_fO2", melt_fO2(mOₑₓ, T, TOₑₓ, Mf, molfracXB, nB_per_mass, IDV), mOₑₓ)

    return (; f, sOₑₓ, mOₑₓ, XCO₂, carbonO, carbO, fₑₓ, fₛ, fₗ, cs,
              total = sOₑₓ + mOₑₓ + carbonO + carbO + fₑₓ)
end

# The one equation: total demand must equal the budget.
unified_residual(f, ctx) = _coh_derive(f, ctx).total - 1

# ------------------------------------------------------------
# Section 6b — regime diagnosis (reporting only)
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# Section 7 — bracketed bisection solver
# ------------------------------------------------------------

# Bisection suits a 1-D monotone bracketed residual: no derivatives, no guess, no divergence.
# Convergence is judged on the RESIDUAL alone — a collapsed bracket with |F| large means no root exists.
function bisect_solve(F::Function, lo::Real, hi::Real; maxiter=200, atol=1e-10, xtol=1e-14, debug=false)
    flo, fhi = F(lo), F(hi)
    # No sign change: the root is outside the bracket, so return the closer endpoint, unconverged.
    if sign(flo) == sign(fhi)
        best = abs(flo) <= abs(fhi) ? lo : hi
        return (converged=false, x=best, resnorm=abs(F(best)), iters=0)
    end
    x, fx = lo, flo
    it = 0
    while it < maxiter
        it += 1
        x  = 0.5*(lo + hi)
        fx = F(x)
        debug && println("  it=$it: x=$x F=$fx  bracket=[$lo, $hi] width=$(hi-lo)")
        (abs(fx) <= atol || (hi - lo) <= xtol*max(1.0, abs(x))) && break
        sign(fx) == sign(flo) ? (lo = x; flo = fx) : (hi = x)
    end
    return (converged = abs(fx) <= atol, x=x, resnorm=abs(fx), iters=it)
end

# ------------------------------------------------------------
# Section 8 — top-level driver
# ------------------------------------------------------------

# Partitions the Oₑₓ budget at (P,T); `out` returns a NamedTuple. `p`: 0 harzburgitic, 1 basaltic.
# `tol` = 1e-8 clears bisection's ~1e-11 floor and stays below the ~0.1-3.5 graphite cliff.
function partition_Oₑₓ(P::Real, T::Real, p::Real, ϕ::Real, TOex::Real, TC::Real, H::Real;
                           nr=100, Rmax=0.4, nrm=15, R0=2e-3, maxiter=500, tol=1e-8, verbose=true, data=nothing, out=false, debug=false)

    (; X, XB, Xox, SymXox, mmXox) = endmember_composition(p)
    Mt = 1e3   # reference cell mass (1 kg in grams, matching mm.* in g/mol)
    (; Mf, Ms, Mc, Mh, TOₑₓ, molfracXB, nB_per_mass, maxXCO₂_raw, Φ) = mass_budget(Mt, ϕ, TC, TOex, H, XB, mmXox)

    if TOₑₓ == 0.0
        verbose && println("Total Oₑₓ budget = 0.0 — nothing to partition.")
        return (converged=true, sOₑₓ=0.0, mOₑₓ=0.0, XCO₂=0.0, iters=0, resnorm=0.0)
    end

    IDV = nancheck("solve_∫ΔVdP", solve_∫ΔVdP([P-0.05P, P, P+0.05P], [T-0.05T, T, T+0.05T])[2,2,1], P, T)

    # Ceilings at Fe³⁺/ΣFe = 1: moles of Fe, half an O each, over the budget. NOT clipped at 100% of
    # TOₑₓ — that is the residual's job, and doing it twice would break fO2 equality.
    maxsOₑₓ = 0.5*(Ms*X[3]/mmXox[3])*mm.O/TOₑₓ

    nOxideTotal = Mf*nB_per_mass
    maxmOₑₓ = 0.5*molfracXB[3]*nOxideTotal*mm.O/TOₑₓ
    maxXCO₂ = clamp(maxXCO₂_raw, 0.0, 1.0)   # all carbon dissolved in the melt; reported, not imposed
    lowclip = 1e-10

    # Bounds are stoichiometric, not tabulation limits — `sfO2` extrapolates past both ends.
    # `tbl_lo`/`tbl_hi` are kept only so the report can flag an answer that relies on extrapolation.
    sfO2, tbl_lo, tbl_hi, tbl_mono, s_floor, qfm_out, smetal = solid_fO2_interpolant(P, T, X, Xox, Ms, TOₑₓ, maxsOₑₓ;
                                                            nr=nr, Rmax=Rmax, nrm=nrm, R0=R0, data=data, verbose=verbose)

    # The solid's FLOOR: every Fe²⁺ reduced to Fe⁰, handing all its iron-bound oxygen back. Taken
    # from the interpolant, not as 2·maxsOₑₓ. The melt has no counterpart — its floor is `lowclip`.
    minsOₑₓ = -s_floor

    ctx = (; sfO2, T, P, IDV, TOₑₓ, Mf, Mc, molfracXB, nB_per_mass, Φ, Mh, nOxideTotal, maxmOₑₓ, maxsOₑₓ, minsOₑₓ, lowclip)

    # Bracket far wider than physical: every demand saturates outside its range, so this only has
    # to straddle the root — it is not a constraint on the answer.
    res  = bisect_solve(f -> unified_residual(f, ctx), -120.0, 40.0; maxiter=maxiter, atol=tol, debug=debug)
    st   = _coh_derive(res.x, ctx)
    sOₑₓ, mOₑₓ, XCO₂ = st.sOₑₓ, st.mOₑₓ, st.XCO₂

    # Oxidation state at the loading actually solved, so it cannot drift from the answer.
    # Metal comes from MAGEMin, not from sOₑₓ, so ferric and metal can now both be non-zero at once.
    nFe_solid = Ms > 0 ? Ms*X[3]/mmXox[3] : 0.0          # moles Fe in the solid residue
    nFe_melt  = molfracXB[3]*nOxideTotal                 # moles Fe in the melt
    nO_solid  = TOₑₓ*st.sOₑₓ/mm.O
    FeR_melt  = nFe_melt  > 0 ? 2*(TOₑₓ*st.mOₑₓ/mm.O)/nFe_melt  : 0.0
    # Native iron from MAGEMin's assemblage, so it includes disproportionated metal.
    # `FeM_net` is the blind coordinate's own answer, a strict LOWER BOUND on FeM_solid.
    metal_wt   = Ms > 0 ? smetal(st.sOₑₓ) : 0.0                        # wt fraction of the solid
    nFe0_solid = metal_wt*Ms/mm.Fe                                     # moles of native iron
    FeM_solid  = nFe_solid > 0 ? nFe0_solid/nFe_solid : 0.0            # Fe⁰/ΣFe, total
    FeM_net    = nFe_solid > 0 ? max(0.0, -nO_solid/nFe_solid) : 0.0   # Fe⁰/ΣFe, net-oxygen only
    # With metal present the net oxygen alone no longer gives Fe³⁺: balance is nFe³⁺ = 2(nO + nFe⁰),
    # which returns the 2:1 ratio at nO = 0 and zero when metal is exactly paid for by released O.
    FeR_solid = nFe_solid > 0 ? max(0.0, 2(nO_solid + nFe0_solid)/nFe_solid) : 0.0
    FeR_net   = nFe_solid > 0 ? max(0.0, 2nO_solid/nFe_solid) : 0.0
    # The two routes, split by the oxygen budget since MAGEMin reports only a total: dissociation
    # frees exactly 1 O per Fe⁰ and that O IS -sOₑₓ, so the rest was made at constant oxygen.
    FeM_diss = min(FeM_net, FeM_solid)                                 # paid for by O leaving the solid
    FeM_disp = max(0.0, FeM_solid - FeM_diss)                          # made at constant oxygen

    eqfO2     = st.f       # the shared fO2 all unpinned sinks equilibrate to
    cohsink   = st.fₑₓ
    shortfall = 1 - st.total

    cs = st.cs
    # Both out of sb24, the dataset that produced `eqfO2`: FMQ off MAGEMin's own dQFM, IW from its
    # wüstite/iron endmembers. Section 4c explains why an empirical curve cannot go here.
    fmq, iw = reference_buffers(P, T; qfm=qfm_out)
    if verbose
        XH2O, XCO2f, XCH4, XCOf, XH2 = cs.XH2O, cs.XCO2, cs.XCH4, cs.XCO, cs.XH2
        rg = classify_regime(st, ctx, res)
        println("========================================================================================================")
        println("==== Solution (P=$(P) GPa | T=$(T) K | ϕ=$(1e2ϕ) % | TOex=$(1e2TOex) wt% | mix=$p | TC=$(1e2TC) wt% | H=$(1e2H) wt%) ====")
        println("========================================================================================================")
        println("   • Converged: $(res.converged)  ($(res.iters) iterations with |F| = $(res.resnorm))")
        (T < 673 || T > 2573 || P > 10.0) && println("     ⚠ (T,P) = ($T K, $P GPa) is outside the Zhang & Duan (2009) COH EOS calibration range (673-2573 K, ≤10 GPa) — fluid speciation is extrapolated here.")
        println()
        println("   • Regime: $(rg.regime)")
        @printf("   • Equilibrium log fO₂ = %.4f  (ΔFMQ = %+.3f, ΔIW = %+.3f)\n", eqfO2, eqfO2-fmq, eqfO2-iw)
        # XCO₂ is still evaluated with no melt, but Φ = 0 nulls it; do not print a phantom melt.
        println(rg.melt_absent ? "   • Melt XCO₂ = — (no melt)" : @sprintf("   • Melt XCO₂ = %.6g", XCO₂))
        carbonstate = if rg.carbonate && rg.graphite ; "graphite/diamond + carbonate coexist (on EMOG(D))"
                      elseif rg.carbonate            ; "solid carbonate present (aC=$(round(rg.aC, sigdigits=4)))"
                      elseif rg.graphite             ; "graphite/diamond present (aC=1)"
                      else                           ; "carbon-undersaturated, no solid C (aC=$(round(rg.aC, sigdigits=4)))"
                      end
        fluidflag = rg.has_fluid ? "$carbonstate   |   COH fluid is net $(rg.fluid_reducing ? "REDUCING (Oₑₓ source)" : "OXIDISING (Oₑₓ sink)")" : "$carbonstate   |   No COH fluid (no free H₂O)"
        if !res.converged
            println("\n   ⚠ CAPACITY EXCEEDED: even with every sink at its cap the assemblage cannot hold the whole")
            println("       TOₑₓ budget — physically this calls for a new Fe³⁺-rich phase to precipitate.")
            # Σ can miss on either side, and the two mean opposite things.
            if shortfall > 1e-6
                println("     ⇒ UNPLACED: $(round(1e2shortfall, digits=4))% of TOₑₓ has no sink to go into.")
            elseif shortfall < -1e-6
                println("     ⇒ OVER-SUBSCRIBED: the sinks demand $(round(-1e2shortfall, digits=4))% MORE oxygen than TOₑₓ holds.")
            end
        end
        println()
        println(" -----------------------------")
        println("   • $fluidflag")
        println(" -----------------------------")
        println("   • Mass balance residual → $(st.total - 1)")
        spin = rg.solid_sat ? " [PINNED at cap]" : rg.solid_metal ? " [PINNED at metal floor — all Fe as Fe⁰, not buffering]" : ""
        mpin = rg.melt_sat  ? " [PINNED at cap]" : rg.melt_floored  ? " [PINNED at floor — holds ~no Fe³⁺, not buffering]" : ""
        # An absent phase has no oxybarometer; printing its sentinel value would read as a real fO₂.
        sread = rg.solid_absent ? "ABSENT (ϕ=1)" : "$(round(st.fₛ,digits=4))$spin"
        mread = rg.melt_absent  ? "ABSENT (ϕ=0)" : "$(round(st.fₗ,digits=4))$mpin"
        println("   • Sanity fO₂ check → solid $sread, melt $mread")
        if rg.solid_absent || rg.melt_absent
            println("     one phase is not present — the remaining sinks carry the whole budget.")
        elseif rg.fO2_equal
            println("     ✓ solid and melt fO₂ agree — the equilibrium the model rests on.")
        elseif rg.pinned
            println("     fO₂ differ by $(round(abs(st.fₛ-st.fₗ),digits=4)) — legitimate: a pinned sink has left the")
            println("     equilibrium and no longer defines an fO₂, so only the unpinned phases share the value above.")
        else
            println("     ⚠ fO₂ differ by $(round(abs(st.fₛ-st.fₗ),digits=4)) with NEITHER sink pinned — this should be")
            println("     impossible (coexisting phases must share one fO₂); treat this result as suspect.")
        end
        println(" -----------------------------")
        # Fe³⁺/ΣFe > 1 is unphysical, so flag it. The 1.05 tolerance is not slop: shares are per
        # CURRENT solid mass but caps per REFERENCE mass, so a floor-pinned solid reads ~1.02.
        unphys(x) = x > 1.05 ? "  ⚠ >1, UNPHYSICAL" : ""
        sfe = rg.solid_absent ? "— (no solid)" : @sprintf("%.6f%s", FeR_solid, unphys(FeR_solid))
        smf = rg.solid_absent ? "— (no solid)" : @sprintf("%.6f  (%.4g wt%% of the solid)", FeM_solid, 1e2metal_wt)
        mfe = rg.melt_absent  ? "— (no melt)"  : @sprintf("%.6f%s", FeR_melt,  unphys(FeR_melt))
        println("   • Oxidation state → solid Fe³⁺/ΣFe = $sfe")
        println("                       solid Fe⁰/ΣFe  = $smf")
        println("                       melt  Fe³⁺/ΣFe = $mfe")
        # Which route made the metal. MAGEMin reports only a total; the split comes from the oxygen
        # budget, and is exact — see the note where FeM_diss/FeM_disp are computed.
        if FeM_solid > 1e-12
            @printf("                         ├ dissociation       %.6f   FeO → Fe⁰ + ½O₂, gives O back to TOₑₓ\n", FeM_diss)
            @printf("                         └ disproportionation %.6f   3Fe²⁺ → 2Fe³⁺ + Fe⁰, oxygen-neutral\n", FeM_disp)
        end
        # A net-oxygen release larger than the metal MAGEMin actually reports is not a split, it is a
        # contradiction between the two — flag it rather than absorb it into the clamp above.
        (FeM_net > FeM_solid + 1e-9) && @printf("                       ⚠ sOₑₓ implies Fe⁰/ΣFe = %.6f but MAGEMin's assemblage holds %.6f — table too coarse near the metal edge? lower Rmax\n", FeM_net, FeM_solid)
        println(" -----------------------------")
        @printf("   • Lever bounds: %.4f ≤ sOₑₓ ≤ %.4f (stoichiometric; the floor is all Fe as Fe⁰), mOₑₓ ≤ %.4f, XCO₂ ≤ %.4g\n",
                -minsOₑₓ, maxsOₑₓ, maxmOₑₓ, maxXCO₂)
        # Only MAGEMin-derived inside [tbl_lo, tbl_hi]; say so rather than let extrapolation pass as data.
        @printf("   • MAGEMin table spans sOₑₓ ∈ [%.4g, %.4g]; solution at %.4g", tbl_lo, tbl_hi, sOₑₓ)
        if rg.solid_absent
            @printf("  — no solid phase (ϕ=1); MAGEMin was not called\n")
        elseif sOₑₓ > tbl_hi
            @printf("  ⚠ ABOVE the table — solid fO₂ extrapolated (no MAGEMin phase changes captured)\n")
        elseif sOₑₓ < tbl_lo
            @printf("  ⚠ BELOW the table — past the most metal-saturated point sampled; solid fO₂ held flat\n")
        else
            @printf("  ✓ within tabulated range\n")
        end
        tbl_mono || println("     ⚠ MAGEMin fO₂ was non-monotonic in Fe³⁺ here (phase change); repaired by isotonic regression.")
        println(" -----------------------------")
        # Guarded: an absent phase has a zero cap, and TC=0 with ϕ=0 makes maxXCO₂ itself 0/0.
        pct(x, cap) = cap > 0 ? "$(round(1e2x/cap, digits=2))%" : "n/a"
        println("   • Distances to cap → sOₑₓ $(pct(sOₑₓ, maxsOₑₓ)), mOₑₓ $(pct(mOₑₓ, maxmOₑₓ)), XCO₂ $(pct(XCO₂, maxXCO₂))")
        if Mh > 0.0
            println(" -----------------------------")
            println("   • COH fluid (mole fractions) → XH₂O=$(round(XH2O,digits=4)), XCO₂=$(round(XCO2f,digits=4)), XCH₄=$(round(XCH4,digits=6)), XCO=$(round(XCOf,digits=6)), XH₂=$(round(XH2,digits=6))")
        end
        println()

        # Shares of TOₑₓ, not of a renormalised pool, with the fluid row signed: renormalising would
        # destroy the one number a host cannot re-derive, how much oxygen H₂O dissociation ADDED.
        released  = max(0.0, -cohsink)   # oxygen handed back by H₂O breakdown, in TOₑₓ units
        pool      = 1.0 + released       # what the rock sinks actually share out
        println("  ════════ PARTITIONED OXYGEN BUDGET ════════")
        println("  (shares of TOₑₓ; they sum to 100 % WITH the signed fluid row)\n")
        # The solid row is signed like the fluid row: below the FeO reference it hands oxygen back.
        @printf("     solid  %-9s       %+10.4f %%   %s\n",
                sOₑₓ < 0 ? "Fe⁰" : "Fe₂O₃", 1e2sOₑₓ,
                sOₑₓ < 0 ? "← SOURCE: FeO primarily dissociating to metallic Fe" :
                sOₑₓ > 0 ? "→ SINK: FeO primarily oxidising to Fe₂O₃"            : "(no net exchange)")
        @printf("     melt   Fe₂O₃            %9.4f %%\n", 1e2mOₑₓ)
        @printf("     melt   CO₂              %9.4f %%\n", 1e2st.carbonO)
        @printf("     solid  carbonate        %9.4f %%\n", 1e2st.carbO)
        @printf("     COH    fluid           %+10.4f %%   %s\n", 1e2cohsink,
                cohsink < 0 ? "← SOURCE: hands oxygen back to the rock" :
                cohsink > 0 ? "→ SINK: draws oxygen from the budget"    : "(no exchange)")
        println("     ---------------------------------")
        # Print the actual miss: at 4 decimals a small one still displays as exactly 100.0000%.
        @printf("     Σ                       %9.4f %%%s\n", 1e2st.total,
                res.converged ? "" : @sprintf("   ⚠ off by %.3g %% — no valid solution, see warning above", 1e2*abs(st.total-1)))
        println()
        if Mh > 0.0
            @printf("     fluid detail →  uptake  %9.4f %%  (O consumed making CO₂/CO)\n", 1e2cs.O_up)
            @printf("                     release %9.4f %%  (O released reducing H₂O to CH₄/H₂)\n", 1e2cs.O_rel)
            if released > 1e-12
                @printf("     ⇒ Oₑₓ ADDED to the tracer inventory by H₂O dissociation: %.4f %% of TOₑₓ\n", 1e2released)
            end
            @printf("     ⇒ the four ROCK sinks share out %.4f %% of TOₑₓ  (rock_pool = %.6f)\n",
                    1e2*(1.0-cohsink), 1.0-cohsink)
            println()
        end

        if Mc > 0.0
            # Melt and fluid carbon come from the joint solve; the solid remainder is a result, not an
            # assumption, because aC=1 is exactly the condition for free graphite/diamond to exist.
            pctC_melt  = 1e2*cs.nC_melt*mm.C/Mc
            pctC_fluid = 1e2*cs.nC_fluid*mm.C/Mc
            pctC_gph   = 1e2*cs.nC_gph*mm.C/Mc
            pctC_carb  = 1e2*cs.nC_carb*mm.C/Mc
            println("  ════════ PARTITIONED CARBON BUDGET ════════\n")
            @printf("     solid  graphite/diamond %9.4f %%\n", pctC_gph)
            @printf("     solid  carbonate        %9.4f %%\n", pctC_carb)
            @printf("     melt   CO₂              %9.4f %%\n", pctC_melt)
            @printf("     COH    fluid            %9.4f %%\n", pctC_fluid)
            println("     ---------------------------------")
            @printf("     Σ                       %9.4f %%\n", pctC_gph + pctC_carb + pctC_melt + pctC_fluid)
        else
            println("  ════════ PARTITIONED CARBON BUDGET ════════\n")
            println("     TC = 0 — no carbon in this tracer.")
        end
        println()
    end

    # Fractions of their own total, ready to scale tracer fields: the five oxygen shares sum to 1
    # exactly when converged, the fC_* carbon shares sum to 1 when TC > 0.
    nCt = Mc > 0 ? Mc/mm.C : 0.0
    fC(x) = nCt > 0 ? x/nCt : 0.0
    # cohsink is SIGNED (negative ⇒ the fluid donates); Oex_added is what H₂O dissociation ADDED and
    # must be booked against hydrogen. rock_pool = 1 - cohsink, NOT 1 + Oex_added, which over-reports.
    if out; return (converged=res.converged, fO2=eqfO2, FMQ=eqfO2-fmq, IW=eqfO2-iw, aC=st.cs.aC, sOₑₓ=sOₑₓ, mOₑₓ=mOₑₓ, XCO₂=XCO₂,
                    carbonO=st.carbonO, carbO=st.carbO, cohsink=cohsink,
                    fluid_O_uptake=cs.O_up, fluid_O_release=cs.O_rel,
                    Oex_added=max(0.0, -cohsink), rock_pool=1.0-cohsink,
                    FeR_solid=FeR_solid, FeR_net=FeR_net, FeM_solid=FeM_solid, FeM_net=FeM_net,
                    FeM_diss=FeM_diss, FeM_disp=FeM_disp, metal_wt=metal_wt,
                    nFe_solid=nFe_solid, FeR_melt=FeR_melt,
                    fC_gph=fC(st.cs.nC_gph), fC_carb=fC(st.cs.nC_carb),
                    fC_melt=fC(st.cs.nC_melt), fC_fluid=fC(st.cs.nC_fluid),
                    iters=res.iters, resnorm=res.resnorm,
                    XH2O=cs.XH2O, XCO2=cs.XCO2, XCH4=cs.XCH4, XCO=cs.XCO, XH2=cs.XH2); end

end
