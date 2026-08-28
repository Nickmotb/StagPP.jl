# =============================================================================================
# misc_corrected.jl — partitioning of excess oxygen (Oₑₓ) between coexisting solid, melt and
#                     C-O-H fluid at fixed pressure and temperature
# =============================================================================================
#
# PURPOSE
#   Given a tracer's bulk excess-oxygen, carbon and free-water budgets at some (P,T), decide how
#   that oxygen distributes over the four reservoirs that can host it, and report the resulting
#   oxygen and carbon budgets in the form a geodynamic model needs to move mass between tracer
#   fields. Entry point: `partition_Oₑₓ`.
#
# THE PHYSICAL PICTURE
#   Four sinks compete for one shared oxygen budget, all of them tied together by a single
#   intensive variable — the oxygen fugacity fO2 that coexisting phases must share at equilibrium:
#
#     solid Fe₂O₃   Fe³⁺ in the solid residue           (MAGEMin, Stixrude & Lithgow-Bertelloni 24)
#     melt  Fe₂O₃   Fe³⁺ dissolved in the melt          (Hirschmann 2022 oxybarometer)
#     melt  CO₂     carbon dissolved in the melt        (EMOG(D) buffer, Stagno & Frost 2010-style)
#     C-O-H fluid   H₂O/CO₂/CH₄/CO/H₂ free fluid        (Holland & Powell 98 + Zhang & Duan 2009)
#
#   Every sink holds more oxygen as fO2 rises, so their demands are monotone in fO2. The fluid is
#   the exception worth understanding: it can be a net oxygen SOURCE. Water arrives already
#   carrying its own oxygen, so reducing H₂O to CH₄/H₂ hands that oxygen back to the rock, while
#   oxidising graphite/diamond to CO₂/CO consumes budget.
#
# METHOD
#   fO2 itself is the single unknown. Each sink states the share of the budget it demands at a
#   trial fO2 — obtained by inverting its own oxybarometer, saturating only at genuine physical
#   limits — and the one equation solved is total mass balance, Σ demands = 1.
#
#   Two consequences worth stating explicitly:
#     • fO2 equality between phases is exact BY CONSTRUCTION, not an equation the solver hunts
#       for: an unsaturated sink is *defined* as the loading whose oxybarometer reads the shared
#       fO2. Equality may fail only for a sink pinned at one of its two physical ends, which is
#       correct — a phase holding no Fe³⁺, or all of it, no longer buffers.
#     • Because every demand is monotone in fO2, the total is monotone, so the root exists and is
#       unique whenever the assemblage has the capacity at all. The solve is therefore plain
#       bisection: no Jacobian, no initial guess, no possibility of divergence.
#
#   The solve is regime-free; `classify_regime` diagnoses the converged state for reporting only,
#   so a label can never drift out of sync with what was actually solved.
#
# REGIMES (diagnostic labels, in the order they are tested)
#   :capacity_exceeded          no fO2 reaches the budget even with every sink full — more oxidant
#                               than the assemblage can hold ⇒ needs a precipitating Fe³⁺ phase
#   :carbon_undersaturated_gap  residual jumps across zero at the graphite/diamond ceiling; an
#                               artefact of the missing aC<1 fluid, not physics
#   :oxidant_exhausted          both Fe³⁺ sinks full; the rock's capacity binds the whole budget
#   :solid_saturated            solid at its Fe³⁺ ceiling, no longer tracking fO2
#   :melt_saturated             melt at its Fe³⁺ ceiling, no longer tracking fO2
#   :solid_depleted             solid holds ~no Fe³⁺, so it no longer buffers
#   :melt_depleted              melt holds ~no Fe³⁺, so it no longer buffers
#   :carbon_saturated           melt XCO₂ pinned at the carbon budget
#   :buffered                   nothing pinned; every sink sits on its own oxybarometer
#
# CONVENTIONS
#   P in GPa (converted to bar internally)
#   T in Kelvin.
#   TOex, TC, H are raw mass fractions of the tracer (0-1, not percentages).
#   sOₑₓ, mOₑₓ, carbonO, and the fluid term are all fractions of the TOₑₓ budget.
#   All fO2 quantities are log₁₀ values.
#
#   Excess oxygen is measured against a REDUCED reference state — Fe as FeO, C as graphite/diamond,
#   H as H₂O — used consistently by every sink. Hence 0.5 O per Fe³⁺, 2 O per dissolved CO₂ carbon,
#   and 0 for water: the oxygen arriving bonded inside H₂O is not part of TOₑₓ, and counting it
#   would double-count it.
#
# DATA SOURCES
#   Solid fO2         MAGEMin `sb24` database, tabulated over solid Fe³⁺/ΣFe and interpolated.
#   Thermodynamics    Holland & Powell (1998), via Perple_X's `hp02ver.dat`.
#   Fluid fugacities  Zhang & Duan (2009) C-O-H equation of state; the coefficients were taken
#                     from a working reference implementation and validated to 0.01% against that
#                     paper's own Table 6 molar volumes. Speciation reproduces the Matveev et al.
#                     (1997) 1273 K / 2.4 GPa benchmark.
#
# KNOWN LIMITATIONS

#  --- Major:
#   • Carbon-undersaturated (aC<1) fluid is NOT implemented. Above the graphite/diamond-saturation
#     ceiling no valid composition exists and the fluid falls back to a flat placeholder; results
#     there are flagged and must not be read as physical. This is the main open gap.
#   • Solid carbonate is not a sink, so carbon that cannot be graphite/diamond has no host.
#   • No precipitating Fe³⁺-rich phase, so a budget larger than the assemblage can hold is reported
#     as `capacity_exceeded` rather than resolved.

#  --- Minor:
#   • C₂H₆ is deliberately excluded from the fluid: its mass-action term scales as P², which
#     amplifies any inconsistency between an external ethane dataset and the data used here.
#   • `logfO2_FMQ` is this tool's own internally consistent FMQ reference, not a universal one —
#     published FMQ curves differ by ~1 log unit depending on the mineral dataset.
#   • Zhang & Duan (2009) is calibrated for 673-2573 K and ≤10 GPa; outside that the fluid
#     speciation is extrapolated and flagged.

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
    nB_per_mass = sum(XB./mmXox)                  # moles of B oxides per unit mass (NOT normalized)
    molMf       = Mf*nB_per_mass                  # moles of the melt's own oxide assemblage

    molCav      = Mc/mm.C                    # moles of available carbon
    maxXCO₂_raw = molCav/(molMf + molCav)    # carbon-limited ceiling on melt XCO₂

    Φ = 2*mm.O*molMf/TOₑₓ   # converts a melt XCO₂ into a fraction of the TOₑₓ budget

    return (; Mf, Ms, Mc, Mh, TOₑₓ, molfracXB, nB_per_mass, molMf, maxXCO₂_raw, Φ)
end

# ------------------------------------------------------------
# Section 2 — solid branch: MAGEMin (Stixrude & Lithgow-Bertelloni 2024)
# ------------------------------------------------------------

# Tabulates solid fO2 against the share of TOₑₓ held as Fe³⁺, returning (interpolant, table reach).
function solid_fO2_interpolant(P, T, X, Xox, SymXox, Ms, TOₑₓ; nr=25, Rmax=0.2, data=nothing)
    Rlist  = LinRange(1e-5, Rmax, nr)   # scan of solid Fe³⁺/ΣFe ratios
    Xdummy = (X=Vector{Float64}(X), Xox=Vector{String}(Xox), mm=get_Xoxmm(Xox))
    Xlist  = Vector{Vector{Float64}}(undef, nr)
    sOwt   = zeros(nr)                  # wt fraction of O in the solid at each Rlist[i]
    for i in 1:nr
        Xl       = oxidize_bulk(X, Rlist[i], Xdummy, wt_out=true, frac=true, FeFormat="FeO_O", SymXox=SymXox)
        Xlist[i] = [getfield(Xl, f) for f in SymXox]
        sOwt[i]  = Xl.O
    end

    flag = isnothing(data)
    flag && (data = Initialize_MAGEMin("sb24", verbose=false))
    out  = multi_point_minimization(10P*ones(nr), k2c(T)*ones(nr), data, X=Xlist, Xoxides=Vector{String}(Xox), name_solvus=true, sys_in="wt", progressbar=false)
    flag && Finalize_MAGEMin(data)

    grid = sOwt .* (Ms/TOₑₓ)   # re-express the x-axis as a fraction of the TOₑₓ budget
    fO2  = [out[i].fO2 for i in eachindex(out)]

    # MAGEMin can fail to converge at some Rlist points; keep only what it actually solved.
    good = isfinite.(fO2)
    if !all(good)
        @warn "solid_fO2_interpolant: MAGEMin returned $(count(!, good))/$nr non-finite fO2 value(s) (Rlist up to Rmax=$(Rlist[end])) — dropping those grid points."
        grid, fO2 = grid[good], fO2[good]
    end

    itp = extrapolate(Interpolations.interpolate((grid,), fO2, Gridded(Linear())), Flat())

    # Below the table, extrapolate log-linearly rather than flat: fO2 → -∞ as Fe³⁺ → 0, and a flat
    # floor would stop the solid expressing reduced states and break fO2 equality. Slope is the
    # dilute-limit exponent 1/a of the same Fe³⁺/Fe²⁺ oxybarometer the melt branch uses.
    g1, f1, slope = grid[1], fO2[1], 1/0.1917
    sfO2(s) = s >= g1 ? itp(s) : f1 + slope*(log10(max(s, 1e-300)) - log10(g1))

    return sfO2, grid[end]
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

# Gibbs energy at temperature T, anchored on the tabulated G0(298.15 K) with Cp supplying the shift.
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
# μC uses whichever carbon polymorph is actually stable at (T,Pbar) — graphite converts to diamond
# above ~2-5 GPa, and using graphite past that exaggerates its pressure-driven destabilization.
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

# Molar volume [cc/mol] of a fluid of composition vx.
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

# One graphite-saturated speciation step at fixed φᵢ: a quadratic in XH2.
# Returns (valid, XH2O, XCO2, XCH4, XCO, XH2). `valid=false` ⇒ no such composition exists here, and
# the returned one is a flat neutral placeholder (flat: no spurious roots; neutral: no tiny divisor).
function _coh_reaction_step(fO₂, Pbar, K1, K2, K3, K4, φCH4, φH2O, φCO2, φH2, φCO)
    XCO2 = K1*fO₂/(φCO2*Pbar)
    XCO  = K2*√fO₂/(φCO*Pbar)
    aq2  = K4*φH2^2*Pbar/φCH4
    aq1  = K3*φH2*√fO₂/φH2O + 1
    c_q  = 1 - XCO2 - XCO
    disc = aq1^2 + 4aq2*c_q
    if disc < 0 || c_q <= 0
        return false, oftype(fO₂,0.5), oftype(fO₂,0.5), oftype(fO₂,0.0), oftype(fO₂,0.0), oftype(fO₂,0.0)
    end
    XH2 = (-aq1 + √disc)/(2aq2)
    (XH2 < 0) && return false, oftype(fO₂,0.5), oftype(fO₂,0.5), oftype(fO₂,0.0), oftype(fO₂,0.0), oftype(fO₂,0.0)
    XH2O = K3*φH2*√fO₂/φH2O * XH2
    XCH4 = aq2*XH2^2
    return true, XH2O, XCO2, XCH4, XCO, XH2
end

# Converges the composition-dependent φᵢ against the reaction system, then returns the speciation.
# The φ loop runs on `ForwardDiff.value(fO₂)` since the EOS root-find is a bisection, so ∂φᵢ/∂fO2 (a
# second-order effect) is dropped; ∂Xᵢ/∂fO2 stays exact via the final step's original fO₂.
function _coh_solve(fO₂, T, P; outer_iters=6)
    Pbar = 1e4*P   # GPa → bar, to match the 1-bar standard state of the thermodynamic data
    K1, K2, K3, K4 = reaction_Ks(T, Pbar)
    fO2_val = ForwardDiff.value(fO₂)
    φ = ones(7)
    for _ in 1:outer_iters
        _, XH2Ov, XCO2v, XCH4v, XCOv, XH2v = _coh_reaction_step(fO2_val, Pbar, K1, K2, K3, K4, φ[1], φ[2], φ[3], φ[4], φ[5])
        vx = zeros(7); vx[1]=XCH4v; vx[2]=XH2Ov; vx[3]=XCO2v; vx[4]=XH2v; vx[5]=XCOv
        s = sum(vx); (s > 0) && (vx ./= s)
        _, fug = zd_fugacity(vx, Pbar, T)
        for i in 1:7
            newφ = fug[i]/(vx[i]*Pbar)
            (isfinite(newφ) && newφ > 0) && (φ[i] = newφ)
        end
    end
    return _coh_reaction_step(fO₂, Pbar, K1, K2, K3, K4, φ[1], φ[2], φ[3], φ[4], φ[5])
end

# fO2 of the fayalite-magnetite-quartz buffer (3Fe₂SiO₄ + O₂ = 2Fe₃O₄ + 3SiO₂), for reporting ΔFMQ.
const FMQ_FA = (G0=-1378945.0, S=151.0, c1=201.1, c2=0.01733,   c3=-1960600.0, c5=-900.9,
                 V0=4.631,  b1=0.505e-4, b5=-0.505e-3, b6=1330000.0, b7=-199.5, b8=4.0)
const FMQ_QZ = (G0=-856424.6,  S=41.5,  c1=110.7, c2=-0.005189, c3=0.0,        c5=-1128.3,
                 V0=2.2688, b1=0.65e-5,  b5=-0.65e-4,  b6=750000.0,  b7=-112.5, b8=4.0)
const FMQ_MT = (G0=-1012288.0, S=146.1, c1=262.5, c2=-0.007204, c3=-1926200.0, c5=-1655.7,
                 V0=4.452,  b1=0.696e-4, b5=-0.696e-3, b6=1850000.0, b7=-277.5, b8=4.0)
function logfO2_FMQ(T, Pbar)
    dG0 = 2hp98_gibbs(T,FMQ_MT) + 3hp98_gibbs(T,FMQ_QZ) - 3hp98_gibbs(T,FMQ_FA) - hp98_gibbs(T,HP98.O2)
    dGP = 2murnaghan_dG(T,Pbar,FMQ_MT) + 3murnaghan_dG(T,Pbar,FMQ_QZ) - 3murnaghan_dG(T,Pbar,FMQ_FA)
    return (dG0 + dGP)/(R*T*log(10))
end

# C-O-H fluid mole fractions at a given fO2, T, P.
function coh_speciation(fO₂, T, P)
    _, XH2O, XCO2, XCH4, XCO, XH2 = _coh_solve(fO₂, T, P)
    return XH2O, XCO2, XCH4, XCO, XH2
end

# Net excess-oxygen draw of the fluid, as a fraction of TOₑₓ. Speciation fixes the fluid's
# composition; hydrogen conservation against the water budget Mh fixes how much fluid exists.
function coh_fOₑᵪ(logfO₂, T, P, Mh, TOₑₓ)
    _, XH2O, XCO2, XCH4, XCO, XH2 = _coh_solve(10^logfO₂, T, P)
    molH_avail       = 2*Mh/mm.H2O
    nH_per_mol_fluid = 2XH2O + 4XCH4 + 2XH2
    n_fluid          = molH_avail/nH_per_mol_fluid

    # Measured against the reduced reference (H as H₂O), so H₂O itself is neutral: CO₂/CO consume
    # excess O (+2/+1), CH₄/H₂ release it (-2/-1). A negative result is physical — a reduced fluid
    # donates oxygen — and the term stays monotone in fO2 either way.
    netO_per_mol_fluid = 2XCO2 + XCO - 2XCH4 - XH2
    raw = n_fluid*netO_per_mol_fluid*mm.O / TOₑₓ

    # |value| > 1 already means this term alone exceeds the whole budget; bound it so a small TOₑₓ
    # with a large Mh cannot produce a magnitude that swamps the bracket search.
    return clamp(raw, oftype(raw, -1e3), oftype(raw, 1e3))
end

# ------------------------------------------------------------
# Section 5 — smooth saturation and finiteness helpers
# ------------------------------------------------------------

@inline function softmax2(a, b, β)
    m = max(a, b)
    return m + log(exp(β*(a-m)) + exp(β*(b-m)))/β
end
@inline softmin2(a, b, β) = -softmax2(-a, -b, β)

# Smooth ceiling: identity well below `cap`, rolling over without a kink as it is approached.
# Caps from above only — a lower rounding would inject a phantom floor into XCO₂, which Φ then
# amplifies into a spurious oxygen demand. `max(0,…)` trims the few-ulp undershoot at x≈0.
@inline softcap(x, cap; sharpness=30.0) = max(zero(x), softmin2(x, cap, sharpness/cap))

# Throws on Inf as well as NaN: an Inf compares silently and would otherwise poison the solve.
@inline function nancheck(label, val, inputs...)
    isfinite(val) || error("$label evaluated to $val (non-finite) — inputs: $inputs")
    return val
end

# ------------------------------------------------------------
# Section 6 — the residual (1 unknown: the shared log₁₀ fO2)
# ------------------------------------------------------------

# Inverts a monotonically increasing oxybarometer for the loading that reads the target fO2,
# saturating at the bounds when the target lies outside the achievable range.
function _invert_demand(g, target, lo, hi; iters=80)
    g(lo) >= target && return lo      # floor already too oxidized for this fO2 → holds nothing
    g(hi) <= target && return hi      # cap still too reduced → this sink is saturated
    for _ in 1:iters
        mid = 0.5*(lo + hi)
        g(mid) < target ? (lo = mid) : (hi = mid)
    end
    return 0.5*(lo + hi)
end

# Maps a trial fO2 to the full state: each sink's demand, its own fO2 reading, and the total.
# Used both while solving and when reporting, so the reported state is the one the residual saw.
function _coh_derive(f, ctx)
    (; sfO2, T, P, IDV, TOₑₓ, Mf, molfracXB, nB_per_mass, Φ, Mh, maxXCO₂, maxmOₑₓ, maxsOₑₓ, lowclip) = ctx

    sOₑₓ = _invert_demand(s -> sfO2(s), f, lowclip, maxsOₑₓ - lowclip)
    mOₑₓ = _invert_demand(m -> melt_fO2(m, T, TOₑₓ, Mf, molfracXB, nB_per_mass, IDV), f,
                          lowclip, maxmOₑₓ - lowclip)

    # Carbon inverts in closed form; softcapped at the carbon budget so entering carbon saturation
    # stays smooth, this being the sink most likely to sit exactly at its cap.
    XCO₂    = softcap(invert_carbon_buffer_fO2(f, P, T), maxXCO₂ - lowclip)
    carbonO = Φ*XCO₂/(1 - XCO₂)

    fₑₓ = nancheck("coh_fOₑᵪ", coh_fOₑᵪ(f, T, P, Mh, TOₑₓ), f, Mh)
    fₛ  = nancheck("sfO2(sOₑₓ)", sfO2(sOₑₓ), sOₑₓ)
    fₗ  = nancheck("melt_fO2", melt_fO2(mOₑₓ, T, TOₑₓ, Mf, molfracXB, nB_per_mass, IDV), mOₑₓ)

    return (; f, sOₑₓ, mOₑₓ, XCO₂, carbonO, fₑₓ, fₛ, fₗ, total = sOₑₓ + mOₑₓ + carbonO + fₑₓ)
end

# The one equation: total demand must equal the budget.
unified_residual(f, ctx) = _coh_derive(f, ctx).total - 1

# ------------------------------------------------------------
# Section 6b — regime diagnosis (reporting only)
# ------------------------------------------------------------

# Classifies the converged state (see REGIMES in the header). Never an input to the solve, so the
# label cannot disagree with what was actually solved.
function classify_regime(st, ctx, res)
    atcap(x, cap)  = cap > 0 && x >= 0.999*cap
    atfloor(x)     = x <= 1.000001*ctx.lowclip
    solid_sat  = atcap(st.sOₑₓ, ctx.maxsOₑₓ);  solid_floored = atfloor(st.sOₑₓ)
    melt_sat   = atcap(st.mOₑₓ, ctx.maxmOₑₓ);  melt_floored  = atfloor(st.mOₑₓ)
    carbon_sat = atcap(st.XCO₂, ctx.maxXCO₂)
    graphite_stable = first(_coh_solve(10.0^st.f, ctx.T, ctx.P))

    regime = if !res.converged
                 res.iters == 0 ? :capacity_exceeded : :carbon_undersaturated_gap
             elseif solid_sat && melt_sat     ; :oxidant_exhausted
             elseif solid_sat                 ; :solid_saturated
             elseif melt_sat                  ; :melt_saturated
             elseif solid_floored             ; :solid_depleted
             elseif melt_floored              ; :melt_depleted
             elseif carbon_sat                ; :carbon_saturated
             else                             ; :buffered
             end
    # fO2 equality may legitimately fail only for a sink pinned at one of its two ends.
    pinned = solid_sat || solid_floored || melt_sat || melt_floored
    return (; regime, solid_sat, melt_sat, carbon_sat, solid_floored, melt_floored, pinned,
              fO2_equal = abs(st.fₛ - st.fₗ) <= 1e-6,
              graphite_stable, has_fluid = ctx.Mh > 0, fluid_reducing = st.fₑₓ < 0)
end

# ------------------------------------------------------------
# Section 7 — bracketed bisection solver
# ------------------------------------------------------------

# Bisection suits a 1-D monotone bracketed residual: no derivatives, no initial guess, no divergence.
# Convergence is judged on the RESIDUAL alone — a collapsed bracket with |F| still large means the
# residual is discontinuous there and no root exists, which must not be reported as success.
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

# Partitions the excess-oxygen budget at (P,T); `out` returns the solution as a NamedTuple.
# `p`: 0 → depleted (harzburgitic) solid, 1 → enriched (basaltic). `tol`: accepted mass-balance
# residual — 1e-8 clears bisection's ~1e-11 floor yet stays far below the ~0.1-3.5 graphite cliff.
function partition_Oₑₓ(P::Real, T::Real, p::Real, ϕ::Real, TOex::Real, TC::Real, H::Real;
                           nr=25, Rmax=0.18, maxiter=500, tol=1e-8, verbose=false, data=nothing, out=false, debug=false)

    (; X, XB, Xox, SymXox, mmXox) = endmember_composition(p)
    Mt = 1e3   # reference cell mass (1 kg in grams, matching mm.* in g/mol)
    (; Mf, Ms, Mc, Mh, TOₑₓ, molfracXB, nB_per_mass, maxXCO₂_raw, Φ) = mass_budget(Mt, ϕ, TC, TOex, H, XB, mmXox)

    if TOₑₓ == 0.0
        verbose && println("Total Oₑₓ budget = 0.0 — nothing to partition.")
        return (converged=true, sOₑₓ=0.0, mOₑₓ=0.0, XCO₂=0.0, iters=0, resnorm=0.0)
    end

    IDV = nancheck("solve_∫ΔVdP", solve_∫ΔVdP([P-0.05P, P, P+0.05P], [T-0.05T, T, T+0.05T])[2,2,1], P, T)

    # Compositional ceilings: every Fe²⁺ in that phase oxidised. NOT clipped at 100% of the budget —
    # that is the mass-balance constraint, already enforced by the residual, and imposing it twice
    # would pin a sink off its oxybarometer for a bookkeeping reason and break fO2 equality.
    molX    = (X./mmXox) ./ sum(X./mmXox)
    maxsOₑₓ = (0.5molX[3]*mm.O)/(sum(X) - 0.5molX[3]*mm.O)*(Ms/TOₑₓ)

    nOxideTotal = Mf*nB_per_mass
    maxmOₑₓ = 0.5*molfracXB[3]*nOxideTotal*mm.O/TOₑₓ
    maxXCO₂ = clamp(maxXCO₂_raw, 0.0, 1.0)   # a mole fraction, so 1.0 here IS a real bound
    lowclip = 1e-10

    sfO2, sOₑₓ_table_max = solid_fO2_interpolant(P, T, X, Xox, SymXox, Ms, TOₑₓ; nr=nr, Rmax=Rmax, data=data)
    maxsOₑₓ = min(maxsOₑₓ, sOₑₓ_table_max)   # never search past what MAGEMin actually tabulated

    ctx = (; sfO2, T, P, IDV, TOₑₓ, Mf, molfracXB, nB_per_mass, Φ, Mh, maxXCO₂, maxmOₑₓ, maxsOₑₓ, lowclip)

    # Bracket far wider than any physical value: every demand saturates outside its own range, so
    # this only has to straddle the root, it is not a constraint on the answer.
    res  = bisect_solve(f -> unified_residual(f, ctx), -120.0, 40.0; maxiter=maxiter, atol=tol, debug=debug)
    st   = _coh_derive(res.x, ctx)
    sOₑₓ, mOₑₓ, XCO₂ = st.sOₑₓ, st.mOₑₓ, st.XCO₂

    eqfO2     = st.f       # the shared fO2 all unpinned sinks equilibrate to
    cohsink   = st.fₑₓ
    shortfall = 1 - st.total

    if verbose
        XH2O, XCO2f, XCH4, XCOf, XH2 = coh_speciation(10^eqfO2, T, P)
        rg = classify_regime(st, ctx, res)
        println("========================================================================================================")
        println("==== Solution (P=$(P) GPa | T=$(T) K | ϕ=$(1e2ϕ) % | TOex=$(1e2TOex) wt% | mix=$p | TC=$(1e2TC) wt% | H=$(1e2H) wt%) ====")
        println("========================================================================================================")
        println("   • Converged: $(res.converged)  ($(res.iters) iterations with |F| = $(res.resnorm))")
        (T < 673 || T > 2573 || P > 10.0) && println("     ⚠ (T,P) = ($T K, $P GPa) is outside the Zhang & Duan (2009) COH EOS calibration range (673-2573 K, ≤10 GPa) — fluid speciation is extrapolated here.")
        println()
        println("   • Regime: $(rg.regime)")
        fmq = logfO2_FMQ(T, 1e4P)
        println("   • Equilibrium log fO₂ = $eqfO2  (ΔFMQ = $(round(eqfO2-fmq, digits=3)))")
        println("   • Melt XCO₂ = $XCO₂")
        fluidflag = rg.has_fluid ? "Graphite/Diamond stable: $(rg.graphite_stable)   |   COH fluid is net $(rg.fluid_reducing ? "REDUCING (Oₑₓ source)" : "OXIDISING (Oₑₓ sink)")" : "No COH fluid (no free H₂O)"
        if !res.converged
            if res.iters == 0
                println("\n   ⚠ CAPACITY EXCEEDED: even with every sink at its cap the assemblage cannot hold the whole")
                println("       TOₑₓ budget — physically this calls for a new Fe³⁺-rich phase to precipitate.")
            else
                println("\n   ⚠ NO ROOT: the mass-balance residual jumps across zero rather than crossing it, so no fO₂")
                println("       balances the budget. This is the graphite/diamond-saturation boundary — past it the COH")
                println("       model falls back to a placeholder composition, and that discontinuity is an artefact of")
                println("       the still-missing carbon-undersaturated (aC<1) fluid regime, not real physics.")
            end
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
        spin = rg.solid_sat ? " [PINNED at cap]" : rg.solid_floored ? " [PINNED at floor — holds ~no Fe³⁺, not buffering]" : ""
        mpin = rg.melt_sat  ? " [PINNED at cap]" : rg.melt_floored  ? " [PINNED at floor — holds ~no Fe³⁺, not buffering]" : ""
        println("   • Sink fO₂ check → solid reads $(round(st.fₛ,digits=4))$spin, melt reads $(round(st.fₗ,digits=4))$mpin")
        if rg.fO2_equal
            println("     ✓ solid and melt fO₂ agree — the equilibrium the model rests on.")
        elseif rg.pinned
            println("     fO₂ differ by $(round(abs(st.fₛ-st.fₗ),digits=4)) — legitimate: a pinned sink has left the")
            println("     equilibrium and no longer defines an fO₂, so only the unpinned phases share the value above.")
        else
            println("     ⚠ fO₂ differ by $(round(abs(st.fₛ-st.fₗ),digits=4)) with NEITHER sink pinned — this should be")
            println("     impossible (coexisting phases must share one fO₂); treat this result as suspect.")
        end
        println(" -----------------------------")
        println("   • Lever upper bounds: sOₑₓ ≤ $maxsOₑₓ (MAGEMin table reached up to $sOₑₓ_table_max), mOₑₓ ≤ $maxmOₑₓ, XCO₂ ≤ $maxXCO₂")
        println(" -----------------------------")
        println("   • Distances to cap → sOₑₓ $(round(1e2sOₑₓ/maxsOₑₓ, digits=2))%, mOₑₓ $(round(1e2mOₑₓ/maxmOₑₓ, digits=2))%, XCO₂ $(round(1e2XCO₂/maxXCO₂, digits=2))%")
        if Mh > 0.0
            println(" -----------------------------")
            lbl = rg.graphite_stable ? "" : "   ⚠ PLACEHOLDER (fO₂ above graphite/diamond saturation — not physical)"
            println("   • COH fluid (mole fractions) → XH₂O=$(round(XH2O,digits=4)), XCO₂=$(round(XCO2f,digits=4)), XCH₄=$(round(XCH4,digits=6)), XCO=$(round(XCOf,digits=6)), XH₂=$(round(XH2,digits=6))$lbl")
        end
        println()

        # Shares are quoted against the pool actually available: a reducing fluid breaks H₂O down and
        # adds that oxygen to TOₑₓ, so it appears in the pool header rather than as a negative row.
        # This keeps every entry positive and the column summing to 100%.
        released  = max(0.0, -cohsink)   # oxygen handed back by H₂O breakdown, in TOₑₓ units
        pool      = 1.0 + released
        poolmass  = TOₑₓ*pool
        println("  ════════ PARTITIONED OXYGEN BUDGET ════════\n")
        if released > 1e-12
            @printf("     available pool = %.6g g  =  TOₑₓ %.6g g  +  %.6g g released by H₂O→CH₄/H₂\n",
                    poolmass, TOₑₓ, TOₑₓ*released)
            println("     ---------------------------------")
        end
        @printf("     solid  Fe₂O₃            %9.4f %%   %12.6g g\n", 1e2sOₑₓ/pool, TOₑₓ*sOₑₓ)
        @printf("     melt   Fe₂O₃            %9.4f %%   %12.6g g\n", 1e2mOₑₓ/pool, TOₑₓ*mOₑₓ)
        @printf("     melt   CO₂              %9.4f %%   %12.6g g\n", 1e2st.carbonO/pool, TOₑₓ*st.carbonO)
        if cohsink > 0
            # Placeholder-derived when graphite is unstable, exactly like the carbon table below.
            @printf("     COH    fluid            %9.4f %%   %12.6g g%s\n", 1e2cohsink/pool, TOₑₓ*cohsink,
                    (rg.has_fluid && !rg.graphite_stable) ? "   ⚠ from PLACEHOLDER fluid — not physical" : "")
        end
        println("     ---------------------------------")
        # Print the actual miss: at 4 decimals a small one still displays as exactly 100.0000%.
        @printf("     Σ                       %9.4f %%%s\n", 1e2st.total/pool + 1e2released/pool,
                res.converged ? "" : @sprintf("   ⚠ off by %.3g %% — no valid solution, see warning above", 1e2*abs(st.total-1)/pool))
        println()

        if Mc > 0.0
            molC_melt  = XCO₂*nOxideTotal/(1-XCO₂)                     # moles of C dissolved as melt CO₂
            molH_avail = 2Mh/mm.H2O
            nH_per_mol_fluid = 2XH2O + 4XCH4 + 2XH2
            n_fluid    = nH_per_mol_fluid > 0 ? molH_avail/nH_per_mol_fluid : 0.0
            molC_fluid = n_fluid*(XCO2f + XCH4 + XCOf)                 # moles of C in the COH fluid
            pctC_melt  = 1e2*molC_melt*mm.C/Mc
            pctC_fluid = 1e2*molC_fluid*mm.C/Mc
            pctC_solid = 100.0 - pctC_melt - pctC_fluid

            # The solid row is a remainder, so it silently assumes graphite/diamond can host it. That
            # fails when no graphite-saturated fluid exists, or when XCO₂ is pinned at the carbon cap
            # (every atom already dissolved) — then the remainder is carbon with nowhere to go.
            gd_absent = rg.carbon_sat || (rg.has_fluid && !rg.graphite_stable)
            println("  ════════ PARTITIONED CARBON BUDGET ════════\n")
            if gd_absent && pctC_solid > 1e-6
                @printf("     UNHOSTED                %9.4f %%   ⚠ graphite/diamond is NOT stable at this fO₂\n", pctC_solid)
            else
                @printf("     solid  graphite/diamond %9.4f %%%s\n", pctC_solid,
                        pctC_solid < -1e-6 ? "   ⚠ NEGATIVE: melt+fluid demand more C than the tracer holds" : "")
            end
            @printf("     melt   CO₂              %9.4f %%\n", pctC_melt)
            @printf("     COH    fluid            %9.4f %%\n", pctC_fluid)
            println("     ---------------------------------")
            @printf("     Σ                       %9.4f %%\n", pctC_solid + pctC_melt + pctC_fluid)
            if gd_absent && pctC_solid > 1e-6
                println("     ⚠ That $(round(pctC_solid,digits=2))% has no valid host in the current model: graphite/diamond cannot exist")
                println("       at this fO₂, and solid carbonate is not implemented yet. It is NOT sitting in the solid as g/d.")
            elseif pctC_solid < -1e-6
                println("     ⚠ No graphite/diamond would remain — the carbon-undersaturated (aC<1) state, not yet modelled.")
            end
        else
            println("  ════════ PARTITIONED CARBON BUDGET ════════\n")
            println("     TC = 0 — no carbon in this tracer.")
        end
        println()
    end

    # sOₑₓ/mOₑₓ/carbonO/cohsink are fractions of TOₑₓ and sum to 1 when converged.
    if out; return (converged=res.converged, fO2=eqfO2, sOₑₓ=sOₑₓ, mOₑₓ=mOₑₓ, XCO₂=XCO₂, carbonO=st.carbonO, cohsink=cohsink, iters=res.iters, resnorm=res.resnorm); end

end
