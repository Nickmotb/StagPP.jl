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
#     solid carb.   carbonate, 2 O per carbon           (EMOG(D) from Holland & Powell 98 data)
#     C-O-H fluid   H₂O/CO₂/CH₄/CO/H₂ free fluid        (Holland & Powell 98 + Zhang & Duan 2009)
#
#   Every sink holds more oxygen as fO2 rises, so their demands are monotone in fO2. The fluid is
#   the exception worth understanding: it can be a net oxygen SOURCE. Water arrives already
#   carrying its own oxygen, so reducing H₂O to CH₄/H₂ hands that oxygen back to the rock, while
#   oxidising graphite/diamond to CO₂/CO consumes budget.
#
#   Carbon is governed by its activity aC, solved jointly with the fluid rather than assumed. Its
#   ceiling is set by whichever solid host is stable, which is what the EMOG(D) buffer decides:
#     f < EMOG(D)   graphite/diamond stable    ⇒ aC ≤ 1
#     f > EMOG(D)   carbonate stable           ⇒ aC ≤ 10^(f_EMOG - f), falling as fO2 rises
#   aC pinned at that ceiling means the solid host is present and holds the carbon left over after
#   melt and fluid; aC below it means no solid carbon at all. Since aC multiplies all three carbon
#   reactions and the melt buffer, it scales the whole carbon system at once, enforcing carbon mass
#   balance exactly and turning both the graphite-out and the graphite→carbonate boundaries into
#   smooth transitions rather than discontinuities.
#
#   Carbonate matters for oxygen, not just carbon: MgCO₃ carries 2 excess O per carbon against the
#   reduced reference, so converting graphite to carbonate is itself a large oxygen sink.
#
# METHOD
#   fO2 itself is the single unknown. Each sink states the share of the budget it demands at a
#   trial fO2 — obtained by inverting its own oxybarometer, saturating only at genuine physical
#   limits — and the one equation solved is total mass balance, Σ demands = 1.
#
#   Two nested bisections, both on monotone functions:
#     outer   fO2      solving oxygen mass balance  (Σ demands = 1)
#     inner   aC       solving carbon mass balance  (melt C + fluid C = TC, capped at aC=1)
#   The inner solve runs with the fugacity coefficients frozen, which are then refreshed from the
#   resulting composition, so the expensive equation of state is still evaluated once per outer
#   iteration rather than once per aC trial.
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
#   :capacity_exceeded       no fO2 reaches the budget even with every sink full — more oxidant than
#                            the assemblage can hold ⇒ needs a precipitating Fe³⁺ phase
#   :oxidant_exhausted       both Fe³⁺ sinks full; the rock's capacity binds the whole budget
#   :solid_saturated         solid at its Fe³⁺ ceiling, no longer tracking fO2
#   :melt_saturated          melt at its Fe³⁺ ceiling, no longer tracking fO2
#   :solid_depleted          solid holds ~no Fe³⁺, so it no longer buffers
#   :melt_depleted           melt holds ~no Fe³⁺, so it no longer buffers
#   :carbonate_graphite_buffered  both solid carbon hosts coexist — sitting on EMOG(D)
#   :carbonate_saturated     carbonate is the solid carbon host; graphite/diamond gone
#   :carbon_undersaturated   no solid carbon at all; every atom dissolved in melt+fluid
#   :buffered                nothing pinned, graphite/diamond present; the "normal" case
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
#   Fluid fugacities  Zhang & Duan (2009) C-O-H equation of state, coefficients taken from a working
#                     reference implementation rather than transcribed from the paper.
#
# BENCHMARKS (what has actually been checked against independent data)
#   ZD09 molar volumes    0.01% mean deviation vs that paper's own Table 6 (pure H₂O, 1200-1900 K).
#   K3, H₂+½O₂=H₂O        within 0.02 log units of tabulated ΔfG°(T) at 1000-1500 K.
#   graphite ⇄ diamond    within 0.14 GPa of Kennedy & Kennedy (1976) over 1000-2100 K — a direct
#                         test of the HP98 carbon data underlying every carbon reaction here.
#   EMOG(D)               ΔFMQ -2.3 to -3.4 over 2-6 GPa, inside the ~FMQ-1 to -3 range reported by
#                         Stagno & Frost (2010).
#   fluid speciation      at matched XO the composition tracks GFluid (Zhang & Duan 2010) Table 2
#                         closely (XH₂O 0.83 vs 0.848, XCH₄ 0.16 vs 0.139 at XO=0.2672).
#
# KNOWN LIMITATIONS

#  --- Major:
#   • ABSOLUTE fO2 OFFSET, ~1.9 log units. Matveev et al. (1997), IW-buffered at QFM-4.5, measured
#     XCH₄=0.814 / XH₂O=0.117 at 1273 K / 2.4 GPa. This model reproduces that composition, but only
#     at ΔFMQ = -6.35 — so the fO2 at which methane takes over is too reducing by ~1.9. Relative
#     behaviour (how species trade off as fO2 moves, and the composition at a given XO) benchmarks
#     well; it is the mapping onto an absolute fO2 scale that is shifted. Partly attributable: K4
#     (C+2H₂=CH₄) drifts ~0.3-0.5 log units from tabulated ΔfG° above 1200 K, and `logfO2_FMQ`
#     carries ~0.72 log units/GPa of pressure dependence against a published ~0.5-0.6, both pushing
#     the same way. Treat ΔFMQ values as indicative to a couple of log units until this is resolved.
#   • No precipitating Fe³⁺-rich phase, so a budget larger than the assemblage can hold is reported
#     as `capacity_exceeded` rather than resolved. Note that in testing every such case was pinned by
#     the sampled Fe³⁺/ΣFe range (`Rmax`), not by stoichiometry — raise `Rmax` before concluding a
#     new phase is needed, since MAGEMin already precipitates whatever Fe³⁺ host is stable.
#   • Melt CO₂ is capped at `XCO2max` (default 0.3) as a stand-in for CO₂ saturation in a silicate
#     melt. Without it the EMOG(D) inversion drives XCO₂→1 near the carbonate boundary, where the
#     melt is no longer a silicate melt and nC_melt = X·nOx/(1-X) diverges, letting a vanishingly
#     small melt swallow the entire carbon budget. This is the least constrained number in the file
#     and the first thing to revisit with a real CO₂ solubility model.
#   • Carbonate is treated as pure MgCO₃ (magnesite) at unit activity. Real carbonate is a
#     magnesite-dolomite-calcite solution, which shifts EMOG(D) and smears it over a finite fO2
#     interval — approximated here by the `wEMOG` ramp rather than modelled.

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

# Tabulates solid fO2 against the share of TOₑₓ held as Fe³⁺, extrapolating physically beyond both
# ends of the table. Returns (interpolant, table lo, table hi) so callers can flag extrapolated use.
# `s_stoich` is the loading at which every Fe²⁺ is oxidised — the true ceiling the upper branch
# approaches, independent of how far MAGEMin was actually sampled.
function solid_fO2_interpolant(P, T, X, Xox, SymXox, Ms, TOₑₓ, s_stoich; nr=25, Rmax=0.2, data=nothing)
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

    # `Flat()` is qualitatively wrong at BOTH ends, so replace it with the oxybarometer's own form,
    # log-linear in Fe³⁺/Fe²⁺ with the dilute-limit exponent 1/a the melt branch uses:
    #   below the table  fO2 → -∞ as Fe³⁺ → 0     (a flat floor would stop the solid expressing
    #                                              reduced states and break fO2 equality)
    #   above the table  fO2 → +∞ as Fe²⁺ → 0     (a flat roof caps the sink at however far MAGEMin
    #                                              happened to be sampled, which is what made
    #                                              `capacity_exceeded` a tabulation artefact)
    # Above the table the ratio is written R/(1-R) with R the oxidised fraction, so it diverges at
    # the stoichiometric limit rather than at the edge of the sampling.
    g1, f1 = grid[1], fO2[1]
    gN, fN = grid[end], fO2[end]
    slope  = 1/0.1917
    RN     = gN/s_stoich
    upper  = RN < 1 - 1e-9                       # is there any room left above the table at all?
    lnN    = upper ? log10(RN/(1 - RN)) : 0.0

    function sfO2(s)
        s < g1  && return f1 + slope*(log10(max(s, 1e-300)) - log10(g1))
        s <= gN && return itp(s)
        upper   || return fN
        R = min(s/s_stoich, 1 - 1e-12)
        return fN + slope*(log10(R/(1 - R)) - lnN)
    end

    return sfO2, grid[1], grid[end]
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

# One speciation step at fixed φᵢ and carbon activity aC: a quadratic in XH2.
# aC multiplies every carbon-bearing reaction (C+O₂=CO₂, C+½O₂=CO, C+2H₂=CH₄); aC=1 is graphite/
# diamond saturation, aC<1 means no free carbon phase remains.
# Returns (valid, XH2O, XCO2, XCH4, XCO, XH2); `valid=false` ⇒ the carbon species alone would exceed
# the fluid, which `_solve_aC` avoids by lowering aC.
function _coh_reaction_step(fO₂, Pbar, K1, K2, K3, K4, φCH4, φH2O, φCO2, φH2, φCO, aC)
    XCO2 = aC*K1*fO₂/(φCO2*Pbar)
    XCO  = aC*K2*√fO₂/(φCO*Pbar)
    aq2  = aC*K4*φH2^2*Pbar/φCH4
    aq1  = K3*φH2*√fO₂/φH2O + 1
    c_q  = 1 - XCO2 - XCO
    disc = aq1^2 + 4aq2*c_q
    if disc < 0 || c_q <= 0
        return false, 0.5, 0.5, 0.0, 0.0, 0.0
    end
    # Stable root of aq2·XH2² + aq1·XH2 - c_q = 0: the textbook form divides by aq2, which cancels
    # catastrophically as aC→0 (aq2→0) and returns 0 instead of the correct linear limit c_q/aq1.
    XH2 = 2c_q/(aq1 + √disc)
    (XH2 < 0) && return false, 0.5, 0.5, 0.0, 0.0, 0.0
    XH2O = K3*φH2*√fO₂/φH2O * XH2
    XCH4 = aq2*XH2^2
    return true, XH2O, XCO2, XCH4, XCO, XH2
end

# Moles of carbon that melt + fluid together demand at a given carbon activity.
# Monotonically increasing in aC, and +Inf where the carbon species would overfill the fluid, which
# is what lets `_solve_aC` bracket it. `Bmelt` is the melt XCO₂ the EMOG(D) buffer wants at aC=1.
function _carbon_demand(aC, fO₂, Pbar, K1, K2, K3, K4, φ, molH, nOxideTotal, Bmelt, XCO2max)
    Xm = softcap(aC*Bmelt, XCO2max)
    nC_melt = Xm*nOxideTotal/(1 - Xm)

    # With no free water there is no fluid phase, so its composition — and its validity — are
    # irrelevant. Gating on them would let a non-existent phase veto the carbon balance.
    molH <= 0 && return nC_melt

    ok, XH2O, XCO2, XCH4, XCO, XH2 = _coh_reaction_step(fO₂, Pbar, K1, K2, K3, K4, φ[1], φ[2], φ[3], φ[4], φ[5], aC)
    ok || return Inf
    nH = 2XH2O + 4XCH4 + 2XH2
    nC_fluid = nH > 0 ? molH*(XCO2 + XCH4 + XCO)/nH : 0.0
    return nC_melt + nC_fluid
end

# Carbon activity satisfying carbon mass balance, capped by whichever solid carbon host is stable
# (`aC_cap`): aC pinned at the cap ⇒ that solid phase is present and holds the remainder; aC below
# it ⇒ no solid carbon at all, every atom dissolved. Bisected in log₁₀ aC as the demand spans decades.
function _solve_aC(fO₂, Pbar, K1, K2, K3, K4, φ, molH, nC_total, nOxideTotal, Bmelt, aC_cap, XCO2max; iters=120)
    nC_total <= 0 && return 0.0                        # no carbon at all ⇒ no carbon species
    d1 = _carbon_demand(aC_cap, fO₂, Pbar, K1, K2, K3, K4, φ, molH, nOxideTotal, Bmelt, XCO2max)
    (isfinite(d1) && d1 <= nC_total) && return aC_cap   # a solid carbon phase is present
    # Bracket must reach far below any plausible aC: the melt's buffer coefficient can exceed 1e40
    # at high fO2, so a floor of only 1e-40 leaves carbon still over-committed at the bracket edge.
    lo, hi = -300.0, log10(aC_cap)                     # log₁₀ aC bracket
    for _ in 1:iters
        mid = 0.5*(lo + hi)
        _carbon_demand(10.0^mid, fO₂, Pbar, K1, K2, K3, K4, φ, molH, nOxideTotal, Bmelt, XCO2max) < nC_total ?
            (lo = mid) : (hi = mid)
    end
    return 10.0^(0.5*(lo + hi))
end

# Solves the fluid, its carbon activity and its oxygen draw together at one fO2.
# The φ loop and the aC solve are nested: φ is frozen while aC is bisected (pure arithmetic), then
# refreshed from the resulting composition, so the expensive EOS is still evaluated only once per
# outer iteration.
function coh_state(f, T, P, Mh, Mc, nOxideTotal, TOₑₓ; outer_iters=6, wEMOG=0.25, XCO2max=0.3)
    Pbar  = 1e4*P                  # GPa → bar, matching the 1-bar standard state of the data
    fO₂   = 10.0^f
    K1, K2, K3, K4 = reaction_Ks(T, Pbar)
    molH     = 2*Mh/mm.H2O
    nC_total = Mc/mm.C
    Bmelt    = invert_carbon_buffer_fO2(f, P, T)       # melt XCO₂ at aC=1

    # Carbon-activity ceiling set by whichever solid host is stable. Below EMOG(D) that is
    # graphite/diamond (aC=1); above it, carbonate, whose own equilibrium fixes aC = 10^(f_EMOG - f)
    # and so drives the carbon activity down as conditions oxidise.
    fE     = logfO2_EMOG(T, Pbar)
    aC_cap = min(1.0, 10.0^(fE - f))

    φ  = ones(7)
    aC = aC_cap
    ok, XH2O, XCO2, XCH4, XCO, XH2 = true, 0.5, 0.5, 0.0, 0.0, 0.0
    for _ in 1:outer_iters
        aC = _solve_aC(fO₂, Pbar, K1, K2, K3, K4, φ, molH, nC_total, nOxideTotal, Bmelt, aC_cap, XCO2max)
        ok, XH2O, XCO2, XCH4, XCO, XH2 = _coh_reaction_step(fO₂, Pbar, K1, K2, K3, K4, φ[1], φ[2], φ[3], φ[4], φ[5], aC)
        vx = zeros(7); vx[1]=XCH4; vx[2]=XH2O; vx[3]=XCO2; vx[4]=XH2; vx[5]=XCO
        s = sum(vx); (s > 0) && (vx ./= s)
        _, fug = zd_fugacity(vx, Pbar, T)
        for i in 1:7
            newφ = fug[i]/(vx[i]*Pbar)
            (isfinite(newφ) && newφ > 0) && (φ[i] = newφ)
        end
    end

    nH       = 2XH2O + 4XCH4 + 2XH2
    n_fluid  = (Mh > 0 && nH > 0) ? molH/nH : 0.0
    nC_fluid = n_fluid*(XCO2 + XCH4 + XCO)

    XCO2_melt = softcap(aC*Bmelt, XCO2max)
    nC_melt   = XCO2_melt*nOxideTotal/(1 - XCO2_melt)

    # Carbon left over once melt and fluid have taken theirs sits in a solid phase, but only while
    # aC is pinned at its ceiling — below the ceiling no solid carbon coexists at all.
    solid_present = nC_total > 0 && aC >= aC_cap*(1 - 1e-9)
    nC_solid = solid_present ? max(0.0, nC_total - nC_melt - nC_fluid) : 0.0

    # Which solid host: graphite/diamond below EMOG(D), carbonate above. Ramped over `wEMOG` log
    # units rather than switched, because the univariant boundary is a step in oxygen demand (2 O per
    # carbon converted) that would otherwise put a jump in the residual. Real carbonate is a
    # magnesite-dolomite-calcite solution, so the transition is genuinely smeared, not just softened.
    x_carb   = 1/(1 + exp(-(f - fE)/wEMOG))
    nC_carb  = x_carb*nC_solid                 # carbonate: MgCO₃ vs MgO+C ⇒ 2 excess O per carbon
    nC_gph   = nC_solid - nC_carb              # graphite/diamond: carries no excess oxygen
    carbO    = 2*nC_carb*mm.O/TOₑₓ

    # Net excess-oxygen draw against the reduced reference (H as H₂O, so H₂O itself is neutral):
    # CO₂/CO consume O (+2/+1), CH₄/H₂ release it (-2/-1). Negative is physical — a reduced fluid
    # donates oxygen back to the rock.
    netO = n_fluid*(2XCO2 + XCO - 2XCH4 - XH2)*mm.O / TOₑₓ
    fₑₓ  = clamp(netO, -1e3, 1e3)   # |·|>1 already exceeds the whole budget; bound for the bracket

    # Thresholded against the total: the logistic ramp has exponential tails, so a bare `> 0` would
    # report a phase present on a numerically negligible fraction far from the boundary.
    ctol = 1e-6*max(nC_total, eps())
    return (; aC, aC_cap, fE, graphite = nC_gph > ctol, carbonate = nC_carb > ctol, valid = ok,
              XH2O, XCO2, XCH4, XCO, XH2, n_fluid, nC_fluid,
              XCO2_melt, nC_melt, nC_solid, nC_carb, nC_gph, carbO, fₑₓ)
end

# --- EMOG(D): the graphite/diamond ⇄ carbonate boundary in a peridotitic assemblage. ---
# 2 MgCO₃ + Mg₂Si₂O₆ = 2 Mg₂SiO₄ + 2C + 2O₂, so log₁₀fO2 = -ΔG/(2·R·T·ln10) at unit activities.
# Below it free carbon is graphite/diamond; above it the stable solid host is carbonate.
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

# Reference buffers for reporting, from Frost (1991) Rev. Mineral. 25, Table 1 (P in bar). These are
# the standard empirical parameterisations the literature quotes ΔFMQ against, so using them keeps
# reported offsets comparable with published values. Reporting only — neither enters the solve.
#
# Computing FMQ instead from the hp02ver mineral data (as `logfO2_FMQ_hp98` below still does) puts it
# 1.3-1.8 log units too oxidised, drifting with T, which mislabels an otherwise correct fO2: the
# fluid speciation reproduces the IW-buffered Matveev et al. (1997) experiment at fO2 = -13.98
# against IW = -13.85 here, i.e. to 0.12 log units on an absolute scale.
logfO2_FMQ(T, Pbar) = -25096.3/T + 8.735 + 0.110*(Pbar - 1)/T
logfO2_IW(T, Pbar)  = -27489.0/T + 6.702 + 0.055*(Pbar - 1)/T

# hp02ver-derived FMQ, kept for comparison against the internally consistent dataset.
# (3Fe₂SiO₄ + O₂ = 2Fe₃O₄ + 3SiO₂)
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
    (; sfO2, T, P, IDV, TOₑₓ, Mf, Mc, molfracXB, nB_per_mass, Φ, Mh, nOxideTotal, maxmOₑₓ, maxsOₑₓ, lowclip) = ctx

    sOₑₓ = _invert_demand(s -> sfO2(s), f, lowclip, maxsOₑₓ - lowclip)
    mOₑₓ = _invert_demand(m -> melt_fO2(m, T, TOₑₓ, Mf, molfracXB, nB_per_mass, IDV), f,
                          lowclip, maxmOₑₓ - lowclip)

    # Fluid, carbon activity and melt CO₂ come out of one joint solve: aC enforces carbon mass
    # balance exactly, so melt XCO₂ needs no ceiling of its own.
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

# Classifies the converged state (see REGIMES in the header). Never an input to the solve, so the
# label cannot disagree with what was actually solved.
function classify_regime(st, ctx, res)
    atcap(x, cap)  = cap > 0 && x >= 0.999*cap
    atfloor(x)     = x <= 1.000001*ctx.lowclip
    solid_sat  = atcap(st.sOₑₓ, ctx.maxsOₑₓ);  solid_floored = atfloor(st.sOₑₓ)
    melt_sat   = atcap(st.mOₑₓ, ctx.maxmOₑₓ);  melt_floored  = atfloor(st.mOₑₓ)
    graphite   = st.cs.graphite       # free graphite/diamond coexists
    carbonate  = st.cs.carbonate      # solid carbonate coexists (above EMOG(D))

    regime = if !res.converged        ; :capacity_exceeded
             elseif solid_sat && melt_sat ; :oxidant_exhausted
             elseif solid_sat             ; :solid_saturated
             elseif melt_sat              ; :melt_saturated
             elseif solid_floored         ; :solid_depleted
             elseif melt_floored          ; :melt_depleted
             elseif carbonate && graphite ; :carbonate_graphite_buffered
             elseif carbonate             ; :carbonate_saturated
             elseif !graphite             ; :carbon_undersaturated
             else                         ; :buffered
             end
    # fO2 equality may legitimately fail only for a sink pinned at one of its two ends.
    pinned = solid_sat || solid_floored || melt_sat || melt_floored
    return (; regime, solid_sat, melt_sat, solid_floored, melt_floored, pinned, graphite, carbonate,
              aC = st.cs.aC, fE = st.cs.fE, fO2_equal = abs(st.fₛ - st.fₗ) <= 1e-6,
              has_fluid = ctx.Mh > 0, fluid_reducing = st.fₑₓ < 0)
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
                           nr=25, Rmax=0.4, maxiter=500, tol=1e-8, verbose=false, data=nothing, out=false, debug=false)

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
    maxXCO₂ = clamp(maxXCO₂_raw, 0.0, 1.0)   # all carbon dissolved in the melt; reported, not imposed
    lowclip = 1e-10

    # The ceiling is the stoichiometric one; the table's reach no longer caps it, since `sfO2`
    # extrapolates above the sampled range along the oxybarometer's own form. `tbl_lo`/`tbl_hi`
    # are kept only so the report can say when the answer relies on that extrapolation.
    sfO2, tbl_lo, tbl_hi = solid_fO2_interpolant(P, T, X, Xox, SymXox, Ms, TOₑₓ, maxsOₑₓ;
                                                 nr=nr, Rmax=Rmax, data=data)

    ctx = (; sfO2, T, P, IDV, TOₑₓ, Mf, Mc, molfracXB, nB_per_mass, Φ, Mh, nOxideTotal, maxmOₑₓ, maxsOₑₓ, lowclip)

    # Bracket far wider than any physical value: every demand saturates outside its own range, so
    # this only has to straddle the root, it is not a constraint on the answer.
    res  = bisect_solve(f -> unified_residual(f, ctx), -120.0, 40.0; maxiter=maxiter, atol=tol, debug=debug)
    st   = _coh_derive(res.x, ctx)
    sOₑₓ, mOₑₓ, XCO₂ = st.sOₑₓ, st.mOₑₓ, st.XCO₂

    eqfO2     = st.f       # the shared fO2 all unpinned sinks equilibrate to
    cohsink   = st.fₑₓ
    shortfall = 1 - st.total

    if verbose
        cs = st.cs
        XH2O, XCO2f, XCH4, XCOf, XH2 = cs.XH2O, cs.XCO2, cs.XCH4, cs.XCO, cs.XH2
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
        spin = rg.solid_sat ? " [PINNED at cap]" : rg.solid_floored ? " [PINNED at floor — holds ~no Fe³⁺, not buffering]" : ""
        mpin = rg.melt_sat  ? " [PINNED at cap]" : rg.melt_floored  ? " [PINNED at floor — holds ~no Fe³⁺, not buffering]" : ""
        println("   • Sanity fO₂ check → solid $(round(st.fₛ,digits=4))$spin, melt $(round(st.fₗ,digits=4))$mpin")
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
        println("   • Lever upper bounds: sOₑₓ ≤ $maxsOₑₓ (stoichiometric), mOₑₓ ≤ $maxmOₑₓ, XCO₂ ≤ $maxXCO₂")
        # The solid oxybarometer is only MAGEMin-derived inside [tbl_lo, tbl_hi]; outside it the
        # value comes from the log-linear extrapolation, so say so rather than let it pass as data.
        @printf("   • MAGEMin table spans sOₑₓ ∈ [%.4g, %.4g]; solution at %.4g", tbl_lo, tbl_hi, sOₑₓ)
        if sOₑₓ > tbl_hi
            @printf("  ⚠ ABOVE the table — solid fO₂ extrapolated (no MAGEMin phase changes captured)\n")
        elseif sOₑₓ < tbl_lo
            @printf("  ⚠ BELOW the table — solid fO₂ extrapolated (holds ~no Fe³⁺)\n")
        else
            @printf("  ✓ within tabulated range\n")
        end
        println(" -----------------------------")
        println("   • Distances to cap → sOₑₓ $(round(1e2sOₑₓ/maxsOₑₓ, digits=2))%, mOₑₓ $(round(1e2mOₑₓ/maxmOₑₓ, digits=2))%, XCO₂ $(round(1e2XCO₂/maxXCO₂, digits=2))%")
        if Mh > 0.0
            println(" -----------------------------")
            println("   • COH fluid (mole fractions) → XH₂O=$(round(XH2O,digits=4)), XCO₂=$(round(XCO2f,digits=4)), XCH₄=$(round(XCH4,digits=6)), XCO=$(round(XCOf,digits=6)), XH₂=$(round(XH2,digits=6))")
        end
        println()

        # Shares are quoted against the pool actually available: a reducing fluid breaks H₂O down and
        # adds that oxygen to TOₑₓ, so it appears in the pool header rather than as a negative row.
        # This keeps every entry positive and the column summing to 100%.
        released  = max(0.0, -cohsink)   # oxygen handed back by H₂O breakdown, in TOₑₓ units
        pool      = 1.0 + released
        println("  ════════ PARTITIONED OXYGEN BUDGET ════════\n")
        if released > 1e-12
            @printf("     available pool = TOₑₓ + %.4f %% of TOₑₓ released by H₂O→CH₄/H₂ (shares are of that pool)\n",
                    1e2released)
            println("     ---------------------------------")
        end
        @printf("     solid  Fe₂O₃            %9.4f %%\n", 1e2sOₑₓ/pool)
        @printf("     melt   Fe₂O₃            %9.4f %%\n", 1e2mOₑₓ/pool)
        @printf("     melt   CO₂              %9.4f %%\n", 1e2st.carbonO/pool)
        @printf("     solid  carbonate        %9.4f %%\n", 1e2st.carbO/pool)
        if cohsink > 0
            @printf("     COH    fluid            %9.4f %%\n", 1e2cohsink/pool)
        end
        println("     ---------------------------------")
        # Print the actual miss: at 4 decimals a small one still displays as exactly 100.0000%.
        @printf("     Σ                       %9.4f %%%s\n", 1e2st.total/pool + 1e2released/pool,
                res.converged ? "" : @sprintf("   ⚠ off by %.3g %% — no valid solution, see warning above", 1e2*abs(st.total-1)/pool))
        println()

        if Mc > 0.0
            # Melt and fluid carbon come straight from the joint solve; the solid remainder is now a
            # genuine result rather than an assumption, because aC=1 is exactly the condition under
            # which free graphite/diamond can exist.
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

    # Both budgets are returned as fractions of their own total, ready to scale tracer fields:
    #   sOₑₓ/mOₑₓ/carbonO/carbO/cohsink  — fractions of TOₑₓ, summing to 1 when converged
    #   fC_*                             — fractions of TC,   summing to 1 when TC > 0
    # (cohsink may be negative: a reducing fluid hands oxygen back rather than consuming it.)
    nCt = Mc > 0 ? Mc/mm.C : 0.0
    fC(x) = nCt > 0 ? x/nCt : 0.0
    if out; return (converged=res.converged, fO2=eqfO2, aC=st.cs.aC, sOₑₓ=sOₑₓ, mOₑₓ=mOₑₓ, XCO₂=XCO₂,
                    carbonO=st.carbonO, carbO=st.carbO, cohsink=cohsink,
                    fC_gph=fC(st.cs.nC_gph), fC_carb=fC(st.cs.nC_carb),
                    fC_melt=fC(st.cs.nC_melt), fC_fluid=fC(st.cs.nC_fluid),
                    iters=res.iters, resnorm=res.resnorm); end

end
