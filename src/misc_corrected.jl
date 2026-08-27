# ============================================================
# misc_corrected.jl — ground-up rewrite of the Oₑₓ partitioning solver
# ============================================================
#
# Stage 1: the solid Fe2O3 / melt Fe2O3 / melt XCO2 fO2-equilibrium core
# (sOₑₓ, mOₑₓ, XCO₂; buffered_residual/carbon_saturated_residual).
# Stage 2: a C-O-H fluid sink layered into that core's mass-balance
# equation (Section 4b) — real thermodynamics (Huizenga 2001's method,
# Holland & Powell 1998 data for Kᵢ(T)), not a placeholder. It needs no
# extra unknown of its own: given the shared fO2 already produced by
# solving Stage 1, the fluid's speciation is a direct, closed-form
# function of (fO2, T, P). Still deferred: real (non-ideal) fugacity
# coefficients — currently ideal mixing (γᵢ=1), a labeled simplification
# matching what Huizenga notes is adequate above ~400°C, not an oversight.
#
# Why this differs from misc.jl, in four points:
#
#   1. Automatic differentiation (ForwardDiff) supplies every derivative
#      except the one branch MAGEMin sits behind (solid fO2 vs. sOₑₓ).
#      misc.jl hand-derived its Jacobian and several entries turned out to
#      reuse the wrong ∂x∂yₙ chain-rule factor, or the wrong sign, or
#      referenced undefined variables outright — AD removes that whole bug
#      class by construction. The catch: to stay ForwardDiff-compatible, no
#      function below may route a differentiated quantity through a
#      `Cbulk` (its fields are hard-typed ::Float64, so it can't hold a
#      Dual number) or through a mutated Float64 buffer (misc.jl's `dummy`
#      array) — both break Dual propagation silently. So the melt/carbon
#      physics is written as plain, non-mutating arithmetic on the unknowns
#      instead of reusing the oxidize_bulk/Cbulk machinery.
#
#   2. Physical bound constraints (mass fractions confined to (0, cap)) are
#      enforced by solving in an unconstrained logit-transformed space, not
#      by ad hoc per-step clipping against the bounds.
#
#   3. The nonlinear solve is Levenberg–Marquardt with an adaptively updated
#      damping parameter λ, not plain Newton. A synthetic 3-equation
#      stand-in test showed why misc.jl likely struggles at "high
#      parameters": as any of sOₑₓ/mOₑₓ/XCO₂ approaches its physical cap,
#      the logit transform saturates, its derivative → 0, and the Jacobian
#      in y-space becomes ill-conditioned — plain (even trust-region
#      capped) Newton stalls indefinitely right there, while LM's λI
#      regularization keeps the step well posed and converges from the same
#      starting point. Pushing TOex/TC/ϕ higher is exactly what pushes one
#      of the three fractions toward its cap, so this is very likely a
#      large part of the reported non-convergence.
#
#   4. The COH fluid (Section 4b) is solved directly, not iteratively.
#      misc.jl introduced fH₂ as a free unknown paired with a pressure-
#      balance equation (`P = COH_P(fH₂,fO2)`) — but that equation also fed
#      it the raw untransformed sigmoid value instead of the real fH₂, its
#      Jacobian entries for that row referenced undefined variables, and no
#      literature method for "P and T both externally fixed" (our actual
#      situation) needs fH₂ as an unknown at all. Huizenga (2001) §2.2 (the
#      standard "externally buffered fO2" method) solves the mole fractions
#      directly from fO2, T, and P once they're known — a closed-form
#      quadratic, no iteration, no extra unknown.

# ------------------------------------------------------------
# Section 1 — bulk composition & mass bookkeeping
# (pure Float64 preprocessing, done once before the solve — no AD needed)
# ------------------------------------------------------------

# Two mantle end-members ("H"arzburgitic and "B"asaltic) mixed by fraction p of B.
function endmember_composition(p::Real)
    XH     = @SVector [0.4347, 0.4597, 0.0835, 0.0090, 0.0100, 0.0001, 0.0030, 0.0] # mass fraction
    XB     = @SVector [0.5097, 0.0988, 0.0718, 0.1268, 0.1698, 0.0225, 0.0007, 0.0] # mass fraction
    Xox    = @SVector ["SiO2", "MgO", "FeO", "CaO", "Al2O3", "Na2O", "Cr2O3", "O"]
    SymXox = @SVector [:SiO2, :MgO, :FeO, :CaO, :Al2O3, :Na2O, :Cr2O3, :O]
    mmXox  = get_Xoxmm(Xox)
    X      = p*XB + (1-p)*XH   # bulk (solid-residue) composition, fed to the MAGEMin table
    return (; X, XB, Xox, SymXox, mmXox)
end

# All mass/mole bookkeeping needed before the solve starts.
# • `TOex` is the total excess-oxygen budget as a wt fraction of the solid tracer
# • `TOₑₓ` is the resulting absolute mass budget (grams — same convention as `Mt`).
function mass_budget(Mt, ϕ, TC, TOex, H, XB, mmXox)
    Mf   = ϕ*Mt          # melt tracer mass
    Ms   = Mt*(1-ϕ)      # solid tracer mass
    Mc   = TC*Mt         # carbon tracer mass
    Mh   = H*Mt          # free-water tracer mass, feeds the COH fluid sink
    TOₑₓ = TOex*Mt       # total excess-oxygen mass budget

    # Moles of the (unoxidized) basaltic end-member, per unit mass — kept UN-normalized on purpose.
    molfracXB   = (XB./mmXox) ./ sum(XB./mmXox)   # mole fractions of the basaltic end-member (sum == 1)
    nB_per_mass = sum(XB./mmXox)                  # moles of B end-member oxides per unit mass (NOT normalized)
    molMf       = Mf*nB_per_mass                  # total moles of the melt tracer's own oxide assemblage

    molCav      = Mc/mm.C                    # moles of available carbon
    maxXCO₂_raw = molCav/(molMf + molCav)    # carbon-limited ceiling on melt XCO₂

    Φ = 2*mm.O*molMf/TOₑₓ   # converts a melt XCO₂ mole fraction into a mass-fraction-of-TOₑₓ (mass-balance eq. 3)

    return (; Mf, Ms, Mc, Mh, TOₑₓ, molfracXB, nB_per_mass, molMf, maxXCO₂_raw, Φ)
end

# ------------------------------------------------------------
# Section 2 — solid branch: MAGEMin (Stixrude & Lithgow-Bertelloni 2024)
# ------------------------------------------------------------
function solid_fO2_interpolant(P, T, X, Xox, SymXox, Ms, TOₑₓ; nr=25, Rmax=0.2, data=nothing)
    Rlist  = LinRange(1e-5, Rmax, nr)   # scan of solid Fe³⁺/ΣFe ratios
    Xdummy = (X=Vector{Float64}(X), Xox=Vector{String}(Xox), mm=get_Xoxmm(Xox))
    Xlist  = Vector{Vector{Float64}}(undef, nr)
    sOwt   = zeros(nr)   # wt fraction of O in the solid, at each Rlist[i]
    for i in 1:nr
        Xl       = oxidize_bulk(X, Rlist[i], Xdummy, wt_out=true, frac=true, FeFormat="FeO_O", SymXox=SymXox)
        Xlist[i] = [getfield(Xl, f) for f in SymXox]
        sOwt[i]  = Xl.O
    end

    flag = isnothing(data)
    flag && (data = Initialize_MAGEMin("sb24", verbose=false))
    out  = multi_point_minimization(10P*ones(nr), k2c(T)*ones(nr), data, X=Xlist, Xoxides=Vector{String}(Xox), name_solvus=true, sys_in="wt", progressbar=false)
    flag && Finalize_MAGEMin(data)

    # Convert the x-axis from "wt fraction O in the solid" to "mass fraction
    # of the TOTAL Oₑₓ budget held in the solid" — the same units as the
    # sOₑₓ unknown the residual below actually solves for.
    grid = sOwt .* (Ms/TOₑₓ)
    fO2  = [out[i].fO2 for i in eachindex(out)]

    # MAGEMin's own Gibbs-energy minimization can fail to converge at some
    # Rlist points (more likely now that Rmax reaches further).
    # Drop any non-finite rows before building the interpolant so the table only ever contains points MAGEMin actually solved.
    good = isfinite.(fO2)
    if !all(good)
        @warn "solid_fO2_interpolant: MAGEMin returned $(count(!, good))/$nr non-finite fO2 value(s) (Rlist up to Rmax=$(Rlist[end])) — dropping those grid points."
        grid, fO2 = grid[good], fO2[good]
    end

    itp = extrapolate(Interpolations.interpolate((grid,), fO2, Gridded(Linear())), Flat())
    # grid[end] is the real edge of what MAGEMin actually explored AND solved successfully.
    return itp, grid[end]
end

# ------------------------------------------------------------
# Section 3 — melt branch: Hirschmann (2022) oxybarometer
# ------------------------------------------------------------
function melt_fO2(mOₑₓ, T, TOₑₓ, Mf, molfracXB, nB_per_mass, IDV)
    # Hirschmann (2022) regression parameters.
    a=0.1917; b=-1.961; c=4158.1; ΔCₚ=33.25; T₀=1673.15
    y1=-520.46; y3=494.39; y4=1838.34; y5=2888.48; y8=-1245.09; y9=-1156.86

    nOxideTotal = Mf*nB_per_mass           # absolute moles of the melt's own (unoxidized) oxide assemblage
    nOₑₓ_added  = TOₑₓ*mOₑₓ/mm.O           # absolute moles of excess O added to the melt
    θ           = nOxideTotal + nOₑₓ_added # renormalization after adding the excess oxygen

    XO     = nOₑₓ_added/θ
    XSiO2  = molfracXB[1]*nOxideTotal/θ
    XMgO   = molfracXB[2]*nOxideTotal/θ
    XFeO   = molfracXB[3]*nOxideTotal/θ
    XCaO   = molfracXB[4]*nOxideTotal/θ
    XAl2O3 = molfracXB[5]*nOxideTotal/θ
    XNa2O  = molfracXB[6]*nOxideTotal/θ

    (XO >= 0.5XFeO) && return oftype(mOₑₓ, 14.0)   # hard limit: more O than there is Fe to host it as Fe2O3

    _T, _ln10 = 1/T, 1/log(10)
    return (log10(XO/(XFeO - 2XO)) - b - c*_T + (ΔCₚ*_ln10/R)*(1 - T₀*_T - log(T/T₀)) + (IDV/(1e-3R))*_T*_ln10
              - _T*(y1*XSiO2 + y3*XMgO + y4*XCaO + y5*XNa2O + y8*XSiO2*XAl2O3 + y9*XSiO2*XMgO)) / a
end

# ------------------------------------------------------------
# Section 4 — melt CO2 branch: EMOG(D) buffer (Stagno & Frost, 2010-style)
# ------------------------------------------------------------
function carbon_buffer_fO2(XCO₂, P, T)
    @assert XCO₂ > 0.0 "XCO₂ must be strictly positive."
    return 5.44 - 21380/T + 0.078*(1e4P - 1)/T + log10(XCO₂)
end

# ------------------------------------------------------------
# Section 4b — C-O-H fluid speciation: Huizenga (2001) method,
# real Holland & Powell (1998) thermodynamics for Kᵢ(T)
# ------------------------------------------------------------
# Standard-state data (ΔfG°, ΔfH° in J/mol at 298.15 K, S° in J/mol/K from HP, 1998)
#       Cp(T) = c1 + c2·T + c3·T⁻² + c5·T⁻⁰·⁵ in [J/mol/K] for the six C-O-H fluid species plus graphite. 
# `hp98_gibbs` below anchors on the tabulated G0 directly and uses Cp only for the T-dependent shift, so
# that small inconsistency never enters the result.
const HP98 = (
    H2O = (G0=-228542.3,   H=-241810.0, S=188.8,  c1=40.1,  c2=8.656e-3,  c3=487500.0,  c5=-251.2),
    CO2 = (G0=-394332.9,   H=-393510.0, S=213.7,  c1=87.8,  c2=-2.644e-3, c3=706400.0,  c5=-998.9),
    CO  = (G0=-137163.7,   H=-110530.0, S=197.67, c1=45.7,  c2=-9.7e-5,   c3=662700.0,  c5=-414.7),
    CH4 = (G0=-50695.63,   H=-74810.0,  S=186.26, c1=150.1, c2=2.062e-3,  c3=3427700.0, c5=-2650.4),
    H2  = (G0=-0.0045494,  H=0.0,       S=130.7,  c1=23.3,  c2=4.627e-3,  c3=0.0,       c5=76.3),
    O2  = (G0=0.0,         H=0.0,       S=205.2,  c1=48.3,  c2=-6.91e-4,  c3=499200.0,  c5=-420.7),
    gph = (G0=-32.79654,   H=0.0,       S=5.85,   c1=51.0,  c2=-4.428e-3, c3=488600.0,  c5=-805.5),
)

# Gibbs energy of a species at temperature T (K). Anchored exactly at the tabulated G0(298.15 K)
# The T-dependent shift comes from standard Cp integration: G(T) = G(Tr) - (T-Tr)·S(Tr) + ΔH(T) - T·ΔS(T)
#       where
#   ΔH(T)=∫Cp dT, ΔS(T)=∫(Cp/T) dT from Tr to T (derived directly from G(T)=H(T)-T·S(T) )
function hp98_gibbs(T, d)
    Tr = 298.15
    ΔH = d.c1*(T-Tr) + 0.5d.c2*(T^2-Tr^2) - d.c3*(1/T-1/Tr) + 2d.c5*(√T-√Tr)
    ΔS = d.c1*log(T/Tr) + d.c2*(T-Tr) - 0.5d.c3*(1/T^2-1/Tr^2) - 2d.c5*(1/√T-1/√Tr)
    return d.G0 - (T-Tr)*d.S + ΔH - T*ΔS
end

# Real, T-dependent equilibrium constants for Huizenga (2001)'s four C-O-H
# reactions, via ΔrG°(T) = -RT·lnK — no more Kᵢ=1 stubs:
#   (1) CO  + 1/2 O2 = CO2      (2) H2  + 1/2 O2 = H2O
#   (3) CH4 + 2 O2   = CO2+2H2O (4) C   + O2     = CO2
# Fugacity coefficients are still assumed ideal (γᵢ=1)
function huizenga_K(T)
    G_H2O = hp98_gibbs(T, HP98.H2O); G_CO2 = hp98_gibbs(T, HP98.CO2)
    G_CO  = hp98_gibbs(T, HP98.CO);  G_CH4 = hp98_gibbs(T, HP98.CH4)
    G_H2  = hp98_gibbs(T, HP98.H2);  G_O2  = hp98_gibbs(T, HP98.O2)
    G_gph = hp98_gibbs(T, HP98.gph)

    ΔG1 = G_CO2 - G_CO - 0.5G_O2
    ΔG2 = G_H2O - G_H2 - 0.5G_O2
    ΔG3 = G_CO2 + 2G_H2O - G_CH4 - 2G_O2
    ΔG4 = G_CO2 - G_gph - G_O2

    return exp(-ΔG1/(R*T)), exp(-ΔG2/(R*T)), exp(-ΔG3/(R*T)), exp(-ΔG4/(R*T))
end

# Shared core: solves the Huizenga (2001) §2.2 quadratic once, returning
# whether a physical graphite-saturated composition exists at this
# fO₂/T/P alongside the mole fractions. `coh_speciation` and `coh_fOₑᵪ`
# both call this rather than each re-deriving (or worse, re-guessing from
# the OTHER one's output) where the validity boundary sits.
function _coh_solve(fO₂, T, P)
    K1, K2, K3, K4 = huizenga_K(T)
    Pbar = 1e4*P                   # P here is GPa (this file's convention throughout); HP98's Kᵢ and
                                    # fO₂ (already log₁₀-bar elsewhere in this file, e.g. carbon_buffer_fO2)
                                    # are referenced to a 1-bar standard state, so Pfluid must be in bar too.
    XCO2 = K4*fO₂/Pbar              # Huizenga eq. 13 (aC=1, γCO2=1)
    XCO  = XCO2/(K1*√fO₂)           # Huizenga eq. 14
    b_H2  = 1/(K2*√fO₂)             # XH2  = b_H2*XH2O          (Huizenga eq. 15)
    b_CH4 = Pbar^2/(K3*fO₂^2)       # XCH4 = b_CH4*XCO2*XH2O^2  (Huizenga eq. 16)

    # Substituting into mass balance ΣXᵢ=1 (Huizenga eq. 9) gives a quadratic
    # in XH2O — re-derived directly here rather than decoding Huizenga's own
    # compressed eq. 17 notation (which reuses "C3" for two different
    # things).
    a_q, b_q, c_q = b_CH4*XCO2, 1+b_H2, XCO2+XCO-1
    disc = b_q^2 - 4a_q*c_q
    # c_q>0 (XCO2+XCO alone already exceed 1) makes BOTH quadratic roots
    # negative even when disc≥0 — checking disc alone isn't enough near the
    # fO2 ceiling for graphite saturation, so check the physical result too.
    #
    # No physical (graphite-saturated) composition exists past this point.
    # The fallback here needs to satisfy two things at once, which two
    # earlier versions each got only one of:
    #   - a TRUE constant (zero dependence on fO2/T/P), so the solver can't
    #     find a spurious self-consistent root sitting inside this region
    #     (an XH2O-only constant let XH2/XCH4 keep varying with fO2 through
    #     b_H2/b_CH4 — not actually flat, so this happened for real);
    #   - a NEUTRAL composition (roughly equal parts H2O/CO2, not "near-pure
    #     CO2 with negligible H2O"), so the H-mass-conservation denominator
    #     in `coh_fOₑᵪ` isn't artificially tiny. A near-zero XH2O sentinel
    #     made the resulting sink swing from 0% (Mh=0) to a solver-breaking
    #     ~1e9 or a flat, Mh-independent cap depending on how it was
    #     special-cased — a discontinuous cliff either way. A neutral
    #     composition instead makes the sink scale smoothly and
    #     continuously with Mh even through this region (still capped
    #     downstream for genuinely large Mh), with no special-casing needed
    #     for Mh=0: it just falls out of the same formula as exactly 0.
    if disc < 0 || c_q >= 0
        return false, oftype(fO₂,0.5), oftype(fO₂,0.5), oftype(fO₂,0.0), oftype(fO₂,0.0), oftype(fO₂,0.0)
    end
    XH2O = (-b_q + √disc)/(2a_q)
    (XH2O < 0) && return false, oftype(fO₂,0.5), oftype(fO₂,0.5), oftype(fO₂,0.0), oftype(fO₂,0.0), oftype(fO₂,0.0)
    XH2  = b_H2*XH2O
    XCH4 = b_CH4*XCO2*XH2O^2
    return true, XH2O, XCO2, XCH4, XCO, XH2
end

# Direct (non-iterative) C-O-H fluid speciation at a GIVEN, externally-set
# fO₂ — Huizenga (2001) §2.2, "externally buffered oxygen fugacity". When no
# physical solution exists (fO2 too high for graphite saturation), returns
# a fixed neutral composition instead — see `_coh_solve` for why.
function coh_speciation(fO₂, T, P)
    _, XH2O, XCO2, XCH4, XCO, XH2 = _coh_solve(fO₂, T, P)
    return XH2O, XCO2, XCH4, XCO, XH2
end

# Mass fraction of the total TOₑₓ budget tied up as O within the speciated
# fluid. The mole fractions above are fully set by (fO₂,T,P) alone — they
# say nothing about how MUCH fluid exists. That's fixed here by hydrogen
# mass conservation against the water budget Mh, the same principle as
# `maxXCO₂_raw` caps carbon: however much fluid exists, its H atoms (2 per
# H2O, 4 per CH4, 2 per H2) can't exceed the H budget Mh actually supplies.
# `logfO₂` is the shared fO2 in the same log₁₀ units used everywhere else
# in this file (`fₛ`/`fₗ`/`fᵪ`). No special-casing for Mh=0 or for the
# invalid-fO2 region: both fall out correctly (0, and a smooth Mh-scaled
# value respectively) from the one formula below, now that the fallback
# composition in `_coh_solve` doesn't create an artificial near-zero
# denominator.
function coh_fOₑᵪ(logfO₂, T, P, Mh, TOₑₓ)
    _, XH2O, XCO2, XCH4, XCO, XH2 = _coh_solve(10^logfO₂, T, P)
    molH_avail       = 2*Mh/mm.H2O
    nH_per_mol_fluid = 2XH2O + 4XCH4 + 2XH2
    n_fluid          = molH_avail/nH_per_mol_fluid
    O_per_mol_fluid  = XH2O + 2XCO2 + XCO
    raw = n_fluid*O_per_mol_fluid*mm.O / TOₑₓ
    # Dividing by TOₑₓ means this can still reach solver-destroying
    # magnitude (~1e9, confirmed) whenever TOₑₓ is small and Mh is large.
    # The other equations here run O(1-100) in log-fO2 units; capping this
    # term keeps it a strong, unmistakable "avoid this" signal in the
    # mass-balance residual without letting its absolute size wreck the LM
    # linear algebra. A capped value >1 already means "this sink alone
    # exceeds the entire TOₑₓ budget" — unphysical regardless of the exact
    # number beyond that point, so clamping loses no real information.
    return min(raw, oftype(raw, 1e3))
end

# ------------------------------------------------------------
# Section 5 — bounded ⇄ unconstrained variable transform
# ------------------------------------------------------------
# The solver operates on y ∈ ℝ; physical quantities are recovered through a
# logistic map into (lo, hi), so a solution can never leave the physically
# valid range and no per-step boundary clipping is needed anywhere.
@inline bounded(y, lo, hi)   = lo + (hi - lo)/(1 + exp(-y))
@inline unbounded(x, lo, hi) = log((x - lo)/(hi - x))

# Finiteness checker and thrower — catches Inf/-Inf too, not just NaN.
# Overflow (Inf) is just as likely as NaN once P/T push Kᵢ(T) or the solid
# table into extreme territory, and it would otherwise pass through
# silently (unlike NaN, `Inf < resnorm` etc. are well-defined comparisons,
# so a plain Inf never trips an error on its own — it just quietly poisons
# the solve). Name kept `nancheck` since every call site already reads that
# way; only what it catches has widened.
@inline function nancheck(label, val, inputs...)
    isfinite(val) || error("$label evaluated to $val (non-finite) — inputs: $inputs")
    return val
end

# ------------------------------------------------------------
# Section 6 — buffered-regime residual (3 unknowns: sOₑₓ, mOₑₓ, XCO₂)
# ------------------------------------------------------------
# The mass-balance equation now always includes the COH fluid sink
# (`coh_fOₑᵪ`), computed directly from the shared fO2 with no extra unknown —
# see Section 4b. When Mh=0 (no water), coh_fOₑᵪ evaluates to exactly 0, so
# this reduces to the pre-COH 3-equation system automatically.
function buffered_residual(y, ctx)
    (; sfO2, T, P, IDV, TOₑₓ, Mf, molfracXB, nB_per_mass, Φ, Mh, lo, hi) = ctx
    sOₑₓ = bounded(y[1], lo[1], hi[1])
    mOₑₓ = bounded(y[2], lo[2], hi[2])
    XCO₂ = bounded(y[3], lo[3], hi[3])

    fₛ = nancheck("sfO2(sOₑₓ)", sfO2(sOₑₓ), sOₑₓ)                                                     # solid fO2, from the MAGEMin table
    fₗ = nancheck("melt_fO2", melt_fO2(mOₑₓ, T, TOₑₓ, Mf, molfracXB, nB_per_mass, IDV), mOₑₓ)          # melt fO2, Hirschmann (2022)
    fᵪ = nancheck("carbon_buffer_fO2", carbon_buffer_fO2(XCO₂, P, T), XCO₂)                            # carbon-buffered fO2, EMOG(D)
    fₑₓ = nancheck("coh_fOₑᵪ", coh_fOₑᵪ(fₛ, T, P, Mh, TOₑₓ), fₛ, Mh)                                          # COH fluid O sink, Huizenga (2001)

    return [ fₛ - fₗ                                    # (1) solid ↔ melt fO2 equilibrium
             fₗ - fᵪ                                    # (2) melt ↔ carbon buffer fO2 equilibrium
             1 - sOₑₓ - mOₑₓ - fₑₓ - Φ*XCO₂/(1 - XCO₂) ] # (3) TOₑₓ mass conservation incl. COH sink
end

# ------------------------------------------------------------
# Section 6b — carbon-saturated ("un-buffered") residual: 2 unknowns
# ------------------------------------------------------------
# If satisfying the EMOG buffer would need more dissolved CO2 than the
# available carbon budget allows, there is no more free carbon left to
# buffer fO2 against: XCO₂ is pinned at its carbon-limited ceiling
# (maxXCO₂) rather than solved for, and the melt↔carbon fO2-equality
# equation no longer applies — only the solid↔melt equilibrium and mass
# balance remain, now with only 2 unknowns (sOₑₓ, mOₑₓ). The COH sink is
# still included in mass balance, same as above.
function carbon_saturated_residual(y, ctx)
    (; sfO2, T, P, IDV, TOₑₓ, Mf, molfracXB, nB_per_mass, Φ, maxXCO₂, Mh, lo, hi) = ctx
    sOₑₓ = bounded(y[1], lo[1], hi[1])
    mOₑₓ = bounded(y[2], lo[2], hi[2])

    fₛ = nancheck("sfO2(sOₑₓ)", sfO2(sOₑₓ), sOₑₓ)
    fₗ = nancheck("melt_fO2", melt_fO2(mOₑₓ, T, TOₑₓ, Mf, molfracXB, nB_per_mass, IDV), mOₑₓ)
    fₑₓ = nancheck("coh_fOₑᵪ", coh_fOₑᵪ(fₛ, T, P, Mh, TOₑₓ), fₛ, Mh)

    return [ fₛ - fₗ                                                 # (1) solid ↔ melt fO2 equilibrium
             1 - sOₑₓ - mOₑₓ - fₑₓ - Φ*maxXCO₂/(1 - maxXCO₂) ]       # (2) mass conservation incl. COH sink, XCO₂ fixed at cap
end

# ------------------------------------------------------------
# Section 7 — generic Levenberg–Marquardt solver
# ------------------------------------------------------------
# Generic in the residual's dimension, so it serves both the 3-equation
# (buffered) and 2-equation (carbon-saturated) systems unchanged. λ grows
# when a step fails to reduce the residual (falling back toward gradient descent,
# which stays well-posed even when the Jacobian is ill-conditioned near a
# saturated bound) and shrinks when a step succeeds (recovering fast
# Newton-like convergence away from any boundary).
function lm_solve(F::Function, y0::AbstractVector{<:Real}; maxiter=500, tol=1e-12, λ0=1e-3, λmax=1e10, λmin=1e-14, debug=false)
    y, Fy   = copy(y0), F(y0)
    resnorm = norm(Fy)
    λ       = λ0
    n       = length(y0)
    Id      = Matrix{Float64}(I, n, n)
    it      = 0
    while resnorm > tol && it < maxiter
        it += 1
        J  = ForwardDiff.jacobian(F, y)
        dy = (J'J + λ*Id) \ (J'Fy)     # Levenberg–Marquardt normal equations
        ytrial  = y .- dy
        Ftrial  = F(ytrial)
        newnorm = norm(Ftrial)
        if debug && (it<=10 || it%50==0)
            println("  it=$it:  y=$y Fy=$Fy ")
            println("           Fy=$Fy")
            println("           (resnorm=$resnorm | λ=$λ)")
            println("           J = $J")
            println("           dy = $dy -> newnorm would be $newnorm \n")
        end
        if newnorm < resnorm
            y, Fy, resnorm = ytrial, Ftrial, newnorm
            λ = max(λ/3, λmin)          # step accepted: trust the (near-)Newton direction more
        else
            λ = min(λ*5, λmax)          # step rejected: fall back toward gradient descent
        end
    end
    return (converged = resnorm<=tol, y=y, resnorm=resnorm, iters=it)
end

# ------------------------------------------------------------
# Section 8 — top-level driver
# ------------------------------------------------------------
# Deliberately narrower than misc.jl's `partition_Oₑₓ` for now: no Fe³⁺/Feᵀ
# ratio inputs, no evolution plotting/diagnostics. The COH fluid sink is
# now real thermodynamics (Section 4b: Huizenga 2001 method, Holland &
# Powell 1998 Kᵢ(T)), not a placeholder — and, because it needs no extra
# unknown of its own (see Section 4b), it's simply always included in
# whichever of the two regimes below actually applies:
#   buffered (3 eq)          — carbon free
#   carbon_saturated (2 eq)  — carbon exhausted, XCO₂ fixed at its ceiling
function partition_Oₑₓ(P::Real, T::Real, p::Real, ϕ::Real, TOex::Real, TC::Real, H::Real;
                           nr=25, Rmax=0.18, maxiter=500, tol=1e-12, verbose=false, data=nothing, out=false, debug=false)

    (; X, XB, Xox, SymXox, mmXox) = endmember_composition(p)
    Mt = 1e3   # reference cell mass (1 kg, expressed in grams to match mm.* in g/mol)
    (; Mf, Ms, Mh, TOₑₓ, molfracXB, nB_per_mass, maxXCO₂_raw, Φ) = mass_budget(Mt, ϕ, TC, TOex, H, XB, mmXox)

    if TOₑₓ == 0.0
        verbose && println("Total Oₑₓ budget = 0.0 — nothing to partition.")
        return (converged=true, regime=:none, sOₑₓ=0.0, mOₑₓ=0.0, XCO₂=0.0, iters=0, resnorm=0.0)
    end

    IDV = nancheck("solve_∫ΔVdP", solve_∫ΔVdP([P-0.05P, P, P+0.05P], [T-0.05T, T, T+0.05T])[2,2,1], P, T)

    # Compositional ceiling for the solid
    molX    = (X./mmXox) ./ sum(X./mmXox)
    maxsOₑₓ = min((0.5molX[3]*mm.O)/(sum(X) - 0.5molX[3]*mm.O)*(Ms/TOₑₓ), 1.0)

    # Compositional ceiling for the melt
    nOxideTotal = Mf*nB_per_mass
    maxmOₑₓ = min(0.5*molfracXB[3]*nOxideTotal*mm.O/TOₑₓ, 1.0)
    maxXCO₂ = clamp(maxXCO₂_raw, 0.0, 1.0)
    lowclip = 1e-10

    sfO2, sOₑₓ_table_max = solid_fO2_interpolant(P, T, X, Xox, SymXox, Ms, TOₑₓ; nr=nr, Rmax=Rmax, data=data)
    maxsOₑₓ = min(maxsOₑₓ, sOₑₓ_table_max)   # never search past what MAGEMin actually tabulated

    lo  = (lowclip, lowclip, lowclip)
    hi  = (maxsOₑₓ - lowclip, maxmOₑₓ - lowclip, maxXCO₂ - lowclip)
    ctx = (; sfO2, T, P, IDV, TOₑₓ, Mf, molfracXB, nB_per_mass, Φ, Mh, lo, hi)

    y0 = [unbounded(5e-4*hi[1], lo[1], hi[1]),   # initial guesses: mostly-reduced solid,
          unbounded(0.5*hi[2],  lo[2], hi[2]),   # half-saturated melt,
          unbounded(0.01*hi[3], lo[3], hi[3])]   # small melt XCO₂

    res = lm_solve(y -> buffered_residual(y, ctx), y0; maxiter=maxiter, tol=tol, debug=debug)

    sOₑₓ = bounded(res.y[1], lo[1], hi[1])
    mOₑₓ = bounded(res.y[2], lo[2], hi[2])
    XCO₂ = bounded(res.y[3], lo[3], hi[3])
    regime = :buffered

    # If the buffered solve fails with XCO₂ pinned at its ceiling, the
    # carbon budget — not the Newton/LM iteration — is the limiting factor:
    # retry in the carbon-saturated regime (XCO₂ fixed at maxXCO₂, melt↔
    # carbon equation dropped, only sOₑₓ/mOₑₓ left to solve for; the COH
    # sink stays active in the mass balance either way).
    ctx2 = nothing
    if !res.converged && XCO₂ > 0.999*hi[3]
        regime = :carbon_saturated
        lo2 = (lo[1], lo[2]); hi2 = (hi[1], hi[2])
        ctx2 = (; sfO2, T, P, IDV, TOₑₓ, Mf, molfracXB, nB_per_mass, Φ, maxXCO₂, Mh, lo=lo2, hi=hi2)
        y0b  = [unbounded(sOₑₓ, lo2[1], hi2[1]), unbounded(mOₑₓ, lo2[2], hi2[2])]
        res  = lm_solve(y -> carbon_saturated_residual(y, ctx2), y0b; maxiter=maxiter, tol=tol, debug=debug)
        sOₑₓ = bounded(res.y[1], lo2[1], hi2[1])
        mOₑₓ = bounded(res.y[2], lo2[2], hi2[2])
        XCO₂ = maxXCO₂
    end

    eqfO2 = sfO2(sOₑₓ)   # solid/melt/carbon-buffer fO2 agree at convergence — report the shared value
    cohsink = coh_fOₑᵪ(eqfO2, T, P, Mh, TOₑₓ)

    if verbose
        Fres = regime==:buffered ? buffered_residual(res.y, ctx) : carbon_saturated_residual(res.y, ctx2)
        XH2O, XCO2f, XCH4, XCOf, XH2 = coh_speciation(10^eqfO2, T, P)
        println("========================================================================================================")
        println("==== Solution (P=$(P) GPa | T=$(T) K | ϕ=$(1e2ϕ) % | TOex=$(1e2TOex) wt% | mix=$p | TC=$(1e2TC) wt% | H=$(1e2H) wt%) ====")
        println("========================================================================================================")
        println("   • Regime: $regime")
        println("   • Converged: $(res.converged)  ($(res.iters) iterations with |F| = $(res.resnorm))")
        println("   • Equilibrium log fO₂ = $eqfO2")
        println("   • Melt XCO₂ = $XCO₂")
        println(" -----------------------------")
        if regime==:buffered
            println("   • Residuals: (1) solid↔melt = $(Fres[1]), (2) melt↔carbon = $(Fres[2]), (3) mass balance incl. COH = $(Fres[3])")
        else
            println("   • Residuals → (1) solid↔melt = $(Fres[1]), (2) mass balance incl. COH [XCO₂ fixed at cap] = $(Fres[2])")
        end
        println(" -----------------------------")
        println("   • Lever upper bounds: sOₑₓ ≤ $maxsOₑₓ (MAGEMin table reached up to $sOₑₓ_table_max), mOₑₓ ≤ $maxmOₑₓ, XCO₂ ≤ $maxXCO₂")
        println(" -----------------------------")
        println("   • Distances → sOₑₓ = $(round(1e2sOₑₓ/maxsOₑₓ, digits=2))%, mOₑₓ = $(round(1e2mOₑₓ/maxmOₑₓ, digits=2))%, XCO₂ = $(round(1e2XCO₂/maxXCO₂, digits=2))% of their respective caps")
        println(" -----------------------------")
        println("   • COH fluid (mole fractions, graphite-saturated, ideal mixing) → XH₂O=$(round(XH2O,digits=4)), XCO₂=$(round(XCO2f,digits=4)), XCH₄=$(round(XCH4,digits=6)), XCO=$(round(XCOf,digits=6)), XH₂=$(round(XH2,digits=6))")
        println("   • COH Oₑₓ sink = $(round(1e2cohsink, digits=4))% of TOₑₓ")
        println(" =============================")
        # A post-processed, subtractive quantity like "melt CO2" can go
        # wildly negative when the solve hasn't converged and cohsink is
        # sitting at its cap — that's not a new number worth displaying,
        # just fallout from an already-flagged `Converged: false`. Print
        # "?" instead of a nonsensical percentage in that case.
        meltCO2 = 1 - sOₑₓ - mOₑₓ - cohsink
        meltCO2_str = (isfinite(meltCO2) && 0.0 <= meltCO2 <= 1.0) ? "$(round(1e2meltCO2, digits=4))%" : "?"
        println("   → TOₑₓ partitioning → [$(round(1e2sOₑₓ, digits=4))% solid, $(round(1e2mOₑₓ, digits=4))% melt as Fe₂O₃, $meltCO2_str melt CO₂, $(round(1e2cohsink, digits=4))% COH fluid]\n")
    end

    if out; return (converged=res.converged, regime=regime, sOₑₓ=sOₑₓ, mOₑₓ=mOₑₓ, XCO₂=XCO₂, cohsink=cohsink, iters=res.iters, resnorm=res.resnorm); end

end
