# TOex — mathematical procedures by event

Each section answers to an event or arrow of [TOex_overview.md](TOex_overview.md). Notation:
`f = log₁₀ fO₂`; masses in grams on a 1 kg reference cell unless stated; `Mₓ` molar masses;
`ϕ` melt fraction, `φᵢ` fugacity coefficients. Fortran names are in `volatiles_mod.f90` unless
another file is given; Julia names are in `src/TOex.jl`.

Contents

- S.1 Reference state and the three ledgers
- S.2 Inputs and gram budgets
- S.3 The one equation, monotonicity, tolerance
- S.4 Solid iron: the tabulated oxybarometer
- S.5 Melt iron and the sink capacities
- S.6 Solid carbon: EMOG buffer and the graphite–carbonate ramp
- S.7 Carbon activity aC
- S.8 Melt CO₂ and its ceiling
- S.9 The C-O-H fluid: speciation, EOS, pressure limit
- S.10 Melt hydrogen: H₂O, H₂ and the hydrogen split
- S.11 The non-ideality loop and its convergence test
- S.12 Bracket, root rule, warm start, failure
- S.13 Reporting: regime, buffers, fluid density
- S.14 Gather, scatter, and the fluid window
- S.15 Solid water storage and water activity
- S.16 Volatile-enhanced melting
- S.17 Hydrous density of solid and melt
- S.18 Darcy transport of the fluid triple
- S.19 Speciated degassing
- S.20 Surface box A and B: Θ, speciation, burial
- S.21 Hydrogen escape
- S.22 The B → A return
- S.23 Ingassing gates
- S.24 Conservation laws and the drift test
- S.25 Constants and sources
- S.26 Root finders: bisection, regula falsi, Illinois

---

## S.1 Reference state and the three ledgers

**Event.** Every budget, demand and transfer in the model.

<img src="md_figs/S01_ledger.png" alt="S.1" width="1100">

*Figure S.1. The ledger of eq. 1.1: O_ex carried per atom of C, H or Fe by each species; positive entries are sinks, negative ones sources.*

The state variables are three extensive amounts per reservoir: n_H (carried as H₂O-equivalent
mass), n_C, and excess oxygen

```julia
O_ex ≡ n_O − n_H/2   (1.1)   # O_ex counts oxygen beyond what H₂O would use, so Fe as FeO, C as graphite and H as water all score zero.
```

measured against a fully reduced assemblage: Fe as FeO, C as graphite or diamond, H as H₂O. Per
carrier this gives CO₂ and carbonate +2 per C, CO +1, graphite 0, CH₄ −2, H₂O 0, H₂ −1 per H₂,
Fe³⁺ +½ per Fe, Fe⁰ −1 per Fe. The negative entries make the fluid and the solid *signed* sinks
that can act as sources. Because H₂O is the reference, the ocean holds no O_ex and water exchange
between ocean and mantle cannot change the total (S.24). Every atom-conserving reaction is the
identity on (n_C, n_H, O_ex); only transport across a boundary moves O_ex.

Real fluid mass is never stored: it is derived from the speciation as `mfluid_per_H`, the ratio
of fluid mass to H₂O-equivalent hydrogen (1 for pure water). Two consumers use it, the Darcy flux
(S.18) and degassing (S.19).

Code: `toex_ref_state`, `TOEX_*` slot constants; Julia `mass_budget`.

## S.2 Inputs and gram budgets

**Event.** `Gram budgets (S.2)`.

<img src="md_figs/S02_budgets.png" alt="S.2" width="1100">

*Figure S.2. The two endmembers per oxide as mass fractions, basalt drawn upward and harzburgite downward from the zero line; the cell mixes them by p (eq. 2.2).*

Inputs: P [GPa] (Pbar = 10⁴ P), T [K], basalt fraction p, melt fraction ϕ, and the tracer mass
fractions TOex, TC, H. On the reference cell Mt = 10³ g:

```julia
Mf = ϕ Mt,  Ms = (1−ϕ) Mt,  Mc = TC Mt,  Mh = H Mt,  TOₑₓ = TOex Mt   (2.1)   # Gram masses of fluid, solid, carbon, hydrogen and excess oxygen from the cell mass Mt and its fractions.
```

Bulk oxides by linear mixing of the two endmembers (EM_Xbs, EM_Xhz, eight oxides including Na₂O
and Cr₂O₃, which sb24 needs):

```julia
X = p XB + (1−p) XH                                                                 (2.2)   # Cell oxide composition as the basalt-fraction p blend of the two endmembers.
nB_per_mass = Σⱼ XBⱼ/Mⱼ,  molfracXBⱼ = (XBⱼ/Mⱼ)/nB_per_mass                         (2.3)   # Moles of oxide per gram of basalt and the basalt oxide mole fractions.
n_ox = Mf nB_per_mass  (melt oxide moles),  nC_total = Mc/M_C,  molH = 2 Mh/M_H2O   (2.4)   # Melt oxide moles, total carbon moles and hydrogen atoms from the gram budgets.
Φ = 2 M_O n_ox   grams of O per unit XCO₂/(1−XCO₂)                                  (2.5)   # Conversion factor from the melt CO₂ ratio to grams of oxygen.
```

Everything downstream is absolute grams; the budget appears only on the right-hand side of the
residual (S.3). A host converts once by mass, grams × Mtot/Mt. TOₑₓ = 0 is a state, not an
absence (Fe³⁺ balanced by CH₄/H₂ or Fe⁰) and is solved like any other budget.

Code: `toex_build_ctx`; Julia `mass_budget`, `endmember_composition`.

## S.3 The one equation, monotonicity, tolerance

**Event.** `Residual F(f) (S.3)`, `Capacities, tolerance (S.3, S.5)`.

<img src="md_figs/S03_residual.png" alt="S.3" width="1100">

*Figure S.3. The oxygen each phase holds at one shared log10 fO₂, stacked, for a wet, carbon-bearing, partly molten cell. Any single phase could hold the budget on its own at some f (dotted: the solid alone at f = −1.1), but all phases read the same f, so the root is where the whole stack meets the budget (f = −2.6).*

```julia
F(f) ≡ sOₑₓ(f) + mOₑₓ(f) + carbonO(f) + carbO(f) + fₑₓ(f) + mH2ₑₓ(f) − TOₑₓ = 0   (3.1)   # The oxygen residual: every phase read at the shared f, summed, minus the budget; the root is F = 0.
```

with sOₑₓ, mOₑₓ the ferric oxygen of solid and melt (S.4, S.5), carbonO the melt CO₂ (S.8),
carbO solid carbonate (S.6), fₑₓ the net fluid oxygen (S.9), mH2ₑₓ ≤ 0 the melt H₂ (S.10).
Posing the problem this way makes fO₂ equality an identity: an unsaturated sink *is* the loading
whose oxybarometer reads f. It can only fail where a sink is pinned at a physical end.

Monotonicity gives existence and uniqueness: sOₑₓ, mOₑₓ are inverses of increasing oxybarometers,
melt CO₂ rises with 10^f, carbonate is a logistic in f, and the fluid shifts oxidised as f rises.
F is non-decreasing, so plain bisection suffices and no boolean may appear inside any demand (a
step of 2 O per carbon lets F straddle zero without a root; hence the ramp of S.6). Where the
melt-H₂ source and carrier switches make F non-monotone, S.12 states the root rule.

Tolerance, grams need a scale and the budget supplies it:

```julia
atol = max( rel · max(|TOₑₓ|, 10⁻⁶ capacity),  10⁻¹² (Mc + Mh + capacity) )   (3.2)   # Absolute tolerance on F, scaled by the budget or the sink capacity, with a gram-level floor.
```

The absolute floor is the inner solves' own gram-level noise; without it a zero budget cannot
converge. Carbon and hydrogen close exactly at every trial f; the outer test is on oxygen alone.

Code: `toex_residual`; Julia `unified_residual`.

## S.4 Solid iron: the tabulated oxybarometer

**Event.** `Solid oxybarometer from table (S.4)`, `Solid Fe3+ demand (S.4, S.5)`.

<img src="md_figs/S04_solid_oxybarometer.png" alt="S.4" width="1100">

*Figure S.4. The tabulated oxybarometer at four (P, T): log10 fO₂ against the signed ferric ratio (metal branch to the left of the FeO reference), and the native-iron fraction of the metal branch. The 60 GPa curve shows metal coexisting with Fe³⁺ (disproportionation).*

The solid's fO₂ at fixed (P, T, p) is a function of its ferric ratio R = Fe³⁺/ΣFe. MAGEMin (sb24,
Stixrude & Lithgow-Bertelloni 2024, corrected local build) is sampled once on a signed R axis and
tabulated; the run reads the table.

Table: `StagfO2.dat`, 100 P × 100 T × 36 R per sector (upper mantle, transition zone, lower
mantle) and composition (ten basalt fractions from 0 to 1); R from −1.96 to 0.95 as 30 asinh-spaced nodes to 0.4 plus 6 logit-spaced
nodes (0.538, 0.671, 0.781, 0.861, 0.916, 0.95). Negative R is the metal branch. `StagMetal.dat`,
same grid, gives native Fe as a weight fraction of the solid. The interpolation coordinate is
u = asinh(R/R₀), R₀ = 2·10⁻³; the R axis is written explicitly in the header and never rebuilt
as an even grid. Metal is clamped in R at both ends; fO₂ is not.

Re-express the axis as oxygen held by the solid:

```julia
sOwtᵢ = X_O(rᵢ) (0.5 rᵢ)/(1 + 0.5 rᵢ),   gridᵢ = sOwtᵢ Ms   [g]   (4.1)   # Grams of excess oxygen the solid holds at table node i, from its ferric ratio r_i.
maxsOₑₓ = 0.5 (Ms X_FeO/M_FeO) M_O   (all Fe²⁺ → Fe³⁺)            (4.2)   # Solid capacity: the oxygen needed to make every Fe²⁺ ferric.
minsOₑₓ = X_Fe(r = −2) (M_O/M_Fe) Ms  (all Fe²⁺ → Fe⁰)            (4.3)   # Solid floor: the oxygen given up when every Fe²⁺ becomes metal.
```

Repair and extrapolation of the sampled curve: isotonic regression (least-squares monotone fit,
pool-adjacent-violators) removes dips, a 10⁻⁶ tilt makes it strictly increasing, and beyond the
last nodes the oxybarometer's own log-linear form in Fe³⁺/Fe²⁺ continues it, slope 1/0.1917,
written R/(1−R) so it diverges at the stoichiometric limit, not at the sampling edge. Below the
signed table's most reduced node the reading is flat: Fe⁰ + ½O₂ = FeO buffers at a finite fO₂.

Demand at a trial f: invert the monotone table for the loading whose reading is f. A cell between
two tabulated compositions shares one fO₂ between them, so its ferric ratio at f is the blend of the
two columns' own ratios weighted by mass fraction and FeO, read directly with no root find (the
reads, the silica-saturation step and their accuracy are in `TOex_numerics.md`, 2.4 and 2.5). If at
the metal floor the rock already reads above f, pin sOₑₓ = −minsOₑₓ (all Fe as Fe⁰); if at the
ferric cap it still reads below f, pin sOₑₓ = +maxsOₑₓ. Iron speciation at the answer:

```julia
FeM_solid = metal_wt Ms/M_Fe / nFe,   FeR_solid = max(0, 2(nO + nFe⁰)/nFe)   (4.4)   # Reported metal fraction and ferric ratio of the solid at the root.
```

and the Fe balance is checked so Fe³⁺/ΣFe never exceeds 1 (audit REQ 14).

Why the axis matters: the oxybarometer is log-linear in R/(1−R); uniform-in-R sampling cannot
reach 0.01 log units at any node count, asinh spacing does at 30. The 2026-09-11 tables from the
stock sb24 build were ~2 log units too oxidised in the lower mantle (merged Fe site entropy);
the corrected build and the extended axis replaced them. Accuracy on [0.4, 0.95]: 0.095 log units
mean, 0.42 worst.

Code: `thermochem_int('fO2'|'metal')` and `thermochem_fO2_columns` in `thermochem_mod.f90`,
`toex_solid_fO2_table` in `volatiles_bind.f90`, `toex_fill_sfO2`, the blend in `toex_derive`; Julia `solid_fO2_interpolant`, `Rhybridspace` in
`sf_H2O.jl`. Generator: `Perple_X_MAGEMin/MAGEMin/s_fO2_dens_metal.jl`.

## S.5 Melt iron and the sink capacities

**Event.** `Melt Fe3+ demand (S.5)`, the "pinned at a cap?" arrows.

<img src="md_figs/S05_melt_oxybarometer.png" alt="S.5" width="1100">

*Figure S.5. The Hirschmann (2022) melt oxybarometer at four (P, T) and the pressure term IDV(P, T) of eq. 5.4.*

Melt oxybarometer, Hirschmann (2022) eq. 21 with the Deng et al. (2020) volume term:

```julia
n_ox = Mf nB_per_mass,  nO⁺ = mOₑₓ/M_O,  nFeO = X_FeO n_ox                       (5.1)   # Melt oxide moles, moles of excess oxygen in the melt and the melt FeO moles.
log₁₀(X_FeO₁.₅/X_FeO) = a f + b + c/T − (ΔCp/(R ln10))(1 − T₀/T − ln(T/T₀))
                        − IDV(P,T)/(RT ln10) + (1/T) Σᵢ yᵢ Xᵢ  (+ cross terms)   (5.2)   # Hirschmann (2022) melt oxybarometer: ferric/ferrous ratio from f, T, P and composition.
X_FeO₁.₅/X_FeO = 2 nO⁺/(nFeO − 2 nO⁺)     (one O⁺ oxidises two Fe²⁺)             (5.3)   # The same ratio written in excess-oxygen moles, since one O⁺ oxidises two Fe²⁺.
IDV(P,T) = ∫₁ᵇᵃʳ^P ΔV dP',  ΔV = ½[V(Mg₁₄Fe₂Si₁₆O₄₉) − V(Mg₁₄Fe₂Si₁₆O₄₈)] > 0    (5.4)   # Pressure term of the melt oxybarometer: integrated volume change of the ferric reaction.
maxmOₑₓ = 0.5 X_FeO n_ox M_O                                                     (5.5)   # Melt capacity: the oxygen needed to make every FeO ferric.
```

Xᵢ on a single-cation basis (X_AlO₁.₅ = 2 n_Al₂O₃/N_cat, etc.); a, b, c, ΔCp, T₀ = 1673 K and
yᵢ from Hirschmann (2022) Table 2; ΔV from the BM4 fits of Deng et al. (2020), positive to
140 GPa. IDV is a 281 × 41 (P, T) table built once per rank (cubic Hermite in P, linear in T);
Julia reads the same table. Anchor: peridotite at 1 bar, 2173 K, Fe³⁺/ΣFe 0.0605 gives −5.966.

Capacities and pinning, the branch table:

| condition | consequence |
|---|---|
| solid at metal floor, sfO₂(−minsOₑₓ) ≥ f | sOₑₓ = −minsOₑₓ, its own reading legitimately ≠ f |
| solid at ferric cap, sfO₂(maxsOₑₓ) ≤ f | sOₑₓ = +maxsOₑₓ, sink full |
| melt floor / cap | same, no metal branch |
| both Fe caps | the rock's capacity binds the whole budget |
| phase absent, max ≤ 0 (ϕ = 0 or 1) | no oxybarometer exists; its fO₂ is suppressed everywhere |

Capacities are stoichiometric and deliberately not clipped at the budget: that clip is the mass
balance itself, already the equation being solved.

Code: `toex_melt_fO2`, `toex_idv_table`; Julia `melt_fO2`, `_invert_demand`.

## S.6 Solid carbon: EMOG buffer and the graphite–carbonate ramp

**Event.** `aC ceiling from EMOG (S.6, S.7)`, `aC at ceiling, precipitate C (S.6)`.

<img src="md_figs/S06_emog_ramp.png" alt="S.6" width="1100">

*Figure S.6. Left: the EMOG(D) buffer against T at five pressures, absolute log10 fO₂. Middle: the same buffer relative to IW, with FMQ − IW dashed, so its place in the buffer landscape is read directly. Right: the graphite–carbonate ramp of eq. 6.6 and the O_ex per solid carbon it implies.*

The solid host of carbon is graphite/diamond below the enstatite-magnesite-olivine-graphite
(EMOG, or EMOD with diamond) buffer and carbonate above it. Gibbs energies from Holland & Powell
(1998):

```julia
ΔH = c₁(T−Tr) + ½c₂(T²−Tr²) − c₃(1/T−1/Tr) + 2c₅(√T−√Tr)                                    (6.1)   # Enthalpy increment of a species from its heat-capacity polynomial between Tr and T.
ΔS = c₁ ln(T/Tr) + c₂(T−Tr) − ½c₃(1/T²−1/Tr²) − 2c₅(1/√T−1/√Tr)                             (6.2)   # Entropy increment of the same species from the same polynomial.
G(T) = G₀ − (T−Tr)S + ΔH − TΔS,   ΔG_P = V₀(T) K/(K'−1) [(1 + K'(Pbar−1)/K)^(1−1/K') − 1]   (6.3)   # Gibbs energy at T and 1 bar, then its pressure correction from a Murnaghan volume integral.
μ_C = min[G_graphite, G_diamond]                                                            (6.4)   # The carbon standard state is whichever of graphite or diamond is stable.
2 MgCO₃ + Mg₂Si₂O₆ = 2 Mg₂SiO₄ + 2 C + 2 O₂,   fE = −ΔG/(2 RT ln10)                         (6.5)   # The EMOG reaction and its equilibrium log fO₂, fE, from ΔG.
```

The host change is ramped, never switched, because a jump of 2 O per carbon would put a step in
F(f):

```julia
x_carb = 1/(1 + exp[−(f − fE)/w]),  w = 0.05 log units                               (6.6)   # Carbonate fraction of the solid carbon as a smooth logistic across the buffer.
nC_carb = x_carb nC_solid,  nC_gph = (1 − x_carb) nC_solid,  carbO = 2 nC_carb M_O   (6.7)   # Solid carbon split into carbonate and graphite, and the oxygen the carbonate holds.
```

nC_solid is what melt and fluid did not take, `max(0, nC_total − nC_melt − nC_fluid)`, with the
max as the ramp (undersaturated it is max(0, 0) because S.7 just made demand equal the budget).
The melt holds no graphite: precipitated reduced carbon goes to the solid slot, and to the melt
slot only when there is no solid.

Code: `toex_reaction_Ks`, `toex_logfO2_EMOG`; Julia `reaction_Ks`, `logfO2_EMOG`.

## S.7 Carbon activity aC

**Event.** `Solve aC (S.7)` and the ceiling arrows.

<img src="md_figs/S07_carbon_activity.png" alt="S.7" width="1100">

*Figure S.7. One trial f, two budgets. The demand D(aC) of melt plus fluid (black; parts in red and teal) rises linearly with aC and diverges where the fluid can hold no more CO₂; the EMOG cap sits at 1 but the highest valid activity is 0.01 (step 1). (a) A small budget meets D below that ceiling: the solve returns that aC and no solid carbon exists. (b) A large budget lies above D at the ceiling: aC is pinned there and the difference precipitates, carbonate or graphite by the ramp of S.6.*

aC is the activity that makes melt CO₂, fluid carbon and solid carbon consistent at a given f; it
is solved inside every oxygen residual. Thermodynamics fixes the ratios of the carbon species at
graphite saturation; aC scales them all down when the budget cannot sustain saturation.

```julia
aC_cap = min(1, 10^(fE − f))       thermodynamic ceiling from EMOG(D)   (7.1)   # Carbon activity can never exceed the value at which graphite or carbonate precipitates.
D(aC) = nC_melt(aC) + nC_fluid(aC)  what melt and fluid demand          (7.2)   # Total carbon the melt and fluid would take at activity aC.
```

Procedure at each trial f:

1. Is the speciation (S.9) valid at aC_cap? If not, walk the ceiling down a decade at a time to
   the first valid activity, then bisect (60 halvings) between it and the last invalid one to
   the validity edge, aC_use. The cap is an upper bound from the buffer, not a guarantee that a
   fluid exists there; bracketing against an invalid endpoint returned garbage. The bisection is
   needed because the fluid's carbon demand diverges at the edge (X_H₂O → 0, n_fluid = molH/n_H);
   stopping at the decade point declared the cell saturated with carbonate beside a water-rich
   fluid below the decarbonation edge (EMOG above CCO), where it should hold no carbonate.
2. Is D(aC_use) ≤ nC_total? Then the cell is saturated: sit at aC_use and let S.6 precipitate the
   remainder. Testing the cap first *is* the decision between the two physical states.
3. Otherwise solve D(aC) = nC_total for log₁₀ aC ∈ [−300, log₁₀ aC_use] with a bracketed
   Illinois (regula falsi, S.26) step on a log residual, log₁₀(D/nC_total): at small aC demand is
   proportional to aC, so log-log is a straight line and ~16 evaluations suffice. The floor is
   −300 because the melt's buffer coefficient can exceed 10⁴⁰ at high f.

D is monotone in aC (XCO₂, XCO ∝ aC; XCH₄ through aC; melt CO₂ through the softcap), which makes
the bracket legal without repair. An invalid speciation returns D = +∞, read as "demand too
high", which keeps D monotone and pushes aC down. With no hydrogen in the fluid D is the melt
term alone, and the hydrogen-free fluid of S.7b is placed after the solve.

Above fE, aC is a formal bridging variable, the graphite-referenced carbon fugacity a
carbonate-buffered rock is consistent with, which lets one set of equations serve on both sides
of the buffer.

Code: `toex_solve_aC`, `toex_carbon_demand`, `toex_illinois`; Julia `_solve_aC`,
`_carbon_demand`.

## S.7b The hydrogen-free CO₂–CO fluid

**Event.** `Solve aC (S.7)`, cells whose fluid holds no hydrogen (dry, or all H in the melt).

Hydrogen fixes how much C–O–H fluid there is (n_fluid = molH_fl/n_H, 9.3). With molH_fl = 0 the
fluid has no carrier, but a CO₂–CO fluid can still exist: its composition follows from aC and f
alone (9.1, XCO₂ and XCO ∝ aC), and closure fixes the activity at which it is a pure phase,

```julia
aC* = 1/(K1 fO₂/(φ_CO2 P) + K2 √fO₂/(φ_CO P)),   XCO₂ = aC* K1 fO₂/(φ_CO2 P),   XCO = 1 − XCO₂   (7b.1)   # Closure activity of a CO₂–CO fluid; φ iterated with ZD09 on the binary.
```

Like water and ice, the phases do not coexist away from the boundary. If aC* < aC_cap the fluid is
stable before any solid saturates: the melt takes its share at aC* and the fluid takes the rest,
nC_fluid = nC_total − nC_melt(aC*), with n_fluid = nC_fluid (one C per molecule). If aC* > aC_cap
the solid saturates first and S.7 stands. aC* = aC_cap is the validity edge of S.7 (above EMOG,
the decarbonation en + mag = fo + CO₂) and, below EMOG, CCO (aC* = 1). With hydrogen present the
fluid is diluted and coexists with carbonate deeper than the edge; without it there is no such
field. The switch is ramped like EMOG (S.6) so F(f) stays continuous,

```julia
nC_fluid = (nC_total − nC_melt(aC_u)) / (1 + exp(−log₁₀(aC_cap/aC*)/w_CO)),   aC_u = min(aC*, aC_cap),   w_CO = 0.05   (7b.2)   # Fluid share of the carbon the melt leaves, ramped across the edge.
```

the remainder going to solid carbon on the S.6 ramp. If the melt holds all carbon below aC*, no
fluid forms. Not applied above fluid_Pmax (S.9). The fluid carries 2XCO₂ + XCO excess oxygen per
carbon (S.1) and no hydrogen, so it passes the window of S.14 at h = 0. It cannot migrate (S.18
needs a hydrogen carrier), stays on its tracers and is re-solved in place, degasses on eruption
(S.19), and leaves the solid's water activity scaling off (S.15, no water in the fluid).

Effect: dry carbon-bearing cells shallower than the edge are buffered at CCO by graphite + fluid
instead of EMOG by graphite + carbonate, 0.3 to 1.6 log units more reduced at 0.5 to 2.5 GPa;
oxidised cells exchange carbonate for CO₂ fluid at unchanged fO₂; deeper than the edge nothing
moves. Wet cells with hydrogen in the fluid are bit-identical.

Code: `toex_co_fluid` and the branch in `toex_coh_state`; Julia `_co_fluid`, `coh_state`.

## S.8 Melt CO₂ and its ceiling

**Event.** Part of the residual (carbonO) and of D(aC).

<img src="md_figs/S08_melt_co2.png" alt="S.8" width="1100">

*Figure S.8. Graphite-saturated melt X_CO₂ against fO₂ at four (P, T), with the soft cap and the EMOG position at each; above EMOG aC falls below 1 and the real X_CO₂ leaves this curve.*

Graphite-saturated CO₂ in silicate melt, a carbon buffer inverted for XCO₂ (Stagno & Frost 2010
form):

```julia
f = 5.44 − 21380/T + 0.078 (10⁴P − 1)/T + log₁₀ XCO₂                                       (8.1)   # Melt CO₂ solubility law linking f, P, T and the melt CO₂ mole fraction.
B_melt(f) = 10^[f − 5.44 + 21380/T − 0.078(10⁴P − 1)/T]                                    (8.2)   # The same law inverted: melt CO₂ ratio per unit carbon activity at this f.
X_m = softcap(aC B_melt, X_max),  nC_melt = X_m n_ox/(1 − X_m),  carbonO = 2 M_O nC_melt   (8.3)   # Melt CO₂ mole fraction with a soft ceiling, the carbon moles it holds and the oxygen it holds.
```

X_max = `toex_XCO2max`, default 0.5: X is a mole fraction on (oxide moles + CO₂), and 0.5 is
stoichiometric dolomite, 47.7 wt% CO₂, the most a carbonate melt can carry. It never engages on
silicate melt (a few wt% CO₂ to 11 GPa, X ≲ 0.05) and binds only on the carbonatite of S.16
(40–45 wt% CO₂, X ≈ 0.42–0.47). The soft cap (S.25) avoids a corner. XCO₂ is gated on the melt
existing; ungated it painted a melt-free lithosphere at the cap. The divisor (1 − X_m) is why the
parameter is clamped into (10⁻⁶, 0.95).

Sources: Stagno & Frost (2010); Stagno et al. (2013); Dasgupta & Hirschmann (2006, 2007) for the
carbonatite range; Brooker & Kjarsgaard (2011) for the miscibility gap below 2.5 GPa; Eguchi &
Dasgupta (2018) as the candidate replacement for (8.1).

Code: `toex_carbon_buffer`; Julia `carbon_buffer_fO2`, `invert_carbon_buffer_fO2`.

## S.9 The C-O-H fluid: speciation, EOS, pressure limit

**Event.** `Speciate C-O-H fluid (S.9)`, `Refresh phi, ZD09 EOS (S.11)`, `Fluid = pure water,
no EOS (S.9)`.

<img src="md_figs/S09_fluid_speciation.png" alt="S.9" width="1100">

*Figure S.9. Left: C-O-H speciation with ZD09 coefficients at 3 GPa, 1500 K along fO₂, at the highest valid aC below the EMOG cap: CH₄ + H₂ below IW, H₂O between, CO₂ + H₂O above EMOG. Right: pure-species fugacity coefficients against P at 1600 K, and the extrapolation beyond 10 GPa.*

Species CH₄, H₂O, CO₂, H₂, CO (C₂H₆ and O₂ carried by the EOS, not by the speciation). With
K₁..K₄ from S.6's Gibbs energies (C + O₂ = CO₂, C + ½O₂ = CO, H₂ + ½O₂ = H₂O, C + 2H₂ = CH₄):

```julia
XCO₂ = aC K₁ fO₂/(φ_CO2 Pbar),   XCO = aC K₂ √fO₂/(φ_CO Pbar)                (9.1)   # Fluid CO₂ and CO mole fractions from carbon activity, fO₂ and the fugacity coefficients.
γ = 1 − XCO₂ − XCO,   α XH₂² + β XH₂ − γ = 0,   XH₂ = 2γ/(β + √(β² + 4αγ))   (9.2)   # The hydrogen species close the sum to one; XH₂ is the positive root of the quadratic.
```

with α from the CH₄ term (aC K₄ φ_H2² Pbar/φ_CH4) and β from the H₂O term (K₃ φ_H2 √fO₂/φ_H2O);
XH₂O and XCH₄ follow and ΣX = 1 identically. The stable form 2γ/(β + √disc) is used because the
textbook root divides by a coefficient that vanishes as aC → 0. Invalid if γ ≤ 0, the
discriminant is negative, XH₂ < 0, or any X is non-finite; the last test is not redundant because
NaN passes every comparison.

Amounts, set by hydrogen and never by carbon:

```julia
nH = 2XH₂O + 4XCH₄ + 2XH₂,   n_fluid = molH_fl/nH         (9.3)   # Hydrogen atoms per fluid molecule, hence the number of fluid moles from the hydrogen it owns.
nC_fluid = n_fluid (XCO₂ + XCH₄ + XCO)                    (9.4)   # Carbon moles carried by the fluid.
fₑₓ = n_fluid (2XCO₂ + XCO − 2XCH₄ − XH₂) M_O    signed   (9.5)   # Net excess oxygen of the fluid in grams: CO₂ and CO count positive, CH₄ and H₂ negative.
```

A hydrogen-free C-O fluid therefore cannot exist (assumption I.3): with water off the carbon
sits in graphite, carbonate and melt. fₑₓ is gated on n_fluid > 0 and never clamped when
non-finite; the gross flows O_up, O_rel are reported separately because a reducing fluid grows
the O_ex inventory and that source is booked against hydrogen.

Fugacity coefficients and fluid density from the Zhang & Duan (2009) EOS:

```julia
ε = Σᵢⱼ xᵢxⱼ K1ᵢⱼ √(εᵢεⱼ),  K1[H₂O,CO₂] = 0.85, K1[CH₄,H₂O] = 0.8, else 1          (9.6)   # Mixture energy parameter of the ZD09 EOS with the two non-unity binary corrections.
Pm = 3.06258806161 σ³ Pbar/ε,  Tm = 154 T/ε;  solve R Tm Z/V = Pm for V            (9.7)   # Reduced pressure and temperature of the mixture, then the volume root of the EOS.
ln γᵢ = Z − 1 − ln Z + S₁ − 2S₂(1 − εᵢ/ε) + 6(1 − Z)(1 − σᵢ/σ),  φᵢ = exp(ln γᵢ)   (9.8)   # Fugacity coefficient of each species from the compressibility and the mixing sums.
ρ_fluid = 1000 M̄/V   [kg m⁻³], the same root                                       (9.9)   # Fluid density from the mean molar mass and the same molar volume.
```

The volume root is bracketed from a seed and warm-started per cell. ZD09 is calibrated to 10 GPa
and 2573 K and is extrapolated to the cell pressure: φ_H2O reaches 10³³ at the CMB, which is the
physics of ln φ ≈ ∫V dP/RT and is balanced against an equally unclamped solid fO₂. The molar
volume error above 10 GPa is +8 % at 10 GPa, +31 % at 100; the root fails for CO above ~77 GPa
and for CH₄/H₂ above ~100. Cells in 10–30 GPa and 30–Pmax are counted per step.

Above `toex_fluid_Pmax` = 50 GPa the fluid is restricted to pure water: no speciation, no EOS,
fₑₓ = 0, nC_fluid = 0, the melt at a_H2O = 1 without H₂, the rest as neutral free water (AQUA
has no H₂O fluid above ~50–60 GPa at mantle T). The Darcy density above Pmax comes from the
pure-water ZD09 root, which never fails to 135 GPa.

Sources: Zhang & Duan (2009) GCA 73:2089; Holland & Powell (1998); Frost & McCammon (2008) for
the buffer landscape the results are checked against.

Code: `toex_coh_reaction_step`, `toex_coh_state`, `toex_zd_fugacity`, `zd_tables`; Julia
`_coh_reaction_step`, `coh_state`, `zd_fugacity`.

## S.10 Melt hydrogen: H₂O, H₂ and the hydrogen split

**Event.** `Hydrogen split melt/fluid (S.10)`, the "melt wants more H than owned?" arrow.

<img src="md_figs/S10_melt_hydrogen.png" alt="S.10" width="1100">

*Figure S.10. Left: melt water from Sun & Yao (2026) against the fluid water activity at four (P, T), with the legacy Katz law dashed. Right: melt H₂ from Hirschmann et al. (2012) against f_H2 and its soft cap.*

Melt H₂O, Sun & Yao (2026) eq. 7 with molecular and hydroxyl water, at the fluid's own P_H2O:

```julia
ln S_H2O = a₀ ln P_H2O + a₁ ln T + a₂ P²/T + (a₃ √P + a₄) X_Na2O/T     (10.1)   # Molecular-water solubility of the melt from water pressure, T and composition.
ln S_OH  = a₅ ln P_H2O + a₆ ln T + (a₇ + a₈ X_Al2O3 + a₉ X_Na2O) P/T   (10.2)   # Hydroxyl solubility of the melt from the same variables.
n_H(melt) = 2 · 10⁻² Mf (S_H2O + S_OH)/M_H2O                           (10.3)   # Hydrogen atoms the melt holds as water at saturation, from the two solubilities.
```

with P_H2O = a_H2O P through the fluid's XH₂O and φ_H2O, so the law's own exponents carry the
activity dependence. The water half is calibrated to 2 GPa; P is clamped there at fixed X_H2O,
so the capacity stops rising with depth (41.7 wt% at a_H2O = 1), the conservative direction.

Melt H₂, Hirschmann et al. (2012) eq. 6, a sixth oxygen sink:

```julia
ln K = ln K₀ + (ΔH/R)(1/T₀ − 1/T) − ΔV (P − P₀)/(RT),  K = X_H2/f_H2,  ΔH = 0       (10.4)   # Temperature and pressure dependence of the melt H₂ solubility constant.
X_H2 = softcap(K f_H2, X_H2max),  n_H2 = X_H2 n_ox/(1 − X_H2),  mH2ₑₓ = −n_H2 M_O   (10.5)   # Melt H₂ mole fraction with a soft ceiling, its moles, and the negative oxygen it counts as.
```

f_H2 = X_H2(fluid) φ_H2 Pbar. This term is the first sensitive to the absolute fO₂ scale (through
H₂ + ½O₂ = H₂O); the internal IW sits ~1.5 log units above the literature IW at 1 bar, which
understates dissolved H₂ by ~5×.

The split, one bracketed Illinois solve (S.26) per pass, branch chosen by a test:

```julia
SATURATED       molH_fl + n_H(melt, λ = 1) = molH_tot,   bracket [0, molH_tot]         (10.6)   # Saturated case: the fluid gets the hydrogen the water-saturated melt cannot hold.
UNDERSATURATED  n_H(melt, λ) = molH_tot, molH_fl = 0,     bracket log₁₀ λ ∈ [−60, 0]   (10.7)   # Undersaturated case: no fluid, the melt holds all hydrogen at a reduced water activity λ.
```

The test is n_H(melt, λ = 1) ≥ molH_tot: the melt's saturated demand already exceeds the
inventory, so no fluid can exist and λ scales both fugacities together (f_H2O/f_H2 = K₃√fO₂
untouched). (10.6) stays linear because molH = 0 is an endpoint log space cannot represent. If
the Illinois core runs out of steps with the bracket intact its last iterate is the answer; only
a non-finite residual falls back, to the melt's saturated content with the fluid taking the rest.
Melt hydrogen is clamped to the cell budget on the converged state only. Melt H₂O carries zero
O_ex and enters the hydrogen balance only; melt H₂ is the sink.

Sources: Sun & Yao (2026); Hirschmann et al. (2012); Katz et al. (2003) with Rozel (2018) as the
legacy `toex_melt_H2O_law = 'katz'`.

Code: `toex_melt_h2o_wt`, `toex_melt_h2`, `toex_solve_hydrogen`; Julia `_solve_hydrogen`,
`melt_H2O_sunyao`.

## S.11 The non-ideality loop and its convergence test

**Event.** `Refresh phi, ZD09 EOS (S.11)` and the "phi, C, H settled?" arrow.

<img src="md_figs/S11_ideality_loop.png" alt="S.11" width="1100">

*Figure S.11. Distance of the fluid speciation from the fully converged answer after k passes of the loop, at three fO₂; geometric convergence and the 10⁻⁸ tolerance.*

Around S.7, S.9, S.10 at fixed f:

```julia
(a) hydrogen split (10.6)/(10.7)   (b) aC (S.7)   (c) speciate (9.1)–(9.2)            # The hydrogen species close the sum to one; XH₂ is the positive root of the quadratic.
(d) x ← x/Σx, φ ← ZD09 (9.6)–(9.8), the one expensive call                   (11.1)   # One non-ideality pass: three cheap sub-solves, then the one costly fugacity refresh.
```

Why a loop at all: the mass-action equations (9.1)–(9.2) give the composition x for a given set
of coefficients φ, and the EOS gives the coefficients φ(x) for a given composition; each is exact
for its input, but the x computed with the previous φ is not the x whose φ was used, so the pair
is not self-consistent until x = x(φ(x)). The loop is successive substitution on that fixed
point: solve x at the current φ, evaluate φ at the new x, repeat until φ (and with it x) stops
moving. φ is frozen during the aC solve so the EOS is evaluated once per pass rather than once per aC
trial (~100× fewer EOS calls). The loop runs `toex_ideality_min` = 2 to `toex_ideality_max` = 30
passes and exits when the joint residual, |Σ fC − 1|, the hydrogen balance and Δφ, is below
`toex_ideality_tol` = 10⁻⁸. Checking φ alone is not enough: the C and H residuals were measured
at the old φ and the loop alternates them. The first pass starts from φ = 1 (ideal). Leaving
φ = 1 would solve every equilibrium as an ideal gas at GPa pressure: at f = −8, 1600 K, 5 GPa it
moves XCO₂ by a factor 1200. Convergence is geometric at ~10⁻² per pass. Only finite positive
φᵢ are accepted; a failed EOS evaluation keeps the previous value.

Closure test of the unified engine (`toex_unified`, since 2026-11-08). The loop's final joint
residual is kept as `res_ideality` = max(r_C, r_H), the relative closure of the carbon and hydrogen
balances on the state the cell is written back with. A cell counts as converged only if the outer
root converged and res_ideality ≤ `toex_closure_tol` = 10⁻⁶. Before, a cell whose hydrogen root had
no sign change (the floor of S.15) was reported converged with r_H of order one. The log line
`TOex: H floor ... closure failures` counts the cells that fail on this test alone.

Code: the loop in `toex_coh_state`, `toex_coh_state_u`; Julia `coh_state`, `coh_state_u`.

## S.12 Bracket, root rule, warm start, failure

**Event.** `Bracket log10 fO2 (S.12)`, `Bracket step (S.12)`, `FAILED: keep last partition
(S.12)`, the "previous fO2 usable?" arrow.

<img src="md_figs/S12_warmstart.png" alt="S.12" width="1100">

*Figure S.12. Left: the residual on the outer bracket [−120, 40] with the root. Right: a ±3 log-unit zoom around the root with the ±0.25 and ±2 warm-start windows.*

Outer bracket [−120, +40] in f, wider than any physical value on purpose: every demand
saturates outside its range, so the bracket only has to straddle the root. No sign change means
no root anywhere: the cell is reported FAILED, the closer endpoint is returned, the caller keeps
the cell's previous partition, counts it and never stops the run. A collapsed bracket with |F|
still large means F is discontinuous there and is also reported as failure, never success. The
unified engine adds two causes: a hydrogen floor in a cell with melt (S.15, eq. 15.3) and a C or H
closure above `toex_closure_tol` (S.11). Both keep the previous partition in the same way.

Step: bisection first, then Illinois once bracketed (S.26), ≤ `toex_maxiter_O2` = 500 iterations. The
outer method stays bracketed because the solid demand is a table, piecewise linear in u, whose
derivative is a step function that can cycle a bare Newton.

Warm start and the root rule. F is not monotone everywhere: in wet, carbon-bearing, molten cells
the melt-H₂ source grows as f falls and the fluid switches carrier (CH₄ ↔ CO₂, graphite in or
out), so two or three discrete roots ~1 log unit apart can satisfy |F| ≤ atol. The rule is
continuity in time: the solve first tries a window of ± `toex_warmstart_narrow` = 0.25 around the
cell's previous fO₂, then ± `toex_warmstart_window` = 2.0, then the full bracket; the reported
root is the one nearest last step's. The previous fO₂ is the mass-weighted fO2POS of the cell's
tracers written back by the solve; the table fills fO2POS only while it is 0 (fresh tracers).
The known cost is that a vanished root goes straight to the full bracket. Runs are therefore
reproducible in time but not at step 1 from a cold start.

Code: `toex_partition`, `toex_solve_cells` (window ladder, counters `narrow/wide/full`); Julia
`bisect_solve`, `safe_root`.

## S.13 Reporting: regime, buffers, fluid density

**Event.** `State at root, regime label (S.13)`.

<img src="md_figs/S13_buffers_profile.png" alt="S.13" width="1100">

*Figure S.13. FMQ, IW and EMOG(D) relative to IW along the mean geotherm of the local run at 2.91 Gyr (core cooling on, corrected tables), against its horizontal-mean fO₂ below the lid.*

The state is re-derived once at the converged f, so what is reported is exactly what the residual
saw. Iron speciation from (4.4) and (5.3); the fluid's mole fractions, M̄ and ρ from (9.9);
`mfluid_per_H`; the melt's water and H₂. The buffer offsets are reported only: with
`fO2buffer_write` on they are evaluated on output steps from each cell's solved T, P and fO₂
(evaluating them every step cost ~23 % of the solve); off by default, never evaluated:

```julia
FMQ: 3 Fe₂SiO₄ + O₂ = 3 SiO₂ + 2 Fe₃O₄,   IW: 4 Fe + 2 O₂ = Fe₄O₄  (SLB 2024, with the Landau
and Hillert–Jarl magnetic terms)                                                                (13.1)   # The two reference buffers from the SLB 2024 endmember Gibbs energies, reported only.
```

dFMQ is meaningless at low T (at 305 K the buffer sits at −80). The regime label is a label
only, tested in order: capacity_exceeded, oxidant_exhausted, solid/melt saturated or depleted,
carbonate_graphite_buffered, carbonate_saturated, carbon_undersaturated, buffered. `fO2_equal`
and `pinned` say whether the identity of S.3 held and, if not, whether it was allowed to fail.
Shares are quoted against the pool actually available (TOₑₓ plus what a reducing fluid released)
so every entry is positive and sums to 100 %.

Validation: the Julia reference and the Fortran port agree to 2·10⁻⁹ log units on 53 cases
(`TOex_bench.txt`).

Code: `toex_report`, `toex_buffers`; Julia `classify_regime`, `partition_Oₑₓ` report block.

## S.14 Gather, scatter, and the fluid window

**Event.** `Gather cell budgets (S.14)`, `Scatter to tracers (S.14)`, `React in transit (S.14)`.

<img src="md_figs/S14_fluid_window.png" alt="S.14" width="1100">

*Figure S.14. The admissible (n_C, n_Oex) region of a parcel with 1 mol H: between CH₄/H₂ below and CO₂ + H₂O above.*

Gather: the cell's tracers are summed into ten slots w₀ (solid O_ex, oxidised C, reduced C; melt
O_ex, oxidised C, reduced C; fluid H, C, O_ex; melt bound H), with masses, melt fraction, basalt
fraction, and the mass-weighted previous fO₂. Composition is fixed per tracer type,
X = p EM_Xbs + (1 − p) EM_Xhz.

Scatter, two rules: rock slots are updated *incrementally* so sub-grid redox heterogeneity
survives (each tracer keeps its offset from the cell mean), fluid slots *proportionally* to the
tracer's hydrogen. Oxidised and reduced carbon are renormalised to conserve carbon per tracer,
so carbon conservation is a property of the writeback, and `toex_fC_dev` is the diagnostic that
shows a wrong split. No O_ex slot may be floored: solid metal, fluid CH₄/H₂ and melt H₂ each make
a legitimately negative entry, and a floor manufactures oxygen one-signed every step.

Solid O_ex is scattered in two components, iron (by the tracer's Fe capacity) and the rest, which
is what stopped a spurious O_ex pile at the CMB (audit REQ 10).

The fluid validity window, asserted after every merge and every subdomain receive:

```julia
−n_H/2 ≤ n_Oex ≤ +2 n_C   (14.1)   # Physical window of the fluid triple: no more reduced than H₂, no more oxidised than CO₂.
```

Upper bound: all C as CO₂ and H no more oxidised than H₂O. Lower bound: CH₄ and H₂ are both
−½ per H, so the most reducing fluid is a function of hydrogen alone. Both bounds are linear and
homogeneous, so the window is closed under summation and a violating sum requires a violating
term; the assertion localises a fault. Position inside the window is the redox state, with no
thermodynamics. A cell whose solution cannot be written back (window violated, or hydrogen that
no host can hold) is stranded to the rock (`toex_strand_unrealisable`).

React in transit (`toex_react_in_transit`, default on): a Darcy parcel entering a cell is
included in that cell's gather and the cell is re-solved, the chromatographic local-equilibrium
picture; without it a passing parcel exchanges nothing with a saturated cell.

Code: `toex_cell_gather`, `toex_cell_scatter`, `toex_apply_slot`, `toex_shares_to_masses`,
`check_window`, `toex_strand_unrealisable`.

## S.15 Solid water storage and water activity

**Event.** `Solid storage S * a_H2O^n (S.15)`, the "water fluid last step?" arrow.

<img src="md_figs/S15_storage.png" alt="S.15" width="1100">

*Figure S.15. Left: S_eff/S against the fluid water activity for three exponents, with the 10⁻³ floor. Right: why the update is implicit, an explicit lag at gain 100 limit-cycles while the implicit step relaxes.*

Storage capacity S(P, T, x_bs) is read from `StagH2O.dat` (the solid's water at saturation from
`scripts/sH2O/sfO2.jl`, ten compositions per sector, linear in x_bs between neighbouring blocks;
sH2O_math.md). It has no water axis. S is saturation against pure
water; when the coexisting fluid is not pure water the effective capacity is (Keppler)

```julia
S_eff = S a_H2O^n,   n = toex_aH2O_n = 1   (15.1)   # Effective solid water capacity, scaled by the fluid water activity to the power n.
```

with a_H2O = XH₂O of the fluid from the last solve, floored at `toex_aH2O_min` = 10⁻³. The update
is implicit because the solid is a high-gain reservoir: the analytic gain of dissolving into the
solid, G = n W (1 − a)/f_H2O (W the cell's bound water, f_H2O the fluid's water), is ~10² and an
explicit lag limit-cycles. The tracer's S is moved by (S_new − w)/(1 + G) toward the target.
Molten tracers in a cell whose fluid carries water take the solved melt content as their S, so
the Darcy sweep of the same step can neither refill nor shed the melt at a_H2O = 1 (S.10). A
minus-ulp melt water is floored before it becomes a negative S.

Code: `toex_solve_cells` post-pass (`toex_aH2O_scaling`), `eval_fO2_traWsol` in `tracers.f90`,
`thermochem_int('solubility')`.

**Unified engine: the solid inside the hydrogen root.** With `toex_unified` the solid is a hydrogen
sink of the hydrogen balance itself (S.10), at the trial f_H2. Its capacity C(P, T, p, a) is read
from `StagH2Oa.dat` at 16 water-activity nodes from 1 down to 10⁻³, log-log between nodes, with the
last slope continued below 10⁻³. The bound water is shared among the cell's solid tracers by their
own capacities at the solved activity a_solid (`toex_apply_bound`).

Hydrogen floor (2026-11-08). Hydrous phases can hold C(a) on a plateau as a falls, so the solid
may need more water at a = 10⁻³ than the cell owns, and the continued last slope does not bring it
down to the budget. In the logged cell (0.27 GPa, 575 K, p = 0.06, no melt) C(10⁻³) = 0.26 wt%
and C(10⁻²⁰) = 0.22 wt% against a water budget of 0.003 wt%. The hydrogen balance has no root there.
The solver returned its lower bracket end and reported a split of 48.7 (solid) and −47.7 (melt)
times the budget as converged. The write-back min(m_H(solid), H_tot − m_H(fluid)) capped the
solid at the budget, so no water was created, but the reported partition was wrong. The floor is
flagged when the root ends on its lower bracket end with the demand there still above the budget:

```julia
floor, Mf = 0:  n_H(solid) = molH_tot − n_H(fluid),   a_solid = 10⁻³ (lowest node)   (15.2)   # No melt: the solid takes what the fluid leaves, closing the balance exactly.
floor, Mf > 0:  cell FAILED, previous partition kept (S.12)                          (15.3)   # With melt the solid–melt split below the table is undetermined.
```

Without melt the solid is the only condensed sink, so (15.2) is exact in mass. Its activity lies
below the table and is not resolved; the lowest node stands in for it, also as the weight of
`toex_apply_bound`. With melt any split satisfies the balance, so the cell is not solved. The floor
occurs in the upper-mantle sector below about 600 to 810 K at P ≤ 2 GPa, and in the cold corners of
the transition-zone and lower-mantle sectors (20 to 60 GPa, 350 to 800 K) where water-rich hydrous
phases are stable. In a 12-step continuation of the newcode2 frame 600 it fired in one cell per
step and never with melt.

Code: `toex_u_xy` (flag `toex_u_hfloor`), `toex_coh_state_u` (fields `h_floor`, `h_floor_melt`),
`toex_solve_cells` (counters); Julia `coh_state_u` (`h_floor`), `partition_Oₑₓ` (`closure_tol`).

## S.16 Volatile-enhanced melting

**Event.** `Melt / erupt on volatile solidus (S.16)` and its two arrows.

<img src="md_figs/S16_volatile_melting_laws.png" alt="S.16" width="1100">

*Figure S.16. The melting law as coded: Katz water depression, carbonated solidi with ledge and taper against data, the carbonatite plateau, and the two switches.*

One switch, `volatile_enhanced_melting`, each agent acting only where its element is present.

Water, Katz, Spiegelman & Langmuir (2003) eq. 16–18, with the melt fraction:

```julia
ΔT_Katz = K X_melt^γ,   X_melt = c_H2O/(D + F(1 − D)),   X_melt ≤ X_sat(P) = χ₁P^λ + χ₂P   (16.1)   # Katz water depression of the solidus from the melt water content, capped at saturation.
```

K = 43, γ = 0.75, D = 0.01, χ₁ = 12, λ = 0.6, χ₂ = 1; c_H2O the bound water, F the cell's melt
fraction from the previous equilibrium (one-pass seed). The dilution term is the brake: at F = 0
it is Katz's first drop, at F = 1 the bulk. Above 7 GPa the depression is damped,
0.15 + 0.85/(1 + (P − 7)/6.8)⁵, our choice, not Katz's; irrelevant at the model's 0.02–0.04 wt%.
The liquidus is lowered by the same law.

Carbonate, Dasgupta & Hirschmann (2010) compiled peridotite solidus, held at 973 °C below 3 GPa:

```julia
T_cs,per: 973 °C at 3 GPa, 1116 at 5, 1197 at 7, 1274 at 10, 1393 at 15, 1477 at 20, 1550 at 25
          (piecewise linear)                                                                      (16.2)   # Carbonated peridotite solidus nodes in P from the Dasgupta & Hirschmann compilation.
T_cs,ecl: 1050 + 36.4 (P − 3) °C to 8.5 GPa, then 1250 °C flat                                    (16.3)   # Carbonated eclogite solidus, linear to 8.5 GPa then held flat.
```

blended from peridotite at x_bs ≤ bs_ref to eclogite at x_bs = 1. The ledge below ~1.9 GPa, where
carbonate is not a phase in lherzolite and CO₂ is a solute:

```julia
u = clamp((P − P_ledge)/w_ledge, 0, 1),  s_l = u²(3 − 2u)                                (16.4)   # Pressure switch: smoothstep from 0 below the carbonate ledge to 1 above it.
s_l ← s_l [1 − smoothstep((P − P_max)/w)],  P_max = 25, w = 7 GPa                        (16.5)   # The same switch tapered back to 0 above the deepest solidus data.
S_c = ½(1 + tanh((log₁₀ c_oxC − log₁₀ C_min)/δ)),  C_min = 10⁻³ wt% C, δ = 0.3 decades   (16.6)   # Concentration switch: tanh ramp in log carbonate content around the dust threshold.
D_c = S_c [ s_l (T_dry − T_cs) + (1 − s_l) 25 K/wt% min(C_CO₂, 2 wt%) ]                  (16.7)   # Carbonate depression: the eutectic drop where the switch is on, the solute drop where it is off.
```

P_ledge = 1.9, w_ledge = 0.6 GPa. (16.5) tapers the eutectic out above the deepest data (32 GPa);
beyond it only the solute term remains, and the physical cutoff at depth is redox (carbonate
reduces to diamond), carried by S_c on c_oxC, not by the solidus.

The eutectic plateau and where the silicate solidus sits:

```julia
F_carb = S_c s_l c_oxC/w_C,   w_C = 12 wt% C                                                       (16.8)    # Carbonatite melt fraction: all carbonate carbon into a melt with 12 wt% C.
T_cs,eff = max( T_dry − D_c − ΔT_Katz,  min(T_sil,w, T_floor) ),  T_floor = 1173 K                 (16.9)    # Effective carbonatite solidus with water added, floored at the volatile-saturated solidus.
T_sil,w = T_dry − ΔT_Katz,   T_sil = max( T_sil,w − S_c s_l ΔT_sil,  T_cs,eff ),  ΔT_sil = 100 K   (16.10)   # Wet silicate solidus, lowered up to 100 K by carbonate, never below the carbonatite solidus.
T_melting(f) = T_cs,eff + (T_sil − T_cs,eff) f/F_carb   for f < F_carb                             (16.11)   # Melting temperature between the two solidi over the carbonatite melt fraction.
```

```julia
lever rule (tested, not ported): T_melting(f) = T_std(f; T_sil,w) − (T_sil,w − T_cs,eff) F_carb/f   for f ≥ F_carb   (16.12)   # CO₂ diluted into the melt lowers the silicate melting point by the eutectic drop scaled by X_C/12 wt%.
```

so 0.5 wt% carbonate carbon gives 4 % carbonatite, not 30 % silicate melt, and carbonate lowers
silicate melting by ~100 K (Dasgupta et al. 2013). The 100 K is a step in c_oxC (S_c is a switch at
10⁻³ wt%), so the silicate productivity above T_sil is the same at 0.03 and 0.5 wt% C. Eq. 16.12
replaces it by the binary lever rule with a liquidus linear in T from (T_cs, 12 wt% C) to
(T_sil,w, 0): f = F_carb (T_sil,w − T_cs)/(T_sil,w − T). It scales with content, but at 15 GPa
T_sil,w − T_cs ≈ 600 K, so 140 K below the dry solidus it gives f = 4 F_carb, more than the shift;
a 20-step continuation of the 2.5 Gyr frame raised the 350–800 km mean melt fraction from 1.6·10⁻³
to 5.2·10⁻³ (overview fig. 8). The lever rule is in a scratch build behind `carb_lever_rule`
(default off) and is not ported. What actually sets the transition-zone melting is the mantle
sitting within ~100 K of the dry solidus there with carbonate present; the liquidus shape between
T_cs and T_sil,w at 10–25 GPa is unmeasured. With c_oxC = 0 the law reduces to Katz bit for
bit. Redox melting and freezing emerge: the law reads c_oxC only, and TOex moves carbon across
EMOG by fO₂ one step earlier.

Measured (2.5 Gyr frame): the constrained law removes the 1000–2500 km eutectic tail and halves
700–1000 km with 0–700 km, the CMB pile and cost unchanged. The 300–700 km melt is the correct
response to the carbon the model places there; a deep oxidised mantle in older frames came from
the pre-correction tables (S.4).

Sources: Dasgupta & Hirschmann (2006, 2007, 2010); Dasgupta et al. (2004, 2007, 2013); Yaxley &
Brey (2004); Kiseeva et al. (2012, 2013); Thomson et al. (2016); Ghosh et al. (2009); Litasov &
Ohtani (2009); Wyllie & Huang (1976); Eggler (1978); Wallace & Green (1988); Falloon & Green
(1989); Hirschmann (2010); Katz et al. (2003); Stagno et al. (2013); Rohrbach & Schmidt (2011).

Code: `volatile_solidus`, `carbonated_solidus`, `Katz2003_*` in `melting_mod.f90`; `Tmelting` in
`melting.f90`. Figure: `TOex_overview.md` fig. 4.

## S.17 Hydrous density of solid and melt

**Event.** `Tracers to grid, densities (S.17)`, `water_dependent_rho_law = 'minimize'`.

<img src="md_figs/S17_rheology_laws.png" alt="S.17" width="1100">

*Figure S.17. The rheological hooks: Karato viscosity law, wet yield envelope, hydrous density tables of solid and melt with a 1-bar literature anchor, and the yield factors.*

Two tables, `StagSdens.dat` and `StagMdens.dat`, nP × nT × nH per endmember, holding the percent
density decrease of rock and of silicate melt relative to the dry material at the same (P, T):

```julia
ρ = ρ_ref(z, C, f) [1 − φ_s Δρ_S(P, T, c_s) − φ_m Δρ_M(P, T, c_m)]   (17.1)   # Bulk density: reference density reduced by the hydrous solid and melt density decreases.
```

φ_s, φ_m the solid and melt mass fractions, c_s the solid's water (bulk minus the melt's share),
c_m the melt's water. First-order form of 1/ρ = Σ wᵢ/ρᵢ, error below 0.3 %.

Solid, ideal mixing with a molar volume of water (Gerya et al. 2004, condensed branch X_liq = 1):

```julia
V = (1 − x) V_A + x V_H2O,  V_A = M_A/ρ_A (MAGEMin sb24),  Δρ_S = 1 − [M_H/V]/ρ_A                   (17.2)   # Solid hydrous density by ideal volume mixing of the anhydrous rock with a water volume.
V_H2O = 10 [V_s + c₁ ΔV_s1 e₁/(1 − e₁)] [(P₀ + φ)/(P + φ)]^(1/5),
        e₁ = exp[−(ΔH_s1 + ΔV_s1 Ψ)/RT],  Ψ = (5/4)(P₀ + φ)^(1/5)[(P + φ)^(4/5) − (P₀ + φ)^(4/5)]   (17.3)   # Molar volume of water in the solid, Gerya et al. (2004) form with its pressure term.
```

Melt, Jing & Karato (2011, 2012) deformable hard-sphere mixture:

```julia
σᵢ = σᵢ⁰ (T/Tᵣ)^ηᵢ exp[ξᵢ α₀ (T − Tᵣ)/3] (V/V₀)^(ξᵢ/3),  V_m = (N_A π/6) Σ Xᵢ σᵢ³,  f = V_m/V           (17.4)   # Species diameters and packing fraction of the Jing & Karato (2011) hard-sphere melt EOS.
P = (RT/V) [ (1 − ξ) Φ(f) − Φ(f₀)(V₀/V)^(1/3) + ξ Φ(f₀)(V_m⁰/V_m)^(5/3) ],  Φ = (1 + f + f²)/(1 − f)³   (17.5)   # The hard-sphere equation of state solved for the melt volume at P and T.
```

Δρ_M = 1 − ρ(X_hyd)/ρ(X_dry), floored at 0 (`melt_floor`): the law's H₂O sphere collapses
under pressure and makes hydrous melt denser than dry above ~10 GPa, which nothing supports
(FPMD water in melt 12.4 / 10.5 / 8.2 cm³/mol at 10 / 20 / 50 GPa, Karki et al. 2021). The
ideal-mixing option (`melt_H2O_model = :ideal`) mixes the water in with the FPMD volume instead.

Alternative law: `water_dependent_rho_law = 'linear'` uses upstream's a0rho_water in
`calc_dens` (no tables).

Sources: Jing & Karato (2011, 2012); Lange & Carmichael (1987); Ochs & Lange (1999); Gerya et al.
(2004); Karki, Ghosh & Karato (2021); Mookherjee et al. (2008); Bajgain et al. (2015); Matsukage
et al. (2005); Sakamaki et al. (2006); Hemley et al. (1987); Stixrude & Lithgow-Bertelloni (2024);
Gerya, Bardi, Karato & Murakami (2026).

Code: `thermochem_int('Sdens'|'Mdens')`, `fx` density branch in `refstate.f90`; generator
`PT_H2O_ρ` in `misc.jl`, `s_fO2_dens_metal.jl`.

## S.18 Darcy transport of the fluid triple

**Event.** Diagram 3.3, every node.

<img src="md_figs/S18_darcy.png" alt="S.18" width="1100">

*Figure S.18. Left: ZD09 fluid density along 1500 K for pure water and two mixtures to the CMB, with pure water checked against IAPWS-95 at 1000, 1500 and 2000 K (dots; the standard is valid to 1 GPa and 1273 K and is extrapolated in the shaded range): the two agree to a few percent to 10 GPa. Right: Darcy velocity per unit face area for the three permeability sectors.*

The sweep marches each column bottom to top. The parcel is a 3-vector (H as H₂O-equivalent, C,
O_ex signed), exchanged between cells and subdomains as such.

Two transfer rules, which one applies is never a judgement call:

- **Dissolution**, pool → bound water: hydrogen only, up to S_eff − w per tracer (S.15). It takes
  zero C and, H₂O being the reference, zero O_ex, so the ledger closes with no extra term; the
  pool's C/H rises, which is correct. Feasibility: only hydrogen already in H₂O form can leave,
  dn_H ≤ n_H + 2 n_Oex; as the pool dries of H₂O, a_H2O → 0 and S_eff → 0 stop it. With carbon on
  the bound is `h2o_form_H()`.
- **Phase split**, pool → {stays, ascends}: congruent, fixed proportions on all three. Darcy moves
  a phase; nothing fractionates elements at a permeability limit.

Exsolution: a tracer above WaterSatCutoff × S_eff sheds its free triple first (the mobile phase),
then bound-derived water as pure H₂O (dC = dO_ex = 0).

Flux law per face, the capacity for this dt:

```julia
allowed_flow = k A g Δρ ρ_fluid dt/η_fluid   (18.1)   # Darcy mass allowed across a face per step from permeability, buoyancy and fluid viscosity.
```

k the sectored constant permeability (`limitflow_constant_permeability`, three sectors), or
k = k₀ (φ_t/φ₀)ⁿ with n = 3 under `limitflow_relperm`, φ_w ≤ 1, φ_weff ≤ 1 − totmelt, and a
percolation threshold; ρ_fluid and η_fluid from the cell's solve (ZD09 density, S.9; IAPWS
pure-water viscosity on the H₂O-equivalent, the weaker half). Masses in the flux law are real
masses through `mfluid_per_H`. What reaches the top is limited by the minimum face capacity on
the path, right for quasi-steady flow.

Stranding: a parcel that has lost its hydrogen carrier cannot migrate; its C and O_ex stay in
the rock of the cell where that happened (in transit). At the surface (erupting tracers, the
ingas-depth release) the triple goes to the surface reservoir instead (S.19, S.20). Free fluid
that a cell cannot hold at the end of the sweep is stranded the same way, so the next step starts
realisable.

Top of column: the whole triple leaves the domain, H to the ocean (water_budget 5), C and O_ex
to box A (carbon_budget 3 and 2, oxygen_budget 3).

The bound-water diffusion path (`smooth_tracer_H2O`) moves pure H₂O and needs no O_ex term.

Sources: Connolly & Podladchikov (2012) for the permeability form (constant sectors used
instead); Rudge (2018) for the Kozeny–Carman exponent; IAPWS (2008).

Code: `WaterSaturationAndVerticalRedistribution`, `parcel_to_cell`, `cell_to_parcel`,
`toex_strand_fluid`, `toex_degas_fluid`, `check_window` in `volatileprocs_mod.f90`.

## S.19 Speciated degassing

**Event.** `Speciated degassing (S.19)`, `Degassed C, O_ex, H into A (S.19)`.

<img src="md_figs/S19_degassing.png" alt="S.19" width="1100">

*Figure S.19. The CO₂/CH₄ endmember split of a degassed parcel against its O_ex/C.*

An erupting tracer releases its melt CO₂ (fraction `degas_C_melt` = 1), its free-fluid carbon
(`degas_C_fluid` = 1) and the O_ex those carry, plus its free water to the ocean. The pair
(C, O_ex) is the redox state of the departing gas in the S.1 currency, so the transfer is exact
without a solver call; the surface box re-speciates at its own conditions (S.20). O_ex keeps its
sign: a reduced gas delivers reducing power. r_deg ≡ O_ex/n_C of this step's erupted parcel, molar,
is the oxidation state of the degassed volatiles per carbon: +2 CO₂, +1 CO, 0 graphite, −2 CH₄;
H₂O is O_ex-neutral so hydrogen enters only through H₂ or CH₄. It is a flux quantity, reset every
step; Θ + 2 of S.20 is the same ratio over box A's whole inventory. r_deg sets the surface
steady state (S.20). Bound water of erupting tracers follows `outgas_fraction_water`
(0 in the reference runs) unless `eruption_zeroH2O`.

The species split reported for a degassed parcel follows the CO₂/CH₄ endmember rule of S.20:
x_CO2 = (r_deg + 2)/4. Counters: `toex_dgC`, `toex_dgO`, `toex_dgH` per step; cumulative
outgassed carbon in carbon_budget(2).

Code: `toex_degas_erupting`, `toex_degas_meltH`, `toex_flush_degas`; `outgas_traceelements` in
`tracers.f90`.

## S.20 Surface box A and B: Θ, speciation, burial

**Event.** Diagram 3.4: `Speciate A`, `Burial A to B`, `B feeds ingassing`.

<img src="md_figs/S20_surface_speciation.png" alt="S.20" width="1100">

*Figure S.20. (a) How the carbon of box A splits between CO₂ and reduced carbon along Θ, with the three regimes named. (b) Free O₂ and reduced hydrogen along Θ: one is zero wherever the other is not. (c) Burial per unit of degassed carbon at f_org = 0.2: below r = 1.6 the availability floor raises the organic share so burial takes all the oxygen and A is left neutral; above it the fixed 20 % organic share leaves free O₂ behind.*

Boxes: A = atmosphere + ocean (C_A, O_A, plus hydrogen slots), B = ingassable crust/sediment
(C_B, O_B). Currency (C, O_ex), not an oxidised/reduced binary, because CH₄ is a reductant.

The deficit coordinate, molar:

```julia
Θ ≡ O_ex/n_C − 2   (20.1)   # Box A redox coordinate: excess oxygen per carbon minus the two oxygens CO₂ needs.
```

Θ + 2 is r_deg (S.19) taken over the box A stock instead of one step's degassed parcel.

Why the threshold sits at 2 per carbon. Water is invisible to the ledger (O_ex = n_O − n_H/2,
S.1), so the hydrogen of box A costs nothing as H₂O and carbon is the only species in the box
that can hold excess oxygen, at most 2 per atom as CO₂. Rearranging (20.1) splits the ledger
into the CO₂ bill and what is left after it:

```julia
O_ex = 2 n_C + Θ n_C                                                                   (20.1a)  # The box A ledger as the cost of all carbon as CO₂ plus the remainder Θ n_C.
```

Θ > 0: the bill is paid, every carbon is CO₂, and the remaining Θ n_C oxygen atoms have no
carrier left, so they are free O₂ = Θ n_C/2 molecules. Θ = 0: exactly CO₂ + H₂O. Θ < 0: the
bill cannot be paid, so |Θ| n_C atoms are missing and some carbon, then some hydrogen, must stay
reduced. Example: n_C = 10, O_ex = 26 gives Θ = 0.6; CO₂ takes 20 atoms, 6 remain, 3 O₂. If
box A ever gains another oxygen sink (a crust or seawater iron reservoir), that term sits between
CO₂ and free O₂ and the threshold moves; there is none, so it is 2. The deficit case is one
equation and two unknowns:

```julia
Θ n_C = −(2 C_red + ½ H_red)   (20.2)   # A negative Θ is a deficit paid by reduced carbon (2 O each) and reduced hydrogen (½ O each).
```

The split is a choice. The minimise-H_red rule (water is the most stable hydrogen carrier at
surface T, carbon pays first):

```julia
C_ox = min(n_C, max(O_ex, 0)/2),  C_red = n_C − C_ox   (20.3)   # Oxidised carbon is what the oxygen can pay for; the rest is reduced.
n_O2 = max(O_ex − 2 n_C, 0)/2                          (20.4)   # Free O₂ molecules are what remains once every carbon is CO₂.
H_red = max(−2 O_ex, 0),  H_wat = n_H − H_red          (20.5)   # Reduced hydrogen appears only when O_ex is negative; the rest is water.
```

closes exactly on all branches; H_red > 0 only where O_ex < 0. Free O₂ and reduced hydrogen
cannot coexist in A, which is why combustion needs no branch: adding CH₄ lowers the free-O₂ pool
by exactly the 2 O₂ that burning it would consume. The endmember rule x = C_red/n_C = −Θ/4 is the
degassing split of S.19.

Burial A → B, two independent scalars (never "move carbon at A's O/C", which would make organic
burial meaningless):

```julia
F_bur = min( C_A (1 − e^(−dt/τ_burial)),  headroom(B) ),   C_B ≤ (w_max_C/100) M     (20.6)   # Burial flux per step: exponential drain of box A carbon, capped by the headroom of box B.
dO = 2 (1 − f_eff) dC (M_O/M_C),   f_eff = min(max(f_org, 1 − r/2), 1),  r = Θ + 2   (20.7)   # Oxygen buried with the carbon; the organic share f_eff is at least what the parcel O/C allows.
```

The residual 2 f_org dC stays in A: that asymmetry is the whole redox engine of the surface
(organic burial makes free O₂). The availability floor 1 − f_eff = C_ox/n_C says carbonate can
only be made from carbon that is CO₂. Burial moves C and O and never H. Steady state
f_org = 1 − r_deg/2: r_deg = +2 (CO₂) → 0, +1.6 → 0.2 (modern Earth), 0 → 1, r_deg < 0 → f_org > 1
impossible, so a mantle degassing net-reducing carbon cannot sustain an oxidising surface without
escape (S.21). When B saturates, burial backs up and C_A rises, a result, not a clamp.

Any O/C a par may set (`initCarbon_oex_per_C`, `ingas_oex_per_C`, `OC_B_prescribed`) is checked
into [0, 2] at init: above 2 the surplus belongs to no carbon atom, below 0 it asserts hydrogen
that is not there. `surface_mode = 'prescribed'` holds B and makes A a diagnostic.

Sources: Claire, Catling & Zahnle (2006); Zahnle, Catling & Claire (2013); Berner (2004) for
f_org; Krissansen-Totton & Catling (2017).

Code: `toex_surface_step` (order: return, burial, speciation, escape), `carbon_budget` slots
5–11 (C.2 of the implementation notes).

## S.21 Hydrogen escape

**Event.** `H escape at diffusion limit (S.21)`.

<img src="md_figs/S21_escape.png" alt="S.21" width="1100">

*Figure S.21. Left: the diffusion-limited drain time against the background pressure for three mean molecular weights. Right: reduced hydrogen remaining and the oxygen credited to A in partial mode.*

What can leave: water is cold-trapped at the tropopause, H₂ and CH₄ are not, so the exposed
reservoir is the reduced hydrogen of (20.5), slot 11, and only that. An underestimate by roughly
the cold-trapped water term, a near-fixed offset, negligible on the reduced branch.

Diffusion limit, Hunten (1973) in Walker's (1977) total-mixing-ratio form:

```julia
Φ_lim = b f_tot(H)/H_a ≈ k f_tot(H),   k = esc_phi_coef = 2.5·10¹³ cm⁻² s⁻¹   (21.1)   # Diffusion-limited hydrogen escape flux, proportional to the hydrogen mixing ratio.
dn_H/dt = −(k A/(N_A N_atm)) n_H,   τ = P N_A/(g μ k 10⁴)                     (21.2)   # The escape as a first-order decay of box A hydrogen with its time constant τ.
```

The area cancels with a hydrostatic background N_atm = P A/(g μ): the drain time depends only on
the gas the hydrogen crosses (P = esc_atm_bar, μ = esc_atm_mmw). Earth: τ ≈ 2.7·10⁴ yr. Modes:
'none', 'instant' (τ → 0), 'partial' (21.2); P → 0 recovers instant, P → ∞ recovers none. Above
f_tot ~ 10⁻² escape becomes energy-limited and (21.1) overestimates it, not capped.

The oxygen credit, per step:

```julia
carbon_budget(11) −= ΔH,   carbon_budget(6) += ΔH M_O/M_H2O   (21.3)   # Bookkeeping of one escape step: hydrogen leaves, its oxygen is credited to box A.
```

one O per two H stays behind, and slot 11 must be decremented by hand or the next speciation
pulls the escaped hydrogen back out of the ocean. Escape is the only irreversible oxidant in the
model; burial moves C and O and never H.

Sources: Hunten (1973); Walker (1977); Catling, Zahnle & McKay (2001).

Code: `toex_surface_step` escape block; `toex_esc_tau`, `toex_esc_frac`, `toex_escH`, `toex_escO`.

## S.22 The B → A return

**Event.** `Return B to A (S.22)`.

<img src="md_figs/S22_return.png" alt="S.22" width="1100">

*Figure S.22. Left: the oxidative-weathering gate g(Θ). Right: the two return legs at their default time scales.*

B inverts exactly into two endmembers, r_B = (O_B/C_B)(M_C/M_O) ∈ [0, 2]:

```julia
C_carb = C_B r_B/2,   C_org = C_B (1 − r_B/2)   (22.1)   # Box B split into carbonate and organic carbon from its O/C ratio r_B.
```

Kerogen CH₂O has O_ex = 0 like graphite, so ignoring organic hydrogen is exact. Two legs, gated
differently, which is the only reason the return can act on redox at all:

- oxidative weathering, CH₂O + O₂ → CO₂ + H₂O: C leaves at O_ex = 0 and arrives at 0; the free O₂
  of A falls by one O₂ per carbon by the coordinate itself. Gate and cap:

  ```julia
  g = max(Θ, 0)/(max(Θ, 0) + Θ_half),   ΔC ≤ free O₂ · ½ (M_C/M_O)   (22.2)
  ```

  rate C_org (1 − e^(−dt/τ_weath_org)) g; exactly 0 whenever Θ ≤ 0 so kerogen weathering never
  conjures reduced hydrogen.
- decarbonation, CaCO₃ + SiO₂ → CaSiO₃ + CO₂: C leaves at +2 and arrives at +2, free O₂ exactly
  invariant; rate C_carb (1 − e^(−dt/τ_decarb)), ungated.

Each leg draws only from its own endmember, which keeps 0 ≤ O_B ≤ 2 C_B (M_O/M_C) without a
clamp. The return runs before burial on the start-of-step B. Both legs are off by default
(`weather_organic`, `weather_decarb`). No thermostat: that would need the A → B rate to depend
on pCO₂ (Walker, Hays & Kasting 1981), not built.

Sources: Urey (1952); Berner, Lasaga & Garrels (1983); Bergman, Lenton & Watson (2004); Lasaga &
Ohmoto (2002); Kerrick & Caldeira (1998); Sleep & Zahnle (2001).

## S.23 Ingassing gates

**Event.** `Ingassing: water, box-B carbon (S.23)`, the "fractured?" and "converging trench?"
arrows, `B feeds ingassing (S.23)`.

<img src="md_figs/S23_ingassing.png" alt="S.23" width="1100">

*Figure S.23. Left: the fracture gate, yield and threshold for a dry and a wet lid. Right: the fraction of box B a fractured cell takes per step against the trench-normal convergence rate, with the floor.*

Where: cells whose top lies within `ingas_depth` (8 km) of the surface.

Fracture gate, shared by water and carbon:

```julia
fractured  ⇔  σ_II ≥ ingas_stress_ratio · min(σ_Y,ductile, σ_Y,brittle)(P, w_H2O)   (23.1)   # A lid cell counts as fractured when its stress reaches a fraction of the lower yield stress.
```

with the water-weakened envelope of `calc_viscosity` (`toex_fracture_yield`), ratio 0.5,
mid-range for dilatancy onset. Strain rate is not a fracture criterion: a hot cell creeping fast
has no crack in it.

Water: erupting lava and fractured solid within the layer hydrate to s = WaterSatCutoff × S
(or `topH2O_val` if `WaterFullSatSurfIngas` is off), and only where `oc_contact` = 1 under
`InitH2OMode = 'Selective'`. Retained melt at depth does not (it refilled every step and drained
the ocean). The ocean supplies the mass; a bone-dry ocean ingasses nothing.

Carbon: the rate must be a convergence, because the trench normal n is an eigenvector of the
surface stress Hessian with no sign. The normal strain rate along n is even in n:

```julia
conv = −n·∇u·n = −(n₁² ∂ₓuₓ + n₁n₂(∂_y uₓ + ∂ₓu_y) + n₂² ∂_y u_y)                               (23.2)   # Trench-normal convergence rate from the velocity gradient projected on the face normal.
ingfrac = max(conv, 0) dt,  as v_n dt/L_n with v_n = max(conv, 0) L_n, L_n = |(n₁ dx, n₂ dy)|   (23.3)   # Fraction of box B taken per step: the convergence Courant number across the cell.
```

the trench-normal Courant number. A floor `ingas_conv_min` = 10⁻⁵ removes the generic
convergence of a closed surface (a pre-subduction surface runs at ~10⁻⁷, a trench at 10⁻⁴–10⁻¹)
at a cost of 0.3 % of the real subducted mass. The cell takes ingfrac of box B's carbon at B's
own O/C, spread over its fractured tracers, capped per tracer at `w_max_C` = 5 wt% C; refusals
are counted and the availability test is `carbon_retrace_tol`-scaled so a one-ulp shortfall does
not drain a cell. The residue is carried signed so B loses exactly what the tracers gain.

Code: `ingas_carbon_on_tracers`, `retrieve_vn_trench`, `toex_fracture_yield`,
`sync_available_carbon_per_cell` in `volatileprocs_mod.f90`; the water block in
`timestep_traceelements_grainsize` (`tracers.f90`).

## S.24 Conservation laws and the drift test

**Event.** Every ledger print; the closed-run test.

<img src="md_figs/S24_conservation.png" alt="S.24" width="1100">

*Figure S.24. Relative change of total carbon and total hydrogen over the local run to 2.92 Gyr (core cooling on).*

- Carbon: n_C conserved absolutely; carbon_budget(4) constant to round-off, always.
- Hydrogen: conserved except escape; water_budget(6) = tracers + ocean + slot 11.
- O_ex: conserved except escape, since O_ex = n_O − n_H/2. Water is O_ex-neutral, so the ocean
  needs no O_ex slot; the only source is escape; every reduced-species conversion is a transfer,
  and a path that moves hydrogen and forgets its oxygen creates O_ex one-signed.

On a closed run with escape inert, total C, H and O_ex are flat to wandering round-off; a
monotone drift is a fault, and its sign says which: O_ex falling is a leak (a converged residual
written into the tracers), O_ex rising is oxygen manufactured downstream (a floored slot). The
invariant is stated on the cell total before and after the writeback. Verified: carbon spread
2·10⁻¹² over 261 steps, hydrogen 7·10⁻¹¹, O_ex held to ten digits for 119 steps.

Code: `ComputeVolatileInventory`, `ComputeTotalOxygenCarbonMass` (the printed intermediate
carbon line excludes fluid carbon; the stored slot does not).

## S.25 Constants and sources

Constants: Mt = 10³ g; outer bracket [−120, +40], ≤ 500 iterations; ideality passes 2–30, tol
10⁻⁸; wEMOG 0.05; toex_XCO2max 0.5; toex_XH2max 0.3; toex_fluid_Pmax 50 GPa; aC bracket
log₁₀ ∈ [−300, log₁₀ aC_use], ceiling walk ≤ 320 decades then 60 edge halvings; λ bracket log₁₀ ∈ [−60, 0]; isotonic
tilt 10⁻⁶; warm-start windows 0.25 / 2.0 log units; tolerance (3.2); softmin2(a, b, β) =
−softmax2(−a, −b, β) for every cap.

Sources by use:

- Fluid EOS and speciation: Zhang & Duan (2009) GCA 73:2089; Holland & Powell (1998) J. Met.
  Geol.; IAPWS (2008).
- Oxygen fugacity and the solid: Stixrude & Lithgow-Bertelloni (2024) sb24; Deng et al. (2020);
  Hillert & Jarl (1978); Frost (1991); Frost & McCammon (2008) Annu. Rev. 36:389.
- Melt and solid storage: Sun & Yao (2026); Katz, Spiegelman & Langmuir (2003) G³ 4:1073; Rozel
  (2018); Hirschmann et al. (2012); Hirschmann (2022); Keppler (a_H2O^n); Connolly (2018, pers.
  comm., solubility tables).
- Two-phase flow: Connolly & Podladchikov (2012); Rudge (2018).
- Melting: the list under S.16.
- Escape and surface redox: Hunten (1973) J. Atmos. Sci. 30:1481; Walker (1977); Catling, Zahnle
  & McKay (2001) Science 293:839; Claire, Catling & Zahnle (2006) Geobiology 4:239; Zahnle,
  Catling & Claire (2013) Chem. Geol. 362:26.
- Weathering and return: Urey (1952); Ebelmen (1845); Walker, Hays & Kasting (1981) JGR 86:9776;
  Berner, Lasaga & Garrels (1983) AJS 283; Berner (2004); Bergman, Lenton & Watson (2004) AJS
  304; Lasaga & Ohmoto (2002) GCA 66:361; Kerrick & Caldeira (1998) Chem. Geol. 145:213; Sleep &
  Zahnle (2001) JGR 106:1373; Krissansen-Totton & Catling (2017) Nat. Comm. 8:15423.
- Hydrous density: Jing & Karato (2011) GCA 75:6780, (2012) GCA 85:357; Lange & Carmichael
  (1987); Ochs & Lange (1999); Xu et al. (2026) JGR; Gerya, Podlesskii, Perchuk & Maresch (2004)
  PCM 31:429; Hemley et al. (1987) Nature 330:737; Karki, Ghosh & Karato (2021) Sci. Rep.
  11:10588; Mookherjee, Stixrude & Karki (2008) Nature 452:983; Bajgain, Ghosh & Karki (2015)
  Nat. Comm. 6:8578; Matsukage, Jing & Karato (2005) Nature 438:488; Sakamaki, Suzuki & Ohtani
  (2006) Nature 439:192; Gerya, Bardi, Karato & Murakami (2026) arXiv:2602.13820.

## S.26 Root finders: bisection, regula falsi, Illinois

**Event.** `Bracket step (S.12)`, `Solve aC (S.7)`, `Hydrogen split (S.10)`.

<img src="md_figs/S26_root_finders.png" alt="S.26" width="1100">

*Figure S.26. (a) A convex residual on a wide bracket, the first secant and the first bisection point. (b) Error per iteration: bisection halves the bracket every step, plain regula falsi keeps one stale endpoint and creeps in from one side, Illinois halves that endpoint's weight and recovers superlinear convergence.*

All three keep a bracket [a, b] with F(a) F(b) < 0, so none can leave the interval or diverge;
they differ in where the next trial point goes.

- **Bisection** takes the midpoint, c = (a + b)/2, and keeps the half in which the sign changes.
  Each step halves the bracket regardless of the function, so the error is 2⁻ᵏ of the initial
  width: 40 steps take [−120, 40] in f to 10⁻¹⁰, no more and no less. It never uses the values of
  F, only their signs, which is why a table-driven F with a step-function derivative cannot
  upset it.
- **Regula falsi** (the secant on the bracket) takes the point where the chord through
  (a, F(a)), (b, F(b)) crosses zero,

  ```julia
  c = (a F(b) − b F(a)) / (F(b) − F(a))   (26.1)
  ```

  and replaces the endpoint on the same side as F(c). Where F is nearly linear this reaches the
  root in a few steps. On a convex or concave F one endpoint is never replaced, the chord always
  crosses on the same side, and the iterates creep toward the root linearly from that side; the
  bracket never collapses. The residual F(f) of S.3 is exactly that shape wherever a sink has
  saturated on one side of the root.
- **Illinois** is regula falsi with one rule added: if the same endpoint is retained twice in a
  row, its function value is halved for the next chord (F(a) ← F(a)/2). The chord then tilts
  toward the stale endpoint, the next point lands on the other side, the bracket collapses from
  both ends, and the convergence order is 1.44 (superlinear). Nothing else changes: the bracket
  is still a bracket, so the safety of bisection is kept.

How TOex uses them. The outer solve on f starts with bisection until the bracket is small
enough to be locally regular, then Illinois to the tolerance (S.12); a collapsed bracket with |F|
still large is reported as failure, not as a root. The aC solve (S.7) is Illinois on a log
residual log₁₀(D/nC_total), because D ∝ aC at small activity makes that residual a straight line
in log₁₀ aC and Illinois needs about 16 evaluations where bisection over 300 decades would need
several hundred. The hydrogen split (S.10) is Illinois on a linear residual for the saturated
branch and on log₁₀ λ for the undersaturated one, and if it runs out of steps with the bracket
intact its last iterate is accepted, since the bracket is still a bracket.

Code: `toex_illinois`, `toex_partition` (outer), `toex_solve_aC`, `toex_solve_hydrogen`; Julia
`bisect_solve`, `safe_root`.
