# TOex — overview

Companion: [TOex_math.md](TOex_math.md) holds the procedure, equations and sources behind every
event named below as `(S.n)`. Read the flowcharts here, follow the section numbers there.

## 1. Introduction

TOex is the volatile-redox module of the StagYY fork in `~/Desktop/watercarb`. It carries three
element ledgers on every tracer and in two surface boxes: hydrogen (as H₂O-equivalent), carbon,
and *excess oxygen* O_ex, the oxygen a parcel holds beyond a fully reduced reference (Fe as FeO,
C as graphite, H as H₂O). One number therefore records the redox state of any reservoir: CO₂ and
carbonate carry +2 O per C, CH₄ −2, H₂ −1, ferric iron +½ per Fe, metallic iron −1.

The core is a one-unknown solver. In a cell at (P, T) with a solid, a melt and a free C-O-H fluid,
all phases share one oxygen fugacity. Instead of solving for every phase's loading, TOex solves
for the shared log₁₀ fO₂ at which the oxygen the sinks demand equals the cell's O_ex budget. Every
other quantity, Fe³⁺/ΣFe of solid and melt, dissolved CO₂ and H₂O, the fluid's speciation and
density, graphite versus carbonate, is read off at that root.

Around the solver, StagYY moves the ledgers: tracers advect them, the melting law reads the
carbonate and water they carry, a Darcy sweep moves the free fluid up columns as an element
triple, erupting tracers degas by speciation into an atmosphere-ocean box, burial and weathering
exchange with a crustal box, hydrogen escapes, and subduction ingasses surface carbon and water
back into the mantle.

The reference implementation is `src/TOex.jl` (StagPP.jl); the Fortran port lives in
`volatiles_mod.f90` (solver), `volatileprocs_mod.f90` (transport, surface, bookkeeping) and
`volatiles_bind.f90` (table hook). A stand-alone copy of the solver is
`Perple_X_MAGEMin/MAGEMin/TOex.jl` in the StagYY tree.

## 2. Goal

Model the coupled water-carbon-oxygen cycle of a convecting mantle with a redox state that is a
*consequence*, not a prescription:

- the mantle's fO₂ follows from its iron and its volatile inventory, cell by cell, every step;
- carbon speciates by that fO₂ (diamond/graphite, carbonate, CO₂ in melt, C-O-H fluid), which
  decides whether it melts, moves and degasses;
- the surface's oxidation state emerges from what degasses, what is buried and what hydrogen is
  lost to space;
- all of it conserves carbon exactly, hydrogen except for escape, and O_ex except for escape.

The reference runs that anchor every volatile experiment are the nine `EarthWater/Reference`
cases (water only, carbon only, both, at three reference viscosities).

## 3. Workflow

Rectangles are events, arrows carry the question that selects them. `(S.n)` is the section of
`TOex_math.md` with the procedure.

### 3.1 The TOex solver (one cell)

Inputs: P, T, basalt fraction p, melt fraction ϕ, and the cell's O_ex, carbon and hydrogen.
Output: fO₂ and the full partition. Fortran `toex_partition`, Julia `partition_Oₑₓ`.

```mermaid
%%{init: {"flowchart": {"wrappingWidth": 2000, "nodeSpacing": 30, "rankSpacing": 45, "htmlLabels": true}, "markdownAutoWrap": false, "themeVariables": {"fontSize": "14px"}, "themeCSS": ".edgeLabel, .edgeLabel p, .edgeLabel span, .edgeLabel div { white-space: nowrap; max-width: none; }"}}%%
flowchart TD
    classDef q fill:#4a4a4a,stroke:#222,color:#fff
    S0["Convert the tracer budgets to grams on a 1 kg cell (S.2)"]
    S1["Set the Fe3+ capacities of solid and melt and the residual tolerance (S.3, S.5)"]
    S2["Build the solid oxybarometer from the fO2 and metal tables (S.4, fig. 7)"]
    S2b["Restrict the fluid to pure water: no speciation, no EOS (S.9)"]
    S3["Bracket log10 fO2 on the warm-start window, then the full range (S.12, fig. 1)"]
    S15["Report the cell FAILED and keep its previous partition (S.12)"]
    S5["Solid Fe3+ demand: invert the oxybarometer at the trial f (S.4, S.5)"]
    S6["Melt Fe3+ demand: invert the Hirschmann oxybarometer at f (S.5)"]
    S7["Split the hydrogen between melt water, melt H2 and the fluid; no fluid if the melt demands more H than owned (S.10)"]
    S8["Carbon-activity ceiling aC_cap from the EMOG buffer at f, lowered until the speciation is valid (S.6, S.7)"]
    S9["Solve aC so that melt plus fluid carbon equals the budget (S.7)"]
    S9b["Hold aC at the ceiling; the remainder precipitates as graphite or carbonate (S.6)"]
    S10["Speciate the C-O-H fluid at this f, aC and phi (S.9)"]
    S11["Refresh the fugacity coefficients with the ZD09 EOS (S.11)"]
    S12["Evaluate the oxygen residual F(f), demands minus budget (S.3)"]
    S13["Take a bisection or Illinois step on the bracket (S.12, S.26)"]
    S14["Re-derive the state at the root and label the regime (S.13, fig. 3)"]

    S0 --> S1
    S1 --> Q1(["solid present?"]):::q
    Q1 -->|yes| S2
    Q1 -->|no| S3
    S1 -->|"P > fluid Pmax? yes"| S2b
    S2b --> S3
    S2 --> S3
    S3 --> Q2(["bracket sign change?"]):::q
    Q2 -->|no| S15
    Q2 -->|yes| S5
    S5 -->|"solid at a cap? (S.5)"| S6
    S6 --> Q3(["any C or H?"]):::q
    Q3 -->|no| S12
    Q3 -->|yes| S7
    S7 --> S8
    S8 --> S9b
    S9b --> Q4(["ceiling demand > budget?"]):::q
    Q4 -->|yes| S9
    Q4 -->|"no: C-saturated"| S10
    S9 --> S10
    S10 --> Q5(["any non-finite fraction?"]):::q
    Q5 -->|"yes: lower the ceiling"| S9
    Q5 -->|no| S11
    S11 --> Q6(["phi, C, H settled?"]):::q
    Q6 -->|no| S7
    Q6 -->|yes| S12
    S12 --> Q7(["residual in tolerance?"]):::q
    Q7 -->|yes| S14
    Q7 -->|no| S13
    S13 --> Q8(["bracket collapsed?"]):::q
    Q8 -->|yes| S15
    Q8 -->|"no: next trial f"| S5
```

At the root the model reports sOₑₓ, mOₑₓ (iron of solid and melt), carbonO (melt CO₂), carbO
(solid carbonate), fₑₓ (net fluid oxygen, signed), mH2ₑₓ (melt H₂), the melt's water, the fluid's
mole fractions, molar mass and density, and a regime label. Carbon and hydrogen close exactly at
every trial f; the outer solve is on oxygen alone.

### 3.2 The StagYY implementation (one timestep)

The solve sits in the melting block of the tracer step (`TraMeltErupt`). Before it, the cell
budgets are produced; after it, what the solve placed is moved.

```mermaid
%%{init: {"flowchart": {"wrappingWidth": 2000, "nodeSpacing": 30, "rankSpacing": 45, "htmlLabels": true}, "markdownAutoWrap": false, "themeVariables": {"fontSize": "14px"}, "themeCSS": ".edgeLabel, .edgeLabel p, .edgeLabel span, .edgeLabel div { white-space: nowrap; max-width: none; }"}}%%
flowchart TD
    classDef q fill:#4a4a4a,stroke:#222,color:#fff
    T1["Stokes solve: velocities, stress and strain rate on the grid"]
    T2["Ingassing: hydrate the lid to the cutoff; box-B carbon enters only at a converging trench (S.23)"]
    T3["Advect the tracers with their H, C and O_ex slots"]
    T4["Melt and erupt on the volatile solidus: Katz water, carbonatite plateau where oxidised C exceeds C_min (S.16, fig. 4)"]
    T5["Update solid water storage S_eff = S times a_H2O to the n, unit activity without a fluid (S.15)"]
    T6["Gather the cell: masses, melt fraction, ten volatile slots; previous fO2 narrows the bracket (S.14)"]
    T7["Call TOex for the cell, section 3.1 (fig. 1)"]
    T8["Scatter the solution: incremental to rock, proportional to fluid; a failed cell keeps its old partition (S.14)"]
    T9["Speciated degassing of erupting tracers into box A and the ocean (S.19)"]
    T10["Darcy sweep of the fluid triple up each column, parcels reaching the top go to the surface, diagram 3.3 (S.18)"]
    T11["Surface step: return, burial, speciation, escape, diagram 3.4 (S.20 to S.22, fig. 2)"]
    T12["Tracers to grid: composition fields, fluid density and viscosity, hydrous density (S.17, figs. 5, 6)"]

    T1 --> QT1(["lid cell fractured?"]):::q
    QT1 -->|yes| T2
    QT1 -->|no| T3
    T2 --> T3
    T3 --> T4
    T4 --> T5
    T5 --> T6
    T6 --> T7
    T7 --> T8
    T8 --> QT2(["tracer erupting?"]):::q
    QT2 -->|yes| T9
    QT2 -->|no| T10
    T9 --> T10
    T10 --> T11
    T11 --> T12
    T12 -->|"next step"| T1
```

Two lags are deliberate and one step long: the melting law reads last step's carbonate fraction,
the solid storage reads last step's water activity. Fluid density and viscosity come from the
solve of the same cell, so they need no grid-to-tracer interpolation.

### 3.3 The Darcy sweep (free fluid)

The fluid is the only reservoir that moves other than by advection. It moves as the element
triple (H as H₂O-equivalent, C, O_ex), never as a redox label. Fortran
`WaterSaturationAndVerticalRedistribution`.

```mermaid
%%{init: {"flowchart": {"wrappingWidth": 2000, "nodeSpacing": 30, "rankSpacing": 45, "htmlLabels": true}, "markdownAutoWrap": false, "themeVariables": {"fontSize": "14px"}, "themeCSS": ".edgeLabel, .edgeLabel p, .edgeLabel span, .edgeLabel div { white-space: nowrap; max-width: none; }"}}%%
flowchart TD
    classDef q fill:#4a4a4a,stroke:#222,color:#fff
    D1["Exsolve: a tracer above the cutoff sheds its free triple as a parcel, then bound water as pure H2O (S.18)"]
    D2["Face check: Darcy capacity k A g drho dt over eta of the face above (S.18)"]
    D6["Blocked part stays as free fluid in this cell, congruent split of the triple (S.18)"]
    D3["Cross the face; dissolve hydrogen only into undersaturated rock, up to S_eff minus w (S.18)"]
    D4["React in transit, if on: re-solve the crossed cell with the parcel included (S.14, fig. 1)"]
    D5["Strand: C and O_ex stay in this cell's rock, the parcel is gone (S.18)"]
    D7["Exit at the top: H to the ocean, C and O_ex to box A (S.20)"]

    D1 --> D2
    D2 --> QD1(["face capacity short?"]):::q
    QD1 -->|"yes: blocked part"| D6
    QD1 -->|"no: passing part"| D3
    D3 --> D4
    D4 --> QD2(["parcel H exhausted?"]):::q
    QD2 -->|yes| D5
    QD2 -->|no| QD3(["above the top?"]):::q
    QD3 -->|yes| D7
    QD3 -->|"no: next face"| D2
```

The sweep runs D1 for every cell of a column, bottom to top; a parcel keeps looping D2 to D4 until
it is blocked, stranded or exits. Every merge and subdomain hand-over asserts the window
−n_H/2 ≤ n_Oex ≤ 2 n_C.

### 3.4 The surface reservoirs

Two boxes, each a (C, O_ex) pair: A is atmosphere plus ocean and receives degassing; B is the
ingassable crust and sediment. Hydrogen in A splits into water and reduced hydrogen by the same
O_ex arithmetic. Fortran `toex_surface_step`, once per step after degassing.

```mermaid
%%{init: {"flowchart": {"wrappingWidth": 2000, "nodeSpacing": 30, "rankSpacing": 45, "htmlLabels": true}, "markdownAutoWrap": false, "themeVariables": {"fontSize": "14px"}, "themeCSS": ".edgeLabel, .edgeLabel p, .edgeLabel span, .edgeLabel div { white-space: nowrap; max-width: none; }"}}%%
flowchart TD
    classDef q fill:#4a4a4a,stroke:#222,color:#fff
    A1["Degassed carbon, O_ex and hydrogen arrive in box A (S.19)"]
    A2["Return B to A: decarbonation of carbonate always, oxidative weathering of organics only with free O2 in A (S.22)"]
    A3["Burial A to B as two scalars, carbon flux and oxygen at 2 (1 - f_eff) per carbon, while B is below capacity (S.20)"]
    A4["Speciate box A: Theta, oxidised carbon, free O2, reduced hydrogen (S.20, fig. 2)"]
    A5["Hydrogen escape from the reduced pool at the diffusion limit; one O per two H stays (S.21)"]
    A6["Box B supplies ingassing at its own O to C ratio (S.23)"]

    A1 --> A2
    A2 --> A3
    A3 --> A4
    A4 --> QA1(["escape on and reduced H present?"]):::q
    QA1 -->|yes| A5
    QA1 -->|"no: nothing leaves"| A6
    A5 --> A6
```

## 4. Figures

Stored in `md_figs/`, generated by `md_figs/scripts/` (see its README). Diagram events cite them as "fig. n".

<a id="fig-1"></a>
**Figure 1. TOex convergence statistics.** Per-step counters of the local run (return legs on, corrected tables, core cooling on; 2.92 Gyr so far): cells solved and failed by the main solve (T7), cells re-solved in transit by the Darcy sweep (D4), root iterations, which bracket the warm start used (S.12), ZD09 root failures and stranded fluid (S.9, S.18), and the solver wall time.

![Figure 1. TOex convergence statistics](md_figs/S12_toex_convergence.png)

<a id="fig-2"></a>
**Figure 2. Surface reservoirs of the local run (return legs on, core cooling on; 2.92 Gyr so far).** Ocean, mantle water and escape (S.21, S.24); carbon in box A, box B and the mantle; the box A ledgers and the free O₂ implied by Θ_A (S.20); exchange fluxes; redox of the boxes; where the volatiles sit (diagram 3.4).

![Figure 2. Surface reservoirs of the fully coupled run](md_figs/S20_surface_reservoirs.png)

<a id="fig-3"></a>
**Figure 3. Radial redox and carbon structure.** Horizontal means of the local run (core cooling on) at 0.73, 1.45 and 2.91 Gyr: fO₂ (S.4, S.13), Fe³⁺/ΣFe (S.4, S.5), the O_ex ledger (S.1), oxidised and reduced solid carbon (S.6, S.7), melt CO₂ (S.8).

![Figure 3. Radial redox and carbon structure](md_figs/S04_radial_redox_carbon.png)

<a id="fig-4"></a>
**Figure 4. Volatile-enhanced melting laws.** The law behind T4 (S.16) with the run parameters: the Katz water depression, the carbonated solidi with ledge and taper against the data, the carbonatite plateau in T(f), and the two gates S_c and s_l.

![Figure 4. Volatile-enhanced melting laws](md_figs/S16_volatile_melting_laws.png)

<a id="fig-5"></a>
**Figure 5. Rheological hooks.** What water does to the solid: the Karato viscosity law, the brittle and ductile yield factors (T1, S.23 fracture test), and the hydrous density tables of solid and melt (T12, S.17).

![Figure 5. Rheological hooks](md_figs/S17_rheology_laws.png)

<a id="fig-6"></a>
**Figure 6. Rheological consequences in the run.** The 2-D frame of the local run (core cooling on) at 2.92 Gyr: viscosity, bulk water, saturation and density residual maps, viscosity against saturation by depth, and radial means of water, saturation and melt.

![Figure 6. Rheological consequences in the run](md_figs/S17_rheology_in_run.png)

<a id="fig-7"></a>
**Figure 7. Where the solid oxybarometer table samples MAGEMin.** The 36 nodes of the signed ferric-ratio axis of StagfO2.dat (S.4) at 5 GPa, 1700 K: on the linear R axis, on the asinh interpolation coordinate, and their spacing. Nodes are dense near the FeO reference where log fO₂ changes by ten units over 10⁻³ in R, sparse along the log-linear tail, and tighten again towards R = 0.95 where the R/(1 − R) form diverges.

![Figure 7. Where the solid oxybarometer table samples MAGEMin](md_figs/S04_R_sampling.png)

<a id="fig-8"></a>
**Figure 8. Carbonate leg of the melting law: fixed shift versus lever rule (scratch build, not ported).** (a) T_melting(f) at 15 GPa for three carbonate contents under the coded law (silicate solidus lowered by a fixed 100 K, solid) and the lever rule of eq. 16.12 (dashed). (b) The same as f(T): the shift gives one silicate productivity whatever c_oxC, the lever rule scales it with c_oxC. (c) Both laws on the horizontal-mean profile of the local run at 2.91 Gyr against the melt fraction it carries. (d) Two 20-step continuations of the 2.5 Gyr production frame with the scratch binary, `carb_lever_rule` off and on. The lever rule with a linear liquidus produces more transition-zone melt than the shift, because the eutectic drop there is ~600 K and even diluted five-fold it exceeds 100 K; see S.16.

![Figure 8. Carbonate leg of the melting law: fixed shift versus lever rule](md_figs/S16_lever_vs_shift.png)

Missing and worth making: a (P, T) map of where the C-O-H fluid exists and where the solve fails
(S.9, S.12); a schematic of the Darcy parcel with its two transfer rules (S.18).

## 5. Data flow through the code

The cell solve of section 3.1 read as code. Each box names a routine and the helpers it calls
(→), or a data structure; each arrow carries the variables handed to the next box. Three levels,
each zooming into one box of the level above.

### 5.1 One cell: tracers in, tracers out

```mermaid
%%{init: {"flowchart": {"wrappingWidth": 2000, "nodeSpacing": 30, "rankSpacing": 45, "htmlLabels": true}, "markdownAutoWrap": false, "themeVariables": {"fontSize": "13px"}, "themeCSS": ".edgeLabel, .edgeLabel p, .edgeLabel span, .edgeLabel div { white-space: nowrap; max-width: none; }"}}%%
flowchart TD
    A["tra(slot, ii, ib): OEXPOS, OxCPOS, ReCPOS, H2OPOS, H2OFREEPOS, CFLUIDPOS, OFLUIDPOS, TEMPPOS, fO2POS, MASSPOS"]
    B["toex_cell_gather"]
    C["toex_partition"]
    D["toex_build_ctx → toex_endmember_composition, toex_mass_budget, toex_solve_IDV, toex_sOwt_coeffs, toex_fill_sfO2"]
    E["toex_bisect_solve → toex_residual → toex_derive (section 5.2)"]
    F["toex_derive, classify_regime"]
    G["res (toex_result_t)"]
    H["toex_shares_to_masses"]
    I["toex_cell_scatter → toex_apply_oex, toex_apply_slot"]
    J["tra(slot, ii, ib)"]

    A -->|"tra(slot)·m·1e-2 per tracer, TEMPPOS·m, fO2POS·m, basalt·m"| B
    B -->|"cl: Pgpa, Tk, pbas, phi, TOex_in, TC_in, H_in, fO2 (warm start)<br/>w0(1:10) kg per slot"| C
    C -->|"P, T, p_mix, phi, TOex, TC, H"| D
    D -->|"ctx: TOex, Mf, Ms, Mc, Mh, nOxideTotal, Phi, maxsOex, minsOex, maxmOex<br/>clipS, clipM, IDV, sfO2, xfrac, nFe_solid, nFe_melt, solid_absent<br/>fluid_ok; maxXCO2"| E
    E -->|"converged, x = log10 fO2, resnorm, iters"| F
    F -->|"st: sOex, mOex, carbonO, carbO, fex, mH2ex, XCO2, cs; regime"| G
    G -->|"res: fO2, aC, sOex, mOex, carbonO, carbO, cohsink, meltH2ex, mC_gph<br/>mC_carb, mC_melt, mC_fluid, mH_fluid, mH_melt, XH2O..XH2, vol_fluid<br/>FeR, FeM, iters, resnorm"| H
    H -->|"w(1:10) kg: S_OEX = sOex + carbO, M_OEX = mOex + carbonO + meltH2ex<br/>F_O = cohsink, S_REC = mC_gph, S_OXC = mC_carb, M_OXC = mC_melt<br/>F_C = mC_fluid, F_H = mH_fluid, M_H = rest"| I
    I -->|"tra(slot) += (w − w0)·100/m per tracer; fO2POS = x"| J
```

`res` also fills the grid diagnostics (`toex_diag`: dFMQ, dIW, melt XCO2 and XH2, Fe³⁺, Fe⁰, COH index) and the per-cell fluid density `toex_rhof` from `vol_fluid`.

### 5.2 Inside toex_derive: the six sinks at one trial f

```mermaid
%%{init: {"flowchart": {"wrappingWidth": 2000, "nodeSpacing": 30, "rankSpacing": 45, "htmlLabels": true}, "markdownAutoWrap": false, "themeVariables": {"fontSize": "13px"}, "themeCSS": ".edgeLabel, .edgeLabel p, .edgeLabel span, .edgeLabel div { white-space: nowrap; max-width: none; }"}}%%
flowchart TD
    A["toex_derive(f, ctx)"]
    B["invert_demand(1) → demand_g → toex_solid_fO2 → toex_R_of_sOwt, toex_solid_fO2_proxy"]
    C["invert_demand(2) → demand_g → toex_melt_fO2"]
    D["toex_coh_state (section 5.3)"]
    E["st (toex_state_t)"]

    A -->|"f, ctx"| B
    B -->|"st%sOex (signed, metal branch < 0), st%fs; f, ctx"| C
    C -->|"st%mOex, st%fl; f, ctx"| D
    D -->|"cs → st%carbonO = Phi·XCO2_melt/(1 − XCO2_melt), st%carbO, st%fex<br/>st%mH2ex, st%XCO2; st%total = sum of the six, st%ok"| E
```

### 5.3 Inside toex_coh_state: carbon, hydrogen and the fugacity loop

```mermaid
%%{init: {"flowchart": {"wrappingWidth": 2000, "nodeSpacing": 30, "rankSpacing": 45, "htmlLabels": true}, "markdownAutoWrap": false, "themeVariables": {"fontSize": "13px"}, "themeCSS": ".edgeLabel, .edgeLabel p, .edgeLabel span, .edgeLabel div { white-space: nowrap; max-width: none; }"}}%%
flowchart TD
    A["toex_coh_state inputs: f, T, P, Mh, Mc, nOxideTotal, Mf, molfracXB"]
    B["toex_reaction_Ks, toex_invert_carbon_buffer_fO2, toex_logfO2_EMOG → hc (toex_hc_t)"]
    C["toex_solve_hydrogen → toex_illinois(2 or 3) → toex_inner_res → toex_melt_H → toex_h2_melt_moles, toex_h2o_melt_molH"]
    D["toex_solve_aC → toex_illinois(1) → toex_inner_res → toex_carbon_demand → toex_coh_reaction_step"]
    E["toex_coh_reaction_step, toex_melt_H"]
    F["toex_zd_fugacity → zd_findroot"]
    G["ideality exit test (S.11); not settled: back to toex_solve_hydrogen"]
    H["toex_coh_state, post-loop: toex_softcap, toex_zd_fugacity (pure water above Pmax), toex_h2_melt_K, toex_melt_h2o_wt"]
    I["cs (toex_coh_t)"]

    A -->|"T, Pbar, f, P"| B
    B -->|"hc, phi, molH_fl = molH_tot"| C
    C -->|"molH_fl, lam; hc, phi"| D
    D -->|"cs%aC, aC_use; fO2, Pbar, K1..K4, phi"| E
    E -->|"vx(1:5) = XCH4, XH2O, XCO2, XH2, XCO normalised; nH2m, nHw"| F
    F -->|"vol, fug(1:5) → phi(i) = fug(i)/(vx(i) Pbar), dphi"| G
    G -->|"settled or nit_max reached: final X's, aC, molH_fl, lam, nH2m, nHw<br/>vol"| H
    H -->|"aC, XH2O..XH2, n_fluid, nC_fluid, XCO2_melt, nC_melt, nC_carb, nC_gph<br/>carbO, fex, meltH2ex, molH_fl, lam, vol_fluid, it_ideality, valid"| I
```

Above `toex_fluid_Pmax` the loop is skipped: the fluid is frozen as pure water, `toex_solve_hydrogen` and `toex_melt_H` run once with `hc%aC_fixed`, and `toex_solve_aC` is called with `fl_ok = .false.` so only the melt demand counts.

Julia mirror in `src/TOex.jl` (StagPP): `partition_Oₑₓ` = `toex_partition`, `mass_budget`, `endmember_composition`, `solve_∫ΔVdP` = `toex_solve_IDV`, `bisect_solve` = `toex_bisect_solve`, `unified_residual` = `toex_residual`, `_coh_derive` = `toex_derive`, `coh_state` = `toex_coh_state`, `_coh_reaction_step` = `toex_coh_reaction_step`, `_carbon_demand` = `toex_carbon_demand`, `zd_fugacity` = `toex_zd_fugacity`, `melt_fO2` = `toex_melt_fO2`, `solid_fO2_interpolant` = `toex_fill_sfO2` plus `toex_solid_fO2`.
