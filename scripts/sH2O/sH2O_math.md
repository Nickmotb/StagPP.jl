# sH2O — the StagYY lookup-table generator: mathematics, sources and choices

This document describes `sfO2.jl` (the caller) and `sfO2_fns.jl` (the kernels) in this folder. They write the six tables
that StagYY reads for water and redox:
- `StagH2O.dat` (solid water storage capacity);
- `StagfO2.dat` and `StagMetal.dat` (the solid oxybarometer and its native iron);
- `StagfO2step.dat` (the silica-saturation step in basalt fraction);
- `StagSdens.dat` and `StagMdens.dat` (hydrous density corrections).

`sH2O_plot.jl` draws every figure in `md_figs/` from a finished run, except the TikZ scheme S.20.

Notation:
- P in GPa and T in K, unless °C is stated.
- Capacity C in wt% of the water-bearing solid, measured at water saturation.
- a = a_H2O, the water activity of the coexisting fluid.
- μ is a molar chemical potential in J/mol.
- R = Fe³⁺/ΣFe, signed: R < 0 is the metal branch.
- x_φ is the mass fraction of phase φ in the solid, with fluid excluded.
- 1 ppm = 10⁻⁴ wt%; composition vectors are in oxide wt% unless given as mole fractions.
- Every equation block is followed by a *Symbols* list defining what it introduces.
- Code names are in `sfO2_fns.jl` unless stated otherwise.

Contents

- S.1 Tables, sectors, compositions and flags
- S.2 Water: fugacity, volume and chemical potential
- S.3 Upper mantle: ume and mbe with excess water
- S.4 Nominally anhydrous minerals: the saturation laws
- S.5 Bridgmanite and the Al_pv flag
- S.6 The MgO-SiO2-H2O model of Komabayashi & Omori (2006)
- S.7 The water-activity map
- S.8 Phase H
- S.9 The DHMS term: the question and the minimisation
- S.10 From the minimisation back to the node and the maps
- S.11 The assembled capacity maps
- S.12 Capacity along geotherms
- S.13 Seams and hand-overs
- S.14 The oxybarometer and native iron tables
- S.15 The silica-saturation step table
- S.16 Hydrous density of solid and melt
- S.17 Numerics and cost
- S.18 Choices and open items
- S.19 Constants and sources

---

## S.1 Tables, sectors, compositions and flags

<img src="md_figs/S01_domain.png" alt="S.1" width="1100">

*Figure S.1. The method that fills each part of every block's P-T domain. The DHMS domain is the part where the MSH re-equilibration of S.9-S.10 is applied. The illustrative geotherms are those of S.12.*

**Sectors.** Each table covers three P-T sectors on an nP × nT grid (production 100 × 100):

| sector | P (GPa) | T (K) |
|---|---|---|
| upper mantle (UM) | 1e-4 to 7 | 500 to 2000 |
| transition zone (TZ) | 7 to 25 | 700 to 4000 |
| lower mantle (LM) | 25 to 135 | 700 to 4000 |

**Compositions.** Every table except the density tables holds the same ten compositions per sector: basalt fraction p = 0, 1, 0.05, 0.1, 0.2, 0.5, 0.65, 0.75, 0.85 and 0.95 (`fO2_Xmix` lists the eight mixes). A mix is p X_B + (1 − p) X_H in wt% of oxides. StagYY reads between neighbouring compositions. Until 2026-09-27, `StagH2O.dat` held only p = 0, 0.3 and 1; S.11 shows why ten are needed. Bulk compositions (wt%) are from Stixrude & Lithgow-Bertelloni (2024):

```julia
X_H = [SiO2 43.43, MgO 45.93, FeO 8.34, CaO 0.90, Al2O3 1.00, Na2O 0.01, Cr2O3 0.30]   (1.1)   # Harzburgite bulk of the sb24 model, wt% of oxides.
X_B = [SiO2 50.42, MgO 9.77, FeO 7.10, CaO 12.54, Al2O3 16.80, Na2O 2.23, Cr2O3 0.07]  (1.2)   # Basalt (MORB) bulk of the sb24 model, wt% of oxides.
```

Iron is handed to MAGEMin as Fe + O at a ferric ratio R. The conversion is `oxidize_bulk`, with n_O = n_Fe(1 + R/2) in the `Fe_O` format. The storage-capacity bulks use R = 0.01 (`Rsol`).

*Symbols.*
- X_H, X_B: oxide wt% vectors of the two end-member bulks.
- p: basalt mass fraction of a mixed bulk, X = p X_B + (1 − p) X_H.
- n_Fe, n_O: moles of Fe and of O handed to MAGEMin with the iron (the Fe + O representation of FeO + Fe2O3).
- R = Fe³⁺/ΣFe: the ferric ratio; negative values put the missing oxygen into native iron.

**Capacity.** The storage capacity of a node is the water the solid holds when it coexists with free water, per unit mass of the solid:

```julia
C = Σ_φ x_φ c_φ,   x_φ = w_φ / Σ_{φ ∉ fluid} w_φ                                   (1.3)   # Capacity: phase-weighted water content of the solid assemblage, fluid excluded.
```

*Symbols.*
- C: capacity (wt% of the water-bearing solid).
- φ: a solid phase of the MAGEMin assemblage.
- w_φ: its weight fraction in the whole system (`ph_frac_wt`, fluid included).
- x_φ: its weight fraction in the solid.
- c_φ: its water content in wt%.

Here c_φ is the water content of phase φ:
- its own stoichiometric water, for a hydrous phase;
- its saturation law (S.4-S.5), for a nominally anhydrous mineral (NAM);
- zero, for the dry list.

MAGEMin's hydrous phases are recognised by H2O > 10⁻⁶ in their composition. A phase that falls in none of these three cases stops the run, so no phase is silently dropped.

**Why the tables have no water axis.** C is defined at water saturation, a = 1 (or the map's a for the DHMS, S.7), with free water in excess. Once free water is present, adding more changes only the amount of fluid, not the solid: the fluid is pure H2O, so it removes nothing from the rock. C is therefore a property of (P, T, p) alone. A rock with bulk water w holds min(w, C) in its solids, and the excess w − C is free fluid. A table in w would repeat C along that axis. StagYY's reader has none since 2026-09-27; a file with the former water-axis header (`nP nT 1`) is rejected at load. Below saturation, the NAM part of the solid's water at a < 1 is handled in StagYY (TOex_math S.15).

**How StagYY uses C in a time step.** Each tracer carries bound water w and free water w_f, and reads C from `StagH2O.dat` at its P, T and p (`eval_fO2_traWsol`). It does so twice per step: after advection, and again after melting and eruption. Figure S.20 shows the order:
1. tracers advect;
2. C is read at the new position, and a tracer whose C fell below w becomes oversaturated;
3. within δ_hyd of the surface, where the rock is fractured, undersaturated tracers take up water from the ocean towards αC_s (scaled by f_oc when the ocean runs short), and the free water of saturated tracers there goes to the ocean;
4. wet rock melts at the lowered solidus; the melt erupts or intrudes, erupted tracers equilibrate with the ocean and crystallise as new crust, and C is read again;
5. the Darcy sweep carries the free water up each column, each tracer sheds at most M_max/n_c, and undersaturated tracers absorb it;
6. tracer values are averaged onto the grid, and set viscosity, yield stress and friction for the next Stokes solve.

<img src="md_figs/S20_routine.png" alt="S.20" width="1100">

*Figure S.20. The water routine of one StagYY time step, on a 3 × 3 patch of cells below the ocean. The figure writes the capacity as S and the composition as X; here they are C and p. M_max = min(k A Δρ g ρ_w Δt / η_w, M_cell) is the Darcy limit on the water that crosses a cell face in one step, and n_c the number of tracers in the cell. Drawn in TikZ (`Drafts/Water/figures/methods/routine_scheme.tex`), not by `sH2O_plot.jl`.*

**Hosts and layout.** The phase contributing the largest x_φ c_φ is stored as the node's host (maps only; it is not written to the table). `StagH2O.dat` layout:
- header `nP nT nEM` and the eight mix fractions;
- one `P0 P1 T0 T1` line per sector;
- then per sector the basalt, harzburgite and mix blocks in header order, as in `StagfO2.dat`. Each block is nP rows × nT columns, 4 decimals, floored at 10⁻¹⁰.

**Flags (`sfO2.jl` settings).**

| flag | production | meaning |
|---|---|---|
| `Al_pv` | true | bridgmanite water from its own Al2O3 (S.5), also the post-perovskite reference (4.4c); false: 30 ppm |
| `dhms` | true | DHMS term (S.6-S.10); false: NAM laws on the sb24 assemblage only |
| `PDmax`, `PHmax` | 50, 60 GPa | no phase D beyond PDmax, no DHMS beyond PHmax (S.8) |
| `TDmax` | 1873 K | no DHMS beyond (K&O calibration limit) |
| `pmbe` | 0.65 | basalt fraction from which a composition runs on mbe instead of ume (S.3) |
| `Phand` | 5-6 GPa | UM hand-over from ume to sb24 + DHMS (S.3) |
| `Pmbe` | 9.5-10 GPa | hand-over from mbe to sb24 (S.3) |
| `melt_H2O_model` | `:ideal` | melt water volume of StagMdens (S.16) |
| `SFO2_DO` | all | environment list of passes: sH2O_um, sH2O_deep, fO2, step, dens |

Code: settings block and `Rhybridspace`, `oxidize_bulk` in `sfO2.jl` / `sfO2_fns.jl`; `assemble!` for (1.3).

## S.2 Water: fugacity, volume and chemical potential

<img src="md_figs/S02_water.png" alt="S.2" width="1100">

*Figure S.2. (a) Pitzer & Sterner fugacity. (b) The difference between the DHMS water potential (Robie + PS to 1 GPa + Brodholt & Wood above) and the potential from Pitzer & Sterner alone. (c) The Brodholt & Wood molar volume.*

Three water models are used, each where its consumer was calibrated:

| use | model | why |
|---|---|---|
| NAM saturation laws (S.4) | Pitzer & Sterner (1994) fugacity | the most widely used; the in-house fits of S.4 were made with it |
| DHMS minimisation (S.6-S.9) | Robie et al. (1978) 1 bar + PS to 1 GPa + Brodholt & Wood (1993) above | the equation of state Komabayashi & Omori (2006) calibrated against |
| density tables (S.16) | Gerya et al. (2004) | the reference generator's choice |

**Pitzer & Sterner.** The residual Helmholtz energy with the ten coefficient rows c_i(T) = Σ_k C_ik T^{e_k}, where the powers are e = (−4, −2, −1, 0, 1, 2), gives

```julia
P(ρ) = RT [ρ + c₁ρ² − ρ²(c₃ + 2c₄ρ + 3c₅ρ² + 4c₆ρ³)/(c₂ + c₃ρ + c₄ρ² + c₅ρ³ + c₆ρ⁴)² + c₇ρ²e^{−c₈ρ} + c₉ρ²e^{−c₁₀ρ}]   (2.1)   # PS equation of state in molar density ρ (mol/cm³).
ln f = ln ρ + c₁ρ + 1/(c₂ + … + c₆ρ⁴) − 1/c₂ − (c₇/c₈)(e^{−c₈ρ} − 1) − (c₉/c₁₀)(e^{−c₁₀ρ} − 1) + P/(ρRT) + ln RT − 1   (2.2)   # Fugacity from the same residual function.
```

*Symbols.*
- P: pressure. Pa inside (2.1)-(2.2); the argument of `fH2O` is in GPa.
- T: temperature (K).
- ρ: molar density of water (mol/cm³).
- R: gas constant, so RT is in cm³·Pa/mol (code `1e6 Rgas T`).
- c₁ … c₁₀: the ten temperature-dependent coefficients, c_i(T) = Σ_k C_ik T^{e_k}, with C_ik the 10 × 6 table `PS_C`.
- f: fugacity; (2.2) gives it in Pa, and `fH2O` divides by 10⁹ to return GPa.

(2.1) is solved for ρ by Newton's method with a bisection fallback (ForwardDiff derivative, bracket 10⁻³ to 10³ mol/cm³), with P capped at 100 GPa. `fH2O` returns f in GPa.

**Brodholt & Wood.** A modified Redlich-Kwong form whose co-volume is linear in V, b = b₀ + b₁V, so V − b = (1 − b₁)V − b₀:

```julia
P(V, T) = RT/((1 − b₁)V − b₀) − a(T)/(√T V ((1 + b₁)V + b₀)) + c/V + d/V² + e/V³ + f/V⁴        (2.3)   # BW H2O equation of state, P in bar, V in cm³/mol.
a(T) = −582468 − 3038.79 T − 9.24574·10⁻³ T² + 3.02674·10⁹/T²                                   (2.4)   # BW attraction term.
```

*Symbols.*
- P: pressure (bar).
- V: molar volume of water (cm³/mol).
- b = b₀ + b₁V: the co-volume.
- a(T): the attraction term (bar cm⁶ K^{1/2}/mol²).
- c, d, e, f: the BW virial constants (not the PS coefficients, and not the fugacity).

The constants are b₀ = −3.90463·10⁻², b₁ = −0.991078, c = 3.64905·10⁴, d = −1.02451·10⁷, e = −1.79681·10⁸, f = 2.18437·10⁹, R = 83.1446 cm³ bar/mol/K. V(P) is found by safeguarded Newton on [4.5, 80] cm³/mol (`bw_V`).

Because b is linear in V, ∫P dV has a closed form F(V) (`bw_F`), so the volume integral needs no quadrature:

```julia
∫_{P₁}^{P₂} V dP = P₂V₂ − P₁V₁ − [F(V₂) − F(V₁)]                                              (2.5)   # Pressure integral of V by parts; F is the antiderivative of P(V).
F(V) = RT ln((1 − b₁)V − b₀)/(1 − b₁) − a/(√T b₀) ln(V/((1 + b₁)V + b₀)) + c ln V − d/V − e/(2V²) − f/(3V³)   (2.6)   # Closed-form antiderivative used by mu_h2o.
μ_H2O(P, T) = G°_Robie(1 bar, T) + RT ln f_PS(min(P, 1 GPa), T) + ∫_{1 GPa}^{P} V_BW dP        (2.7)   # Water potential of the DHMS minimisation at a = 1.
```

*Symbols.*
- P₁ = 1 GPa (10⁴ bar) and P₂ = P: the integration limits, with V₁ = V_BW(P₁, T) and V₂ = V_BW(P₂, T).
- F(V): the antiderivative ∫P dV of (2.3) at fixed T.
- G°_Robie(1 bar, T): the Gibbs energy of water at 1 bar (Robie et al. 1978 heat-capacity row, `KO_W`), in J/mol.
- f_PS: the PS fugacity in bar.
- V_BW: the BW volume.

The factor 0.1 converts cm³·bar to J.

(2.5) matches numerical quadrature to 10⁻⁵ J. The switch at 1 GPa uses PS where BW is least accurate, at low density. Panel (b) of the figure is why each consumer keeps its own model: at 60 GPa the two potentials differ by ~50 kJ/mol, several times the free-energy margins that place the DHMS reactions.

Code: `fH2O`, `BW`, `bw_a`, `bw_P`, `bw_F`, `bw_V`, `mu_h2o`, `KO_W` (Robie 1-bar row).

## S.3 Upper mantle: ume and mbe with excess water

<img src="md_figs/S03_upper_mantle.png" alt="S.3" width="1100">

*Figure S.3. Upper-mantle capacity with host phases (top) and the share held by hydrous phases (bottom), for harzburgite, the 20 % basalt mix and basalt. Dotted lines: the 5-6 GPa hand-over of the ume compositions.*

**Databases.** The UM compositions are minimised in MAGEMin (Riel et al. 2022) with databases that carry the crustal and serpentinite hydrous phases, which sb24 lacks:
- `ume` for p < `pmbe`. This is the ultramafic set of Evans & Frost (2021), extended with amphibole, plagioclase and augite from Green et al. (2016).
- `mbe` for p ≥ `pmbe`. This is the metabasite set of Green et al. (2016), extended.

`mbe` has no Cr2O3, so that oxide is dropped there. In production `pmbe` = 0.65: the mixes from p = 0.65 and basalt run on mbe.

**Mixed compositions (test, 40 × 40 UM nodes, p = 0-1 in steps of 0.05, both databases at every p).** Neither database carries every hydrous phase:
- ume has no lawsonite at any p; mbe has no antigorite at any p.
- mbe forms lawsonite only from p ≈ 0.85-0.9 (25 % of the nodes at p = 0.9, 29 % at basalt).
- The two agree best at p = 0.6-0.7: median |log10 ratio| 0.02-0.03, more than ×2 apart at 7-9 % of the nodes. Elsewhere they differ by ×2 or more at 22-33 % of the nodes.

`pmbe` = 0.65 puts the database switch where the two agree, and lawsonite grows in over p = 0.85-1 as mbe forms it. Under the former rule (`pmbe` = 1) the 0.95 mix ran on ume without lawsonite, and cold UM capacity stepped up between p = 0.95 and 1.

Each bulk carries excess water (80 wt% for the peridotites, 800 wt% for basalt, on the oxide scale) under the `aH2O` buffer, so a pure-water fluid is always present. Melt (`liq`) and the DEW aqueous fluid are removed. The capacity is (1.3) over the solid.

Why melt is off: hydrous phases above the wet solidus are removed by StagYY's own melting routines, and the melt's solubility is solved separately (TOex_math S.10, S.16). Keeping the melt here would double-count that water.

**ume hand-over (`Phand` = 5-6 GPa).** `ume` has no phase A, so above ~5.5 GPa at cold T it stores water in brucite + opx instead. From 5 GPa the sb24 + DHMS value (S.10) is computed on the same nodes, and the table blends linearly in wt%:

```julia
C = (1 − w) C_ume + w C_sb24+DHMS,   w = min((P − 5)/(6 − 5), 1)                           (3.1)   # UM hand-over of the peridotite blocks; pure sb24 + DHMS from 6 GPa to the 7 GPa seam.
```

*Symbols.*
- C_ume: capacity from the ume minimisation (1.3).
- C_sb24+DHMS: capacity from sb24 with the NAM laws and the DHMS term (S.10) at the same node.
- w: the linear weight in P (GPa).

The mbe hand-over uses the same form with 9.5 and 10 GPa. The host switches at w = 0.5.

**mbe compositions (`Pmbe` = 9.5-10 GPa).**
- They stay on `mbe` through the whole UM sector and into the TZ sector, then blend into sb24 over 9.5-10 GPa with the same linear rule.
- `mbe` keeps lawsonite past 14 GPa at cold T, so the hand-over, not the database, sets where cold basalt releases its lawsonite water. The window is set by experiment: lawsonite is stable in MORB + H2O to 10 GPa below 700 °C (Okamoto & Maruyama 1999).
- Two caveats:
  - above ~9 GPa `mbe` holds coesite instead of stishovite and its garnet has no majorite;
  - its lawsonite-out temperature is ~80-200 K hotter than the experiments (830 °C against 750 °C at 5.5 GPa, 930 °C against < 900 °C at 8 GPa).
- The seam at 7 GPa is `mbe` on both sides at the same P, so it is exact.

Code: `um_run` and the `:sH2O_um` block of `sfO2.jl`; the `mbe` rows at the end of the `:sH2O_deep` block.

## S.4 Nominally anhydrous minerals: the saturation laws

<img src="md_figs/S04_nam_laws.png" alt="S.4" width="1100">

*Figure S.4. Every NAM law at a = 1 against P at four temperatures. Grey band: the P range of the law's box; the T range is printed. Dashed lines are akimotoite (= ringwoodite/21) and the Al-bearing post-perovskite term; the post-perovskite panel applies (4.4b) to the 30 ppm bridgmanite reference.*

Most laws have the Keppler form, with f the PS fugacity (S.2):

```julia
c = A f^n exp(−(ΔH + P ΔV)/RT)                                                              (4.1)   # NAM water solubility at a = 1: Keppler power law in fugacity with an Arrhenius term.
```

*Symbols.*
- c: water content (wt%; the laws give ppm, times 10⁻⁴).
- A: pre-factor.
- f: water fugacity (S.2).
- n: fugacity exponent. n = 1 for molecular or two-proton defects; 0.5 for a single hydroxyl defect.
- ΔH: reaction enthalpy.
- ΔV: reaction volume.
- R: gas constant.
- P, T: pressure and temperature.

A law is evaluated at (P, T) clamped to its calibration box, so it is held flat outside that box. The boxes follow the calibration constraints of the source data (decision recorded in S.18). At a < 1 the solid's content is C·aⁿ, applied in StagYY (TOex_math S.15).

| phase(s) | law and constants | source | box P (GPa), T (K) |
|---|---|---|---|
| olivine | c = 10⁻⁴ exp(0.5·0.6447 ln f + 4905.5403/T) wt%, f in GPa | Dong et al. (2021) fit | 0.01-14.2, 1200-2400 |
| opx, pure En | A = 0.01354 ppm/bar, n = 1, ΔH = −4563 J/mol, ΔV = 12.1 cm³/mol | Mierdel et al. (2007) (En term, Rauch & Keppler 2002) | 0.01-14, 500-1500 |
| opx, Al term | A = 0.042, n = 0.5, ΔH = −79685, ΔV = 11.3; weight min(X_Al2O3/0.1005, 1) | Mierdel et al. (2007) | 0.01-14, 1073-1500 |
| cpx | x_jd·c_jd + (1 − x_jd)·c_di, x_jd the jadeite end-member fraction | | |
| diopside | A = 2.15, n = 0.5, ΔH = 0, ΔV = 7.43 | Bromiley et al. (2004) | 0-10, 873-973 |
| jadeite | A = 7.144, n = 0.5, ΔH = 0, ΔV = 8.019 | Bromiley & Keppler (2004) | 0-10, 873-973 |
| garnet, majorite | A = 0.866, n = 0.5, ΔH = 2580, ΔV = 5.71 | Lu & Keppler (1997) form | 0.01-10, 700-1500 |
| coesite, quartz | c = max(−49 + 6P + 0.06T, 0) ppm | linear fit | 0-10, 800-1073 |
| wadsleyite | in-house fit of Dong + Litasov data (22 points), n = 1 | this work | data extent 14-20, 1473-2173 |
| ringwoodite | in-house fit of Dong et al. data (6 points), n = 1 | this work | 19-23, 1723-2000 |
| akimotoite | ringwoodite / 21 | Bolfan-Casanova et al. (2000), D(rw/ak) = 21 | as ringwoodite |
| HP clinoenstatite | in-house fit (6 points), n = 0.5, ΔV in Pa form | this work | 8-13.4, 1373-1673 |
| stishovite | A = 24 ppm/GPa^0.5, n = 0.5, ΔH = −3065, ΔV = 4.29 (Al-free) | "Chen et al. 2021", as labelled in `sf_H2O.jl` | 0.01-60, 800-1500 |
| CaCl2-type silica | in-house fit (11 points), n = 0.5 | this work | 25-100, 2900-3500 |
| α-PbO2 silica | in-house fit (5 points), n = 0.5 | this work | 128-144, 3650-4100 |
| post-perovskite | c_ref / D(P, T), D = C_pv/C_ppv fitted to the Al-free and Al-bearing grids; weight X_Al2O3 (4.4) | Townsend et al. (2016) | 90-135, 1000-4000 |
| bridgmanite | S.5 | Litasov et al. (2003) | - |
| ferropericlase | 50 ppm | constant | - |
| CF phase | 30 ppm | constant | - |
| Ca-perovskite | 0 | constant | - |
| cristobalite | 1 ppm | constant | - |

Units: A is in ppm/bar^n with f in bar, ΔH in J/mol and ΔV in cm³/mol (entered as m³/mol times P in Pa).

Table symbols:
- X_Al2O3: the Al2O3 mole fraction of the opx on the oxide basis (`Comp`); 0.1005 is that of Al-saturated opx.
- x_jd: the jadeite (or jdm) end-member fraction of the cpx solution.
- In the coesite row, P is in GPa, T in K, and the result in ppm.
- In the olivine row, f is in GPa.

**In-house fits.** These use the fugacity form with P in GPa:

```julia
c = p₁ f(P, T)^n exp(−(p₂ + p₃ P)/RT)                                                        (4.2)   # In-house fits (nam_fits), least squares (LsqFit) on the tabulated data.
```

*Symbols.*
- p₁: pre-factor (wt%/GPa^n; ppm/bar^0.5 for HP clinoenstatite).
- p₂: enthalpy (J/mol).
- p₃: volume term (J/mol/GPa, i.e. 10³ cm³/mol; m³/mol in the Pa form of HP clinoenstatite).
- f: fugacity in GPa.

The fitted parameters (p₁, p₂ in J/mol, p₃ in J/mol/GPa) are:

| fit | p₁ | p₂ | p₃ |
|---|---|---|---|
| wadsleyite | 5.90·10⁻⁵ | −90047 | 10679 |
| ringwoodite | 6.47·10⁻⁵ | −56616 | 8647 |
| CaCl2-st | 9.63·10⁻³ | 10423 | 3035 |
| α-PbO2 | 3.51·10⁻³ | 284300 | 54.3 |
| HP clinoenstatite (Pa form, ppm) | 0.1758 | −37167 | −2.40·10⁻⁷ m³/mol |

The CaCl2 fit starts at [10⁻², 10⁴, 3·10³]; from ones(3) it lands in a negative-A minimum.

**Stishovite branches.** Below the CaCl2-type transition, which is tabulated from Umemoto et al. (2006) and extended linearly, stishovite uses the low-P law. Above it, it uses the CaCl2 fit:

```julia
P_b(T) = linear interpolation of the Umemoto (T, P) points;   c_st = P ≤ P_b ? c_Chen : c_CaCl2   (4.3)   # Stishovite law switches at the extrapolated CaCl2-type boundary.
```

*Symbols.* P_b(T) is the CaCl2-type transition pressure (GPa) at T (K). c_Chen is (4.1) with the stishovite row. c_CaCl2 is (4.2) with the CaCl2 fit.

The switch is a step of about two decades in the basalt lower mantle (Fig. S.10). It is kept: the CaCl2 phase is measured to hold far more water, and no data constrain a transition width (S.18).

**Post-perovskite.** Townsend et al. (2016) give the bridgmanite/post-perovskite partition coefficient D = C_pv/C_ppv from ab initio free energies, for an Al-free and an Al-bearing system. Their D grids are fitted with

```julia
ln D(P, T) = (p₁ + p₂P)/T − (p₃ + p₄P + p₅T + p₆)                                                  (4.4a)   # Fit to the digitised Townsend grids, one parameter set per system.
c_ppv = c_ref [x_Al / D_Al(P, T) + (1 − x_Al) / D_noAl(P, T)]                                          (4.4b)   # Post-perovskite water: coexisting-bridgmanite reference divided by the partition coefficient.
c_ref = Al_pv ? c_pv(a_ppv) (5.1) : 30 ppm                                                            (4.4c)   # Reference: the bridgmanite law at the ppv's own Al2O3, or the 30 ppm constant.
```

*Symbols.*
- D_Al, D_noAl: the fitted partition coefficients (P in GPa, T in K, clamped to 90-135 GPa and 1000-4000 K).
- p₁ … p₆: their fit parameters.
- x_Al: the Al2O3 mole fraction of the post-perovskite (oxide basis).
- a_ppv: its Al2O3 in wt%.
- c_ref: the water of a bridgmanite of the same Al content.

Near the geotherm, D_noAl ≈ 3-5 (Al-free post-perovskite holds 1/3-1/5 of the bridgmanite's water) and D_Al ≈ 0.3 (it holds ~3×), as Townsend et al. find. Fixed on 2026-09-27: the law previously multiplied by D.

**Checks against the literature.** These come from the 2026-09-27 review; the boxes were compared against Guimond et al. (2023, Table 1):
- **opx Al term:** its data span 1.5-3.5 GPa and 1073-1373 K. With ΔH = −80 kJ/mol it grows steeply below its lowest run temperature (1.9 wt% at 1 GPa/1000 K, 12 wt% at 800 K if extrapolated), so its T box starts at 1073 K, not 500 K like the En term. On the production maps the floor changes only UM nodes below 1073 K: 415 nodes by more than 1 %. The 17 that lose more than 0.1 wt% are opx-hosted, near the surface, and held up to 1.9 wt% before.
- **Garnet:** its data span 5-9 GPa and 1373-1473 K. The law varies only ×2 over 700-1600 K.
- **Diopside and jadeite:** the 873-973 K box holds almost every mantle node at its 973 K value (single line in the figure).
- **Wadsleyite and ringwoodite:** the boxes are the data extent, so the laws are flat below 1473 K and 1723 K.
- **Constants:**
  - ferropericlase 50 ppm lies inside the FTIR range of 10-80 ppm;
  - CF phase 30 ppm is below the measured < 100 ppm (FTIR, 26-32 GPa);
  - akimotoite = ringwoodite/21 is the measured partitioning (Bolfan-Casanova et al. 2000; the end member holds ~280 ppm);
  - Ca-perovskite 0 is a conservative floor against disputed hundreds-to-thousands of ppm.

Code: `nam_fits`, `nam_laws` (scalar and threaded vector methods), `NAMLAW`, `DRY`, the `:cpx`, `:opx`, `:pv`, `:ppv` branches of `assemble!`.

## S.5 Bridgmanite and the Al_pv flag

<img src="md_figs/S05_bridgmanite.png" alt="S.5" width="1100">

*Figure S.5. (a) The Al_pv law with the Litasov et al. (2003) data; bands show the Al2O3 range of the sb24 bridgmanite of harzburgite, basalt and the 20 % mix. (b) Al2O3 in sb24 bridgmanite with P. (c) The resulting bridgmanite water. (d) Rock capacity along 2200 K with the flag off (thin) and on (thick).*

**The law.** With `Al_pv = true`, bridgmanite water is a function of its own Al2O3 content a (wt%), read from the sb24 phase composition at each node. It is linear between three anchors and flat beyond them:

```julia
c_pv(a) = linear through (0 wt%, 100 ppm), (5.5 wt%, 1400 ppm), (15 wt%, 75 ppm); a clamped to [0, 15]   (5.1)   # Bridgmanite water against its Al2O3 (BDG_AL, BDG_W).
```

*Symbols.* a is the Al2O3 content of the bridgmanite (wt%, `Comp_wt` × 100). c_pv is its water (wt%; 100 ppm = 0.01 wt%). The anchors are `BDG_AL` = (0, 5.5, 15) and `BDG_W` = (0.01, 0.14, 0.0075), and the law is `bdg_h2o`.

The anchors are the three compositions of one study (Litasov et al. 2003, SIMS):

| composition | Al2O3 (wt%) | H2O (ppm) |
|---|---|---|
| pure MgSiO3 | 0 | ~100 |
| Al-Mg-pv | 4-7 | 1000-1500 |
| peridotite Al-Fe-Mg-pv | 5-6 | 1400-1800 |
| MORB pv | 13-17 | 40-110 |

The rise and fall follows the substitution mechanism:
- At low Al, Al³⁺ enters with an oxygen vacancy (MgAlO2.5), which takes up OH.
- At high Al the charge-coupled Al-Al and Fe³⁺-Al substitutions dominate, and these take up almost none.

**Why no Fe³⁺ factor.** A factor (1 − f_Fe³⁺Al) was considered and dropped. At R = 0.01, sb24 puts little Fe³⁺ into hot basalt bridgmanite, so that factor would leave it wet, against the measured 40-110 ppm.

**Why no pressure term (decision).** The vacancy component is reported to fall with pressure, but its depth scale is not constrained well enough to fix. Composition alone sets the law.

Behaviour on the maps, at 2200 K:
- Al in sb24 bridgmanite rises through the top of the LM as garnet dissolves, then stays flat: harzburgite ~1.3 wt%, the 20 % basalt mix ~4.5 wt% from 27 GPa, hot basalt 21-25 wt%.
- The water follows: harzburgite ~400 ppm, the 20 % mix ~1170 ppm (near the peak), hot basalt 75 ppm.
- Cold basalt climbs from 8 to 25 wt% Al2O3 over 30-75 GPa and passes the falling limb (8-15 wt%). That limb lies between anchors and is the least constrained part.
- Rock capacity along 2200 K: harzburgite LM rises from 0.0034 to 0.032 wt%, the 20 % mix from 0.003 to 0.089 wt%, basalt by 10-50 % (stishovite dominates).
- The Al-free end (100 ppm) also serves the Al-free bridgmanite of the MSH re-equilibration (S.10).

With `Al_pv = false` every bridgmanite gets 30 ppm. That is the previous behaviour, reproduced bit for bit.

The literature spread is wide: from ~10 ppm without Al (Bolfan-Casanova et al. 2000, 2003), through ~1000 ppm (Fu et al. 2019) and 150-300 ppm ("nearly dry", Liu et al. 2021), to a 2026 Cryo-SIMS report of up to 5500 ppm (Semba et al., not verified here). The law sits inside that spread and is anchored to one internally consistent data set.

Code: `BDG_AL`, `BDG_W`, the `law === :pv && alpv` branch of `assemble!`, `nam_laws(...; alpv)`.

## S.6 The MgO-SiO2-H2O model of Komabayashi & Omori (2006)

<img src="md_figs/S06_ko_benchmark.png" alt="S.6" width="1100">

*Figure S.6. (a) Our minimiser in MgO-SiO2-H2O at Mg/Si = 1.4 (K&O's phA + En bulk), water-saturated at the fitted activity, with K&O's 50 digitised reaction checkpoints by class. (b) ΔP of our boundary at the checkpoint's T and (c) ΔT at its P, ours minus K&O's; for dehydration reactions, filled at K&O's own a(H2O) and open at our activity map (S.7), joined per checkpoint. Shaded band: the 0.2 GPa / 20 K tolerance of the benchmark counts below.*

**Phases.** The dense hydrous magnesium silicates (DHMS) and their anhydrous partners are taken from K&O Tables 3-4. There are 16 rows plus phase H (S.8), and each row carries its MgO, SiO2 and H2O per formula (m, s, h):

| phases | type |
|---|---|
| phA, cHu, phE, phD, shB, br, phH | DHMS (water counted as hydrous) |
| hyWad, hyRin | hydrous NAM (water through the wa and ri laws, S.10) |
| fo, hen, wd, rg, aki, pv, per, stv | anhydrous |

**Gibbs energy.** The apparent Gibbs energy at (P, T) builds on H, S and Cp from 298.15 K. The volume is expanded at 1 bar and integrated with a temperature-dependent Murnaghan equation:

```julia
Cp = a + bT + cT⁻² + dT² + eT⁻³ + fT^{−1/2} + gT⁻¹                                             (6.1)   # K&O heat capacity; H(T) and S(T) are its closed-form integrals.
V₁(T) = V₀ exp[α₀(T − T₀) + (α₁/2)(T² − T₀²) + α₂ ln(T/T₀) − α₃(1/T − 1/T₀)]                    (6.2)   # 1-bar volume from the K&O thermal expansion polynomial.
K(T) = 1/(β₀ + β₁T + β₂T² + β₃T³),   K′(T) = K′₀ + K′_T (T − T₀) ln(T/T₀)                      (6.3)   # Isothermal bulk modulus (bar) and its pressure derivative.
G(P, T) = H(T) − T S(T) + V₁K/(K′ − 1) [(1 + K′(P − 1)/K)^{1 − 1/K′} − 1]                        (6.4)   # Apparent Gibbs energy; Murnaghan ∫V dP from 1 bar (P in bar).
```

with the closed-form integrals of (6.1):

```julia
H(T) = ΔfH + a(T − T₀) + (b/2)(T² − T₀²) − c(1/T − 1/T₀) + (d/3)(T³ − T₀³) − (e/2)(T⁻² − T₀⁻²) + 2f(√T − √T₀) + g ln(T/T₀)   (6.5)   # Enthalpy from ∫Cp dT.
S(T) = S° + a ln(T/T₀) + b(T − T₀) − (c/2)(T⁻² − T₀⁻²) + (d/2)(T² − T₀²) − (e/3)(T⁻³ − T₀⁻³) − 2f(T^{−1/2} − T₀^{−1/2}) − g(1/T − 1/T₀)   (6.6)   # Entropy from ∫Cp/T dT.
```

*Symbols.*
- ΔfH: formation enthalpy from the elements at T₀ = 298.15 K and 1 bar (J/mol).
- S°: third-law entropy at T₀ (J/mol/K).
- a … g: the K&O heat-capacity coefficients (J/mol/K with the powers of T shown). These a … g are unrelated to the letters of S.2.
- V₀: the volume at T₀ and 1 bar (J/bar; ×10 gives cm³/mol).
- α₀ … α₃: thermal-expansion coefficients.
- V₁(T): the 1-bar volume at T.
- β₀ … β₃: compressibility coefficients (1/bar).
- K(T): the isothermal bulk modulus (bar).
- K′₀, K′_T: the pressure derivative of K and its temperature coefficient.
- P: pressure in bar.
- G: the apparent Gibbs energy (J/mol).
- m, s, h: moles of MgO, SiO2 and H2O per formula unit.

Printed coefficients are rescaled on entry (`ko_row`). A row with V₀ = 0 (water) returns H − TS.

**Benchmark.** Reaction positions are recomputed with our G and our μ_H2O:
- Ohtani's two calibration points of R87 are reproduced within 0.01 GPa.
- Fluid-absent reactions: median |ΔP| = 0.04 GPa.
- Dry reactions (Fabrichnaya 2004 data): median 0.07 GPa.
- Dehydration reactions at K&O's own activities: 31 of 33 fall within 0.2 GPa or 20 K.
- With our activity map (S.7): 28 of 33, median ΔP 0.08 GPa and median ΔT 7 K.

The water contents of K&O Fig 7 at Mg/Si = 1.4 come out as 1.22, 1.36, 0.73 and 1.74 wt%, against 1.25, 1.37, 0.75 and 1.78. The 0.86 field is missed because the activity is extrapolated there; no capacity effect is involved.

Code: `ko_row`, `KO`, `KO_NAME`, `KO_DHMS`, `KO_LAW`, `g_ko`; the benchmark in `sH2O_plot.jl` (S06 block).

## S.7 The water-activity map

<img src="md_figs/S07_activity.png" alt="S.7" width="1100">

*Figure S.7. a_H2O of the DHMS minimisation on a log colour scale; contours and anchors are labelled in a: a = 1 from Komabayashi et al. (2005) (circles), the activities Komabayashi & Omori (2006) fitted to fluid-present experiments (diamonds), and the activities that put the Nishi et al. (2014) phase D and phase H dehydration lines on this model (squares). Dashed box: the region where the map varies (held flat outside).*

**Why a map.** K&O reproduce fluid-present experiments with water activities below 1: the fluid carries dissolved silicate, and the activity falls with P and T. We make this a smooth field. A thin-plate spline in ln a is fitted exactly through the anchors, in scaled coordinates x = P/5 and y = T(°C)/200:

```julia
ln a(x, y) = c₀ + c₁x + c₂y + Σᵢ wᵢ φ(rᵢ),   φ(r) = r² ln r,   rᵢ = |(x, y) − (xᵢ, yᵢ)|             (7.1)   # Thin-plate spline through the anchors (ah2o_fit); a linear system with the three affine side conditions.
a = exp(min(ln a(clamp(P, 5, 60)/5, clamp(T°C, 500, 1600)/200), 0))                            (7.2)   # Evaluation: held flat outside 5-60 GPa and 500-1600 °C, capped at 1 (ln_ah2o).
```

*Symbols.*
- a: water activity of the fluid.
- (x, y): scaled coordinates of the evaluation point.
- (xᵢ, yᵢ): the scaled coordinates of anchor i, whose target is ln aᵢ.
- rᵢ: the scaled distance to anchor i.
- φ: the thin-plate basis.
- wᵢ: the anchor weights.
- c₀, c₁, c₂: the affine part.

The weights and affine coefficients solve

```julia
[Φ  Q; Qᵀ  0] [w; c] = [ln a_anchor; 0],   Φ_ij = φ(|p_i − p_j|),   Q = [1  x  y]                       (7.1a)   # Interpolating thin-plate spline: exact at the anchors, minimum bending energy.
```

Here p_i = (x_i, y_i), Φ is the n × n matrix of basis values between anchors, Q the n × 3 affine matrix, w and c = (c₀, c₁, c₂) the unknowns, and n = 32 the anchor count. Clamping (7.2) holds a flat outside the box; the cap ln a ≤ 0 stops the spline overshooting above pure water.

**Anchors (32).**
- **K&O fitted activities:** 18 labelled points, 12-31 GPa, a = 0.012 to 0.85.
- **a = 1 anchors:** 8 points on the cold side: 5-9 GPa up to 800 °C (Komabayashi et al. 2005) and 11.7-25 GPa below 750 °C.
- **Six Nishi et al. (2014) anchors,** computed so that their measured dehydration lines lie on this model:

```julia
a_D(P, T) = exp{−[(236 G_pv + 80 G_stv − 200 G_phD)/226 + μ_H2O]/RT}    at (30, 1187), (40, 1105), (48.6, 1035) GPa, °C   (7.3)   # Phase D dehydration (K&O R88) placed on Nishi's line.
a_H(P, T) = exp{[G_phH − G_pv − μ_H2O]/RT}                              at (40, 997), (50, 1030), (57.9, 1056)            (7.4)   # Phase H -> bridgmanite + H2O placed on Nishi's line.
```

*Symbols.*
- G_pv, G_stv, G_phD, G_phH: apparent Gibbs energies (6.4) at (P, T).
- μ_H2O: the a = 1 potential (2.7).

(7.3) is K&O's reaction R88, 200 phD = 236 pv + 80 stv + 226 H2O, written per mole of H2O. (7.4) is MgSiO2(OH)2 = MgSiO3 + H2O. Each gives the activity at which that dehydration is at equilibrium at Nishi's measured (P, T).

These give a = 0.053, 0.011 and 0.0029 (phase D) and 0.0077, 0.0024 and 0.0012 (phase H). The phase H anchors depend on phase H itself, so the map and phase H are iterated three times (S.8).

**Tension at 30-37 GPa.** K&O's 0.03 at 31 GPa/1040 °C and Nishi's 0.053 at 30 GPa/1187 °C pull the spline in opposite directions. Phase D dehydration in Si-rich bulk dips ~100 K below Nishi's line at 33-37 GPa.

Code: `KO_AH2O`, `ah2o_fit`, `ln_ah2o`, `msh_setup`.

## S.8 Phase H

<img src="md_figs/S08_phaseH.png" alt="S.8" width="1100">

*Figure S.8. Our minimiser in Nishi et al.'s two bulks with their boundaries (dashed) and phase H observations (squares); Ohtani et al. (2014)'s diamond-cell phase-H-out point (triangle).*

Phase H, MgSiO2(OH)2 (m = s = h = 1), is not in K&O. It is calibrated in the K&O manner, from phase-equilibrium brackets:

1. **Heat capacity and thermal properties.**
   - Cp = Cp_pv + Cp_br − Cp_per (oxide sum).
   - α(T) is phase D's.
   - K_T(T) is phase D's shape, rescaled to K₀ = 196.07 GPa at 298 K, the value fitted to the measured cells of Nishi 2014 and Ohtani 2014.
2. **ΔfH and S.** These are exact on two points of Nishi's phD = phH + St line: (45.4 GPa, 600 °C) and (48.6 GPa, 1035 °C). The reaction is written for K&O's phase D composition:

```julia
phD = 1.18 phH + 0.40 stv + 0.05 H2O                                                           (8.1)   # Balanced with K&O's phase D (m 1.18, s 1.58, h 1.13).
ΔG_8.1(P_k, T_k) = 0,  k = 1, 2   ⇒   S = (g₂ − g₁)/(T₂ − T₁),  ΔfH = T₁S − g₁                     (8.2)   # G is linear in (ΔfH, S): two points fix both exactly.
```

*Symbols.*
- (P_k, T_k): the two calibration points.
- g_k: ΔG of (8.1) evaluated with ΔfH = S = 0 for phase H, divided by 1.18. Since G_phH = G_phH(ΔfH = S = 0) + ΔfH − TS, the condition ΔG = 0 at both points is ΔfH − T_k S = −g_k.
- The H2O term uses μ = μ_H2O + RT ln a.

3. **V₀ and K′.** These are fitted (LsqFit) so that phase H first appears against the phH-free assemblage of its own bulk (Mg/Si = 1). The targets are three points of Nishi's phH-in curve: (28, 630), (31.5, 800) and (34, 900) GPa, °C.
   The misfit at each target is G_phH − μ − (n_i G*_i + n_j G*_j), where (i, j, n_i, n_j) is the phH-free minimum of (9.2) at Mg/Si = 1 and G* is (9.1). V₀ (J/bar) and K′ enter G_phH through (6.2)-(6.4). ΔfH and S are re-solved from (8.2) inside every evaluation.
4. **Iteration.** The a_H2O anchors (7.4) and steps 2-3 are iterated three times; the result converges to the printed digits.

Result:
- ΔfH = −1721.19 kJ/mol and S = 71.57 J/mol/K.
- V₀ = 3.2293 J/bar (32.29 cm³/mol; Bindi et al. 2014 measure 35.03 at ambient) and K′ = 8.04.
- V at 42 GPa/1273 K is 28.74 cm³/mol against 29.79 measured (−3.5 %).

**Why an effective volume.** With the measured volume and K&O's phase D equation of state, ΔV of phD → phH + St is +1.3 cm³ at 45 GPa, so the reaction runs backwards. The fitted volume absorbs the extrapolation of K&O's phase D and brucite equations of state beyond their 35 GPa calibration. Phase H's volume is used only here, never in a density table.

**Benchmark.**
- phH-in in its own bulk lies at 28.16, 30.01, 34.78 and 40.15 GPa, against 28, 31.5, 34 and 39.5.
- The phase D and phase H dehydration lines and phD = phH + St are exact (they are calibration inputs).
- Ohtani et al.'s phase D + pv + pc stays on the phase D side (+9.9 and +5.4 kJ).
- 8 of Nishi's 9 in-situ observations fall inside the phase H field; 34.4 GPa/900 °C misses by 0.4 GPa.

**Limits.**
- **PDmax = 50 GPa.** Phase D is excluded beyond this (G = 10¹² in the minimisation); Nishi sees none above 48.6 GPa. Without the cut, the extrapolated phase D returns at ~55 GPa and 950-1050 °C in Si-rich bulk.
- **PHmax = 60 GPa.** No DHMS beyond this: Ohtani et al. (2014) find MgSiO4H2 phase H stable to 60 GPa (and to 1610 K), and Tsuchiya & Umemoto (2019) compute its dissociation into bridgmanite + ice VII at ~60 GPa near 1000 K (52 GPa without zero-point energy, Tsuchiya 2013). Cold slabs therefore release their DHMS water (~11 wt% in harzburgite) in one layer at 60 GPa, as the dissociation does. Al-bearing phase H, which survives deeper, is out of scope with the other Al phases (S.10).
- **TDmax = 1873 K,** K&O's calibration limit. No DHMS survive near it, so this cut leaves no step.

**Dehydration temperature: a dataset choice.** Phase H's dehydration is calibrated to Nishi et al. (2014, multi-anvil, ~1000-1060 °C). Ohtani et al. (2014, diamond cell) place phase-H-out ~350 K hotter at 50 GPa. Anchoring on Ohtani instead would widen the phase H field by that much. The choice follows Nishi because the phase D, phase D → H and phase H lines then come from one consistent study.

Code: `ko_phH`, `msh_setup`, the `PDmax` test in `dhms!`.

## S.9 The DHMS term: the question and the minimisation

**The question.** sb24 has no hydrous phases. At a cold node it therefore returns an anhydrous assemblage, and the NAM laws give that assemblage's water. Harzburgite at 18 GPa and 1100 K, for example, is ringwoodite + stishovite + garnet holding 1.21 wt%. Given free water, the Mg-Si silicates of that rock would react to phase D + superhydrous B and hold ~8 wt%.

The DHMS term asks, at every node separately: *if the Mg-Si part of this rock could react with free water, what would it become, and how much water would it then hold?*

It is answered in seven steps. Steps 1-3 and 6-7 are bookkeeping on the sb24 result (S.10); steps 4-5 are the minimisation (this section).

```text
one node: (P, T) of one sector, one of the ten compositions
 1. sb24 minimisation of the dry bulk at R = 0.01          → solid phases φ and their fractions x_φ
 2. water of that assemblage, (1.3)                        → C_sb24 (the value with dhms = false)
 3. split the solid into the MSH part and the rest, (10.1) → N_M, N_S (moles), m_MS (mass), w_MS (its water)
 4. water potential of the free fluid at (P, T)            → μ = μ_H2O (2.7) + RT ln a (7.2)
 5. minimise G over the MgO-SiO2 phases at that μ, (9.2)   → the best pair of phases (i, j) and amounts n_i, n_j
 6. no DHMS in the pair → keep C_sb24 unchanged
    a DHMS in the pair  → rebuild the node's water, (10.2)-(10.5)
 7. store C, the hydrous share and the host at the node    → maps → hand-over blends (S.3) → StagH2O.dat
```

Nodes never exchange information: every node is minimised and rebuilt from its own (P, T, bulk) only. The worked example at the end of S.10 follows two nodes through all seven steps.

**Step 4: why water is not a constraint.** A capacity is measured at water saturation, so free water is always present in excess. The rock can take as much water as it likes, and every mole it takes is drawn from the fluid at the fluid's chemical potential μ. Water is therefore exchanged with a reservoir, not distributed, and it drops out of the mass balance. Each phase is instead charged for the water it takes. That charge is its transformed energy:

```julia
G*_φ = G_φ(P, T) − h_φ μ                                                                          (9.1)   # Legendre transform: energy of a phase after releasing its water to the reservoir at μ.
```

A hydrous phase wins only if its G is lower than its anhydrous competitors by more than the value of the water it binds, h_φ μ. The activity a < 1 of the map (S.7) lowers μ and so makes the hydrous phases harder to stabilise; that is how K&O's dehydration temperatures are reproduced.

**Step 5: the minimisation.** What is left is a two-component problem. Distribute the rock's N_M moles of (Mg,Fe)O and N_S moles of SiO2 among the 17 phases so that the total transformed energy is lowest:

```julia
min Σ_φ n_φ G*_φ   s.t.   Σ n_φ m_φ = N_M,   Σ n_φ s_φ = N_S,   n_φ ≥ 0                         (9.2)   # Linear programme over the 17 phases in (MgO, SiO2).
```

*Symbols.*
- φ: the 17 K&O phases (S.6, S.8).
- G_φ: the apparent Gibbs energy (6.4).
- h_φ, m_φ, s_φ: H2O, MgO and SiO2 per formula unit.
- μ: the water potential of the fluid, μ_H2O (2.7) plus RT ln a (7.2).
- G*_φ: the transformed energy (J/mol).
- n_φ ≥ 0: the moles of phase φ, per gram of the node's solid.
- N_M, N_S: the moles of (Mg,Fe)O and SiO2 to be distributed, per gram of solid (10.1).

(9.2) is linear in the unknowns n_φ, with two equality constraints. The lowest value of a linear function over such a set is reached at a corner, and a corner has at most two nonzero n_φ. This is the phase rule: two components at fixed P, T and μ give at most two phases.

`msh_min` therefore tries every pair of phases, 136 pairs for 17 phases. For each pair it finds the unique amounts that reproduce the bulk, and it keeps the cheapest pair whose amounts are both non-negative:

```julia
det = m_i s_j − m_j s_i;   n_i = (N_M s_j − N_S m_j)/det;   n_j = (m_i N_S − s_i N_M)/det            (9.3)   # The 2 × 2 solve for pair (i, j); pairs with det = 0 are skipped.
keep if n_i, n_j ≥ −10⁻¹²; answer = argmin n_i G*_i + n_j G*_j                                    (9.4)   # Feasible pairs only; single-phase answers appear as pairs with one n = 0.
```

*Symbols.*
- (i, j): a pair of phases.
- det: the determinant of their 2 × 2 composition matrix.
- (n_i, n_j): the amounts that reproduce (N_M, N_S) exactly.
- The −10⁻¹² tolerance admits round-off; negative amounts inside it are set to zero.

A pair with a negative amount cannot make the bulk (it would need a negative quantity of a phase), so it is discarded. Every pair that survives is a possible assemblage, and the cheapest one is the answer. Nothing is iterated and there is no tolerance on the result: the answer is exact for the node's (N_M, N_S, μ). Phase D is removed above PDmax by setting its G* to 10¹².

<img src="md_figs/S09_msh_hull.png" alt="S.9" width="1100">

*Figure S.9. G* per mole of oxide against SiO2/(MgO + SiO2) for all 17 phases, relative to the per-stv chord, at two conditions. Orange: DHMS. The lower convex hull is the minimiser's answer for every bulk; the vertical lines are two bulk compositions.*

**The same answer as a picture (Fig. S.9).** Plot every phase at x = s/(m + s), its Si fraction, and y = G*/(m + s), its transformed energy per mole of oxide. Any mixture of two phases lies on the straight line between them. The cheapest way to make a bulk of Si fraction x_b is therefore the lowest line that passes over x_b. That line is a segment of the lower convex hull of the points, and its two ends are the answer pair.

- A phase above the hull is never stable at that (P, T).
- A DHMS on the hull is stable over the bulk range of its two hull segments.
- As P and T change, the points move, and a phase reaches or leaves the hull. On the maps this shows as a phase boundary.

Code: `msh_min`, the `G` tuple in `dhms!`.

## S.10 From the minimisation back to the node and the maps

<img src="md_figs/S10_capacity.png" alt="S.10" width="1100">

*Figure S.10. Capacity of the three blocks (columns) in the three sectors (rows, pressure increasing downwards), with host phases, computed on 300 × 300 nodes per sector (the shipped table has 100 × 100). Contours: 0.01, 0.1, 1 wt%. Dashed: 50 and 60 GPa (phase D and DHMS limits); dotted: the hand-over windows.*

<img src="md_figs/S11_dhms_gain.png" alt="S.11" width="1100">

*Figure S.11. log10 of the capacity with the DHMS term over the capacity without it (`dhms = false`), per block and sector.*

**Steps 1-2: the dry node.** sb24 is minimised at the node for the dry bulk of the block (R = 0.01). `assemble!` then gives every solid phase its water by (1.3) and sums it to C_sb24. With `dhms = false` this is the final value of the node.

**Step 3: the MSH part.** The phases that the MSH system can represent are those of `MSH_SET`: ol, wa, ri, opx, hpcpx, ak, pv, mw, st, coe and qtz. In the same pass, `assemble!` sums what they hold:

```julia
N_M = Σ x_φ (w_MgO/40.30 + w_Fe/M_Fe),  N_S = Σ x_φ w_SiO2/60.08,  m_MS = Σ x_φ,  w_MS = Σ x_φ c_φ   (10.1)   # Moles of (Mg,Fe)O and SiO2 per gram of solid, their mass and their NAM water.
```

*Symbols.*
- The sums run over the MSH phases of the node only.
- x_φ: the weight fraction of phase φ in the solid.
- w_MgO, w_Fe, w_SiO2: the weight fractions of those components in phase φ (`Comp_wt`).
- M_Fe: the molar mass of the iron component, 55.85 for Fe and 71.85 for FeO.
- N_M, N_S: moles per gram of solid.
- m_MS: the mass fraction of the MSH part.
- c_φ: the NAM water of phase φ (S.4-S.5, wt%).
- w_MS: the MSH part's share of C_sb24 (wt% of solid).

Fe is counted with Mg, because K&O find that Fe, Al and Ca do not move the DHMS reactions above 14 GPa. Everything else (garnet, clinopyroxene, Ca-perovskite, the CF phase, native iron) is the *rest*. It is left exactly as sb24 found it, with its NAM water C_sb24 − w_MS.

**Steps 4-5.** μ is evaluated at the node, and (9.2) is solved for (N_M, N_S). The answer is one pair (i, j) with amounts (n_i, n_j) in moles per gram of solid.

**Step 6: keep or rebuild.** The node is rebuilt only if the pair contains a DHMS: phA, cHu, phE, phD, shB, br or phH.

- **No DHMS in the pair.** The node keeps C_sb24, its hydrous share and its host unchanged. This holds even when the anhydrous pair differs from sb24's assemblage. sb24 carries Fe, Al, Ca and Na, so its anhydrous answer is the better one; the MSH model is only asked whether water changes it.
- **A DHMS in the pair.** The MSH part of the node is replaced by the pair. The water is then recounted:

```julia
W   = Σ_{DHMS} n_φ h_φ 18.015                                        [g H2O per g anhydrous solid]   (10.2)   # Water bound in the DHMS of the re-equilibrated subsystem.
M_b = (m_MS − 60.08 N_S)/N_M                                                                        (10.3)   # Effective (Mg,Fe)O molar mass: carries Fe and the minor oxides of the subsystem.
w_N = Σ_{anhydrous} n_φ (m_φ M_b + 60.08 s_φ) c_law(φ)                                              (10.4)   # NAM water of the re-equilibrated anhydrous phases (KO_LAW: hen -> opx below 8 GPa, hpx above).
C   = (100 W + w_N + C_sb24 − w_MS)/(1 + W)                                                          (10.5)   # New capacity per gram of hydrous solid; the non-MSH phases keep their NAM water.
C_hyd = 100 W/(1 + W)                                                                                (10.6)   # Part held by hydrous phases.
```

*Symbols.*
- n_φ, h_φ: the pair's amounts (mol per g of solid) and water per formula.
- 18.015: the molar mass of H2O (g/mol).
- W: grams of water bound in the DHMS per gram of the dry solid.
- M_b: the mean molar mass of the MSH part's (Mg,Fe)O (g/mol).
- m_φ M_b + 60.08 s_φ: the dry molar mass of phase φ built from that (Mg,Fe)O.
- c_law(φ): the NAM law of an anhydrous pair member (`KO_LAW`: fo → ol, wd → wa, rg → ri, aki → ak, pv → pv, per → mw, stv → st, hen → opx below 8 GPa and hpx above; hyWad and hyRin → wa and ri).
- w_N: the NAM water of the anhydrous pair member (wt% of solid).
- C_sb24: the node's capacity before the rebuild (step 2).
- C: the new capacity (wt% of the hydrous solid).
- C_hyd: the part of C held in DHMS.

How (10.5) is built, per gram of dry solid:
1. **Remove** the MSH part's old NAM water, w_MS. It no longer exists: those phases have reacted.
2. **Add** the water of the pair: 100 W in the DHMS, and w_N in an anhydrous partner (for example periclase beside phase H).
3. **Keep** the rest's water, C_sb24 − w_MS.
4. **Divide** by the new mass, 1 + W. M_b is defined so that the pair's dry mass equals m_MS: N_M M_b + 60.08 N_S = m_MS. The rebuilt solid therefore weighs the old 1 g plus the bound water W, and (10.5) is in wt% of the hydrous solid, as (1.3) is.

hyWad and hyRin are treated as NAMs: they carry the wadsleyite and ringwoodite laws, not their small K&O h, because those laws already describe hydrous wadsleyite and ringwoodite at saturation.

**Step 7: the host and the hydrous share.** The host is the member with the largest water contribution: the DHMS in practice. The hydrous share is C_hyd/C. Both are stored for the maps only.

**Where each node goes.** Nodes are gathered into one list per pass: all (i, j, sector, composition) of TZ + LM in the `:sH2O_deep` pass, the hand-over rows in the `:sH2O_um` pass. MAGEMin, `assemble!` and `dhms!` each run once over the whole list, threaded over nodes. The results are then written back to their own node, `S[i, j, sector, composition] = C[k]`:

| sector and composition | how the node's value is made |
|---|---|
| UM, P < 5 GPa, all compositions | ume or mbe alone, which carry their own hydrous phases (S.3); no DHMS term |
| UM 5-6 GPa, ume compositions (p < `pmbe`) | linear blend of ume and sb24 + DHMS (3.1); from 6 GPa sb24 + DHMS alone |
| UM and TZ below 9.5 GPa, mbe compositions (p ≥ `pmbe`) | mbe alone; 9.5-10 GPa blend into sb24 + DHMS |
| TZ and LM, all other nodes | sb24 + DHMS, steps 1-7 |
| beyond 60 GPa, 1873 K, or N_M = 0 | sb24 + NAM laws only (steps 1-2) |

A boundary on the maps is where the winning pair changes between neighbouring nodes. Its position is resolved to the grid spacing: 0.07 GPa and 15 K in the UM, 0.18 GPa and 33 K in the TZ, 1.1 GPa and 33 K in the LM. The figures smooth these boundaries for display only; the table holds the node values.

**Worked example.** Harzburgite, two nodes, as printed by the code:

| step | 18 GPa, 1100 K (cold TZ) | 45 GPa, 1200 K (cold LM) |
|---|---|---|
| 1. sb24 solid | ri 0.903, st 0.041, gtmj 0.056 | pv 0.745, mw 0.231, capv 0.019, fee 0.005 |
| 2. C_sb24 | 1.213 wt% (ringwoodite law) | 0.032 wt% (bridgmanite law) |
| 3. MSH part | ri + st: N_M = 0.01234, N_S = 0.00685 mol/g (Mg/Si 1.80), m_MS = 0.944, w_MS = 1.212 | pv + mw: N_M = 0.01248, N_S = 0.00707, m_MS = 0.976, w_MS = 0.032 |
| 4. water potential | μ_H2O = −172.2 kJ, a = 0.88 | μ_H2O = +59.1 kJ, a = 0.0045 (RT ln a = −53.8 kJ) |
| 5. best pairs (Σ n G*, kJ per g) | **phD + shB −10.639**; phA + phD −10.620; phD + br −10.610 | **phH + per −5.097**; phD + per −5.065; pv + per −5.056 |
| 6. rebuild | W = 0.086: 0.00257 phD + 0.00093 shB mol/g; w_N = 0 | W = 0.127: 0.00707 phH + 0.00540 per mol/g; per holds w_N ≈ 0.001 |
| C = (100W + w_N + C_sb24 − w_MS)/(1 + W) | (8.58 + 0 + 0.0008)/1.086 = **7.90 wt%**, host phD | (12.74 + 0.001 + 0)/1.127 = **11.30 wt%**, host phH |

Reading the example:
- At 18 GPa, all 1.21 wt% of ringwoodite water is replaced by DHMS water. Only garnet's 0.0008 wt% survives from sb24.
- The three best pairs lie within 30 J per gram of rock. A small change in P, T or bulk Mg/Si moves the node to phase A + phase D. This is why the TZ maps are a mosaic of DHMS fields.
- At 45 GPa, the dry answer pv + per is only 41 J/g above phH + per. That margin is what the activity map (S.7) and the phase H calibration (S.8) set.

**Scope of the subsystem.**
- **Domain:** P ≤ PHmax (60 GPa), T ≤ TDmax (1873 K) and N_M > 0.
- **Not in the MSH part:** garnet, majorite, clinopyroxene, Ca-perovskite, the CF phase and every Al phase.
  - Al-bearing DHMS (phase egg, δ-AlOOH) are out of scope for this project, so basalt's DHMS are an MSH proxy built from the Mg of its pyroxene, akimotoite and bridgmanite together with its stishovite.
  - Majorite was tested and not adopted. Moving garnet's MgSiO3 excess into the MSH part changes harzburgite by < +0.05 wt% and the p = 0.3 mix by +0.5 to +1.3 wt% (phase D, 16-24 GPa). Basalt goes from ~0.02 to 4-5 wt%, but that is an upper bound: the test treats the majorite component as free and counts Ca- and Na-majorite as MgSiO3.

Code: `MSH_SET`, `assemble!(...; sub)`, `dhms!`; the `:sH2O_um` and `:sH2O_deep` blocks of `sfO2.jl`.

## S.11 The assembled capacity maps

<img src="md_figs/S12_hydrous_share.png" alt="S.12" width="1100">

*Figure S.12. Share of the capacity held by hydrous phases.*

<img src="md_figs/S19_composition.png" alt="S.19" width="1100">

*Figure S.19. Capacity against basalt fraction at six nodes: the ten table compositions (dots), with DHMS (solid) and without (dashed). Grey band: the silica-saturation step between pa and pb (StagfO2step, R = 0.01).*

<img src="md_figs/S14_T_profiles.png" alt="S.14" width="1100">

*Figure S.14. Capacity against T at fixed P, with DHMS (solid) and without (dashed).*

Production run, 100 × 100 per sector, ten compositions, `Al_pv = true`. The map figures (S.10-S.12, S.14) show harzburgite, basalt and the 20 % basalt mix. Median capacities per host phase, in wt%:

| block | UM | TZ | LM |
|---|---|---|---|
| harzburgite | atg 12.4, phA 5.6 (from 6 GPa), chl 0.66, opx 0.024, ol 0.05 | br 18, phE 10.7, phD 7.9, phA 5.7, shB 4.3; ri 1.05, wa 0.30, ol 0.03, gtmj 0.019 | phH 11.3 (31-60 GPa, < ~1300 K), phD 8.2; pv 0.032, ppv 0.010 |
| 20 % mix (p = 0.2) | atg 9.7, chl 2.8, phA 4.3, amp 0.26 | br 15.3, phE 8.3, phD 6.7, phA 4.4, shB 3.4; hyWad 1.3, ri 0.80, wa 0.29, gtmj 0.020 | phH 10.9, phD 7.9; pv 0.089, ppv 0.034 |
| basalt | law 4.5, chl 3.1, amp 0.46, g 0.016 | law 2.0 (to 9.5 GPa), phD 2.0 (cold, from ~13 GPa); gtmj 0.026 | phH 4.7, phD 4.1; st 0.40, apbo 0.28, pv 0.006 |

- **Share held by hydrous phases (Fig. S.12).** Close to 1 wherever a hydrous phase is stable, and 0 elsewhere. Mixed values appear only at field edges: the UM actinolite/talc and amphibole fields, and the cold TZ edge where hyWad sits beside a DHMS.
- **Hot mantle.** It is carried entirely by the NAM laws:
  - harzburgite: up to 0.4 wt% in UM olivine, up to 1.4 wt% ringwoodite in the TZ, 0.032 wt% in the LM;
  - basalt: 0.006-0.8 wt% in the LM, set by the silica phases (the CaCl2-type stishovite switch is the step at 55-110 GPa in Fig. S.10).
- **DHMS gain (Fig. S.11), median and range:**
  - UM above 6 GPa, ume compositions: ×48 (harzburgite) to ×43 (p = 0.5), to ×72; the mbe compositions (p ≥ 0.65) take no DHMS term in the UM;
  - cold TZ: ×8-14 (to ×463);
  - LM: ×355 (harzburgite), ×154 (20 % mix), ×12 (basalt), with the largest ×577.
- **The gain is zero** above ~1700 K and above 60 GPa, so the term never alters the hot mantle (Fig. S.14: solid and dashed lines coincide above each dehydration step). - **Nodes that lose capacity with the DHMS term.** 14 nodes of the mixes lose capacity (none in harzburgite or basalt): 12 in the TZ at 13-20 GPa, 2 in the LM at 27-28 GPa, all at 1467-1633 K on the hot edge of the DHMS field. There a trace of DHMS (< 0.14 wt%, mostly < 0.01) triggers the rebuild (S.10 step 6), and the MSH part's NAM water is recomputed on K&O phases:
  - forsterite instead of sb24's wadsleyite at 13.4 GPa and 1467 K (K&O is Fe-free, so its olivine-wadsleyite boundary lies deeper), capacity ×0.18-0.30 at p = 0.5-0.75;
  - hyWad instead of ringwoodite at 19-20 GPa (×0.73-0.98 at p = 0.05-0.65);
  - the Al-free bridgmanite law instead of the Al-bearing one at 27-28 GPa (×0.12-0.19 at p = 0.65).
  A fix is proposed in S.18.
- **Range:** no node is NaN or below the 10⁻¹⁰ floor; the minimum is 4.3·10⁻⁴ wt% (UM, P ≈ 0).

**Composition dependence (Fig. S.19).** Capacity is not linear in p: it steps at silica saturation (S.15) and wherever the MSH pair changes (S.9), and in the UM wherever a hydrous phase needs the Ca and Al of basalt or the Mg of peridotite. A test minimised the capacity at p = 0-1 in steps of 0.05 (40 × 40 nodes per sector) and compared it with a linear read between table compositions at the held-out p:

| read | UM: > ×2 off | TZ: > ×2 off | LM: > ×2 off | LM median / p95 error (wt%) |
|---|---|---|---|---|
| 3 compositions (0, 0.3, 1), the former table | 24 % | 8.7 % | 31 % | 0.020 / 0.27 |
| the ten table compositions | 6.0 % | 1.9 % | 2.1 % | 0.004 / 0.021 |

The UM columns use the production database rule (`pmbe` = 0.65). The largest single UM error of the three-composition read is 7.2 wt%: the serpentinite and lawsonite water was spread linearly across the whole mix range.

## S.12 Capacity along geotherms

<img src="md_figs/S13_geotherms.png" alt="S.13" width="1100">

*Figure S.13. Capacity along three illustrative geotherms (bottom), with DHMS (solid) and without (dashed).*

The geotherms are schematic, chosen to cross every regime:

```julia
adiabat:    T = 1600 + 12P (P ≤ 25),  1900 + 600(P − 25)/110                                         (12.1)   # Mantle adiabat, ~1900 K at 25 GPa, 2500 K at 135 GPa.
warm slab:  T = 700 + 60P (≤ 10),  1300 + 15(P − 10) (≤ 25),  1525 + 6(P − 25)                        (12.2)   # Warm slab top.
cold slab:  T = 500 + 45P (≤ 10),  950 + 15(P − 10) (≤ 25),  1175 + 6(P − 25)                         (12.3)   # Cold slab top.
```

*Symbols.* T is in K and P in GPa; each path is piecewise linear with the breaks shown. Capacities along a path are bilinear interpolations on the table nodes of the sector containing P.

Along the three paths (harzburgite / 20 % mix / basalt, wt%):

| P (GPa) | mantle adiabat | warm slab | cold slab |
|---|---|---|---|
| 3 | 0.022 / 0.032 / 0.014 (opx, omph) | 7.1 / 9.6 / 3.5 (atg, law) | 12.4 / 10.0 / 5.8 (atg, law) |
| 6.5 | 0.032 / 0.025 / 0.017 | 0.10 / 0.08 / 0.61 | 5.6 / 4.3 / 3.9 (phA, law) |
| 12 | 0.087 / 0.063 / 0.021 | 0.51 / 0.39 / 0.018 | 10.8 / 8.4 / 0.013 (phE, gtmj) |
| 18 | 0.67 / 0.56 / 0.027 (wa) | 1.2 / 0.94 / 0.026 (ri) | 7.9 / 6.8 / 0.37 (phD) |
| 23 | 0.99 / 0.75 / 0.027 (ri) | 4.3 / 3.4 / 0.027 (shB) | 7.9 / 6.4 / 0.020 (phD) |
| 30 | 0.032 / 0.089 / 0.022 | 1.4 / 0.98 / 0.047 (shB) | 8.2 / 8.0 / 4.1 (phD) |
| 45 | 0.032 / 0.089 / 0.011 | 0.032 / 0.089 / 0.022 | 11.3 / 11.0 / 4.4 (phH / phD) |
| 58-62 | 0.032 / 0.089 / 0.006 | 0.032 / 0.089 / 0.006 | 0.032 / 0.088-0.089 / 0.47 (pv, st) |
| 90 | 0.032 / 0.089 / 0.44 (st) | 0.032 / 0.089 / 0.45 | 0.032 / 0.089 / 0.46 |
| 125 | 0.012 / 0.037 / 0.28 (ppv, apbo) | 0.011 / 0.038 / 0.29 | 0.010 / 0.043 / 0.30 |

Reading along the paths:
- **Adiabat.** The transition-zone maximum (1 wt% in harzburgite ringwoodite) drops by a factor of 30 into bridgmanite at 25 GPa. That drop is the classic dehydration-melting contrast at the top of the lower mantle.
- **Warm slab.** It loses its serpentinite and lawsonite water by 6.5 GPa. It regains DHMS only in the TZ (shB) and loses them again by ~33 GPa.
- **Cold slab.** It keeps 5-12 wt% in the peridotite blocks from the surface to ~52 GPa. Phase H then dehydrates on this path at ~1340 K, just below the 60 GPa cut. Basalt carries lawsonite to 9.5 GPa, is nearly dry to ~16 GPa, and carries phase D, then phase H, from ~25 GPa.

## S.13 Seams and hand-overs

<img src="md_figs/S15_seams.png" alt="S.15" width="1100">

*Figure S.15. Capacity against P through 0-30 GPa at four temperatures; vertical lines at the 7 and 25 GPa seams; grey: each block's hand-over window.*

- **7 GPa seam** (UM last row against TZ first row, nodes whose T agree to 9 K). The capacity ratio has median 1.000 in every composition, 5-95 % range 0.94-1.01. The largest excursion is 0.12 decades, in the 0.95 mix (0.08 in basalt), on the lawsonite-out front, where the two sectors' T grids fall on opposite sides of the reaction.
- **25 GPa seam** (same T grid on both sides). The ratio is 1.000 at every node of every composition.
- **Hand-over windows** (Fig. S.15):
  - ume compositions, 5-6 GPa: the capacity moves linearly from ume's value to the sb24 + DHMS value, which carries phase A from 6 GPa at cold T;
  - basalt, 9.5-10 GPa: the lawsonite water (3.9 wt% at 770 K) is released over the window, and sb24 basalt continues nearly dry (majoritic garnet + cpx) until phase D forms from ~13 GPa below ~1100-1200 K.
- **Other steps are reactions, not seams:**
  - phase A-out at 900-960 K (6-10 GPa);
  - the phase D cut at 50 GPa (one T node of notch in the dehydration edge);
  - the DHMS cut at 60 GPa (from ~11 wt% in the peridotites, or ~5 wt% in basalt, to NAM values, below ~1300 K);
  - the stishovite switch in basalt.

## S.14 The oxybarometer and native iron tables

<img src="md_figs/S16_fO2_metal.png" alt="S.16" width="1100">

*Figure S.16. log10 fO2 of harzburgite and basalt at R ≈ 0.03 per sector (top two rows) and fO2 (solid) and native iron (dashed) against R at three nodes for four compositions (bottom).*

This pass reproduces the 2026-09-19 reference generator (`watercarb/Perple_X_MAGEMin/MAGEMin/s_fO2_dens_metal.jl`) operation for operation.

**Grid and compositions.** sb24 is minimised on the (P, T) grid of every sector, for 36 values of R and ten compositions. The R axis is 30 asinh-spaced nodes to the knee R = 0.4 followed by 6 logit-spaced nodes to 0.95:

```julia
R_i = R₀ sinh(u_i),  u ∈ linspace(asinh(−1.96/R₀), asinh(0.4/R₀), 30),  R₀ = 2·10⁻³               (14.1)   # Asinh part of the signed ferric axis (Rasinhspace).
R_i = 1/(1 + e^{−v_i}),  v ∈ linspace(logit 0.4, logit 0.95, 7)[2:end]                                (14.2)   # Logit tail toward the stoichiometric limit (Rlogspace).
```

*Symbols.*
- R_i: the ith node of the signed ferric axis.
- R₀: the asinh scale (the axis is near-linear for |R| ≪ R₀ and logarithmic beyond).
- u_i, v_i: the uniformly spaced transformed coordinates.
- logit x = ln(x/(1 − x)).
- linspace(a, b, n): n equally spaced values from a to b.

The oxybarometer is log-linear in R/(1 − R). Asinh spacing reaches 0.01 log units at 30 nodes, which no uniform grid does (TOex_math S.4).

**Post-processing.**
- fO2 is taken as `out.fO2`.
- Native iron is the summed weight fraction of feg, fea and fee (`metal_wt_frac`), clipped at 0.
- Non-finite values are filled from the mean of their finite P-T neighbours, sweep by sweep within each (R, composition, sector) slice, as the shipped tables were. In production this fills ten ±Inf nodes (TZ, R = 0.4, mixes 0.75 and 0.85).
- Along R, every finite fO2 column is forced non-decreasing by isotonic regression (pool-adjacent-violators), with a 10⁻⁶ relative tilt so it is strictly increasing (`enforce_R_monotone!`). In production 180 932 of the 300 000 columns need it, mostly for sub-10⁻³ dips.
- Non-finite nodes are reported per sector and R.

**Build.** The tables need the corrected local sb24 build: Fe site species, SLB property modifier and volume interaction. `~/.julia/artifacts/Overrides.toml` points `MAGEMin_jll` at it, so the "from MAGEMin_jll" message is expected. The stock build is ~2 log units too oxidised above ~50 GPa.

**Order and chunking.** Points are minimised in chunks of 200 000 in column-major (i, j, r, k, q) order. The chunk only splits the work; each point is independent.

Written layout:
- header `nP nT nR nEM` plus the eight mix fractions;
- the three sector lines;
- the R list;
- then, per R node and sector, the blocks basalt, harzburgite and mixes, as nP × nT;
- 4 decimals (fO2) or 6 (metal).

Code: `:fO2` block of `sfO2.jl`; `Rhybridspace`, `enforce_R_monotone!`, `metal_wt_frac`, `fmtv`.

## S.15 The silica-saturation step table

<img src="md_figs/S17_step.png" alt="S.17" width="1100">

*Figure S.17. The two edges of the silica-saturation step at R = 0.01: the last basalt fraction with an olivine-family phase or ferropericlase (pa) and the first with free silica (pb). Grey: no bracket in [0.3, 1].*

Capacity and fO2 are step-shaped in basalt fraction p at silica saturation, while StagYY reads the compositions linearly. `StagfO2step.dat` gives the reader the two edges of the step. Each node is on a coarse grid (40 P × 16 T per sector) for each of the ten R values of `step_R`: −1, −0.2, 0, 0.01, 0.03, 0.08, 0.15, 0.3, 0.5 and 0.8.

```julia
pa = sup{p ∈ [0.3, 1] : an ol/wa/ri/mw phase is present},   pb = inf{p ∈ [0.3, 1] : coe/st/qtz/seif/apbo present}   (15.1)   # The two step edges.
```

*Symbols.*
- p: the basalt mass fraction of the bulk p X_B + (1 − p) X_H.
- pa: the largest p at which a phase of the olivine family (ol, wa, ri) or ferropericlase (mw) is still present.
- pb: the smallest p at which a free silica phase (coesite, stishovite, quartz, seifertite or α-PbO2) is present.

Both are evaluated at the node's (P, T, R).

Each edge is found by bisection with 10 halvings (resolution 7·10⁻⁴):
- Only nodes whose end members (p = 0.3 and 1) bracket a change are bisected.
- The other nodes are NaN.
- Both searches share the two end-member minimisations.

Code: `:step` block of `sfO2.jl`.

## S.16 Hydrous density of solid and melt

<img src="md_figs/S18_density.png" alt="S.18" width="1100">

*Figure S.18. Top: percent density decrease of harzburgite solid, basalt solid and harzburgite melt at 3.9 wt% H2O. Bottom left: the melt water partial molar volume of the ideal model (used) and of the Jing & Karato hard-sphere law. Bottom right: all four corrections against water at 20 GPa, 2000 K.*

`StagSdens.dat` and `StagMdens.dat` hold the percent density decrease of hydrous solid and melt against the dry material at the same (P, T). The grid is 50 P (10⁻⁴-50 GPa) × 50 T (500-3000 K) × 50 H2O (0-8 wt%), with two slots: basalt and harzburgite. Their use in StagYY is TOex_math S.17.

**Solid.** Ideal volume mixing of the anhydrous rock with the Gerya et al. (2004) water volume, condensed branch:

```julia
V = (1 − x) M_A/ρ_A + x V_H2O^Gerya,   Δρ_S = 100 [1 − (M_H/V)/ρ_A]                                  (16.1)   # ρ_A: sb24 anhydrous density, x: water mole fraction, M_H the hydrous molar mass.
```

*Symbols.*
- w: the water in wt%.
- x: the water mole fraction, x = (w/M_H2O)/(w/M_H2O + (100 − w)/M_A).
- M_A: the mean molar mass of the anhydrous oxide mixture (g/mol per oxide mole).
- M_H: the hydrous mean molar mass, M_H = (1 − x)M_A + x M_H2O.
- ρ_A: the anhydrous sb24 density (g/cm³, MAGEMin `rho`, interpolated on a 50 × 50 grid).
- V: the hydrous molar volume (cm³/mol).
- V_H2O^Gerya: the water volume 10 ∂G_s/∂P (J/bar to cm³/mol), with G_s the Gerya et al. (2004) water Gibbs energy at X_liq = 1 (TOex_math 17.3).
- Δρ_S: the percent density decrease.

**Melt.** The dry melt density comes from the Jing & Karato (2011) hard-sphere law (TOex_math 17.4-17.5). Water enters by ideal mixing, with a partial molar volume from first-principles molecular dynamics:

```julia
V̄_H2O(P, T) = V₁(T) · V_BM3(P; V₀ = 25 cm³, K₀ = 1.91 GPa, K′ = 5.80)/25                            (16.2)   # FPMD compression curve, scaled to the 1-bar Jing & Karato water volume V1(T).
V₁(T) = [M_H/ρ_JK(1 bar, X_hyd) − (1 − x) M_A/ρ_JK(1 bar, X_dry)]/x                                  (16.3)   # 1-bar partial volume of water in the JK melt.
Δρ_M = 100 [1 − M_H/((1 − x) M_A/ρ_dry + x V̄_H2O)/ρ_dry],  floored at 0                               (16.4)   # Melt density decrease, melt_H2O_model = :ideal.
```

*Symbols.*
- V̄_H2O: the partial molar volume of water in the melt (cm³/mol).
- V₁(T): its 1-bar value from the Jing & Karato law.
- V_BM3(P; V₀, K₀, K′): the third-order Birch-Murnaghan volume, solving P = (3/2)K₀[(V₀/V)^{7/3} − (V₀/V)^{5/3}][1 + (3/4)(K′ − 4)((V₀/V)^{2/3} − 1)] by bisection.
- ρ_JK(P, X): the Jing & Karato density of composition X; X_dry and X_hyd are the dry and hydrous melt.
- ρ_dry: the dry melt density at (P, T).
- Δρ_M: the percent density decrease of the melt.

**Why `:ideal`.** The alternative `:JK` evaluates the hydrous melt with the same hard-sphere law. Its water sphere (ξ_H2O = −0.35) collapses under pressure:
- V̄_H2O falls to ~4 cm³/mol at 5 GPa and below zero above ~25 GPa (Fig. S.18);
- hydrous melt therefore becomes denser than dry melt from ~5 GPa, and the floor sets 89 % of the table to zero.

The FPMD volumes (Karki, Ghosh & Karato 2021: 12.4, 10.5 and 8.2 cm³/mol at 10, 20 and 50 GPa) and the sink-float estimates (Sakamaki et al. 2006, Matsukage et al. 2005) keep water light. `:ideal` gives 9.7 cm³/mol at 20 GPa and 7.8 at 40 GPa.

The shipped `StagMdens.dat` was made with `:ideal` and is reproduced byte for byte; the reference script's settings said `:JK`.

Code: `PT_H2O_ρ`, `Gerya_solid_interpolator`, `Gerya_solid_H2O_density_correction`, `Jing_Karato_silicate_melt_density`, `melt_H2O_partial_volume`, `bm3_volume`.

## S.17 Numerics and cost

**Optimisation round** (no number changed):
- NAM laws evaluated on a threaded node list, and fugacities shared between laws with the same clamp (opx, cpx pairs).
- The 10.8 M-tuple index of the fO2 pass (~430 MB) replaced by `CartesianIndices` slices in the same column-major order.
- The step bisection restricted to bracketed nodes, and the two end-member minimisations shared by both searches.

MAGEMin itself dominates the cost: each point is independent, so thread count and chunk size do not change results.

**Cost** of the production run (24-core Apple silicon, 16 threads, 128 GB; peak resident memory ~20 GB, from the 200 000-point fO2 chunks):

| pass | points | wall time |
|---|---|---|
| NAM fits + DHMS setup (phase H iteration) | - | 3.4 s |
| UM (ume/mbe, 10 compositions) | 100 000 | 393 s * |
| UM hand-over (sb24 + DHMS) | 14 500 | not timed separately |
| TZ + LM (sb24 + DHMS) + mbe below 10 GPa | 208 500 | 445 s * |
| fO2 + metal | 10 800 000 | 13 079 s (3.6 h, ~830 points/s) |
| step table | 19 200 nodes × bisection | 704 s |
| density tables | 2 × 2 500 sb24 + 250 000 melt | 10 s |
| total | | 3 h 55 min |

\* Ten-composition rerun of 2026-09-27 on 12 threads, alongside the two variant runs on 6 threads each; the other rows are the full production run on 16 threads. The sH2O-only passes (`SFO2_DO=sH2O_um,sH2O_deep`) take ~15 min this way.

## S.18 Choices and open items

**Decisions taken (with the reasons given above):**

1. **Water models.** PS for the NAM laws, BW for the DHMS, Gerya for densities (S.2).
2. **Melt off in the UM databases.** Hydrous phases above the solidus are handled by StagYY's melting routines (S.3).
3. **Laws held flat outside their boxes.** The boxes cover the calibration constraints (S.4). The opx-Al term is held at 1073 K below its lowest run temperature (Mierdel et al. 2007).
4. **Stishovite switch kept as a step** at the extrapolated CaCl2 boundary (S.4).
5. **Transition-zone laws measured against melt rather than free fluid** are accepted: at those temperatures melt solubility is high and little fluid is present.
6. **Bridgmanite law from Al2O3 only** (`Al_pv`, S.5); no depth term; no Fe³⁺ factor.
7. **DHMS as a point-wise capacity** on the sb24 assemblage (S.10), MSH only.
8. **Al phases out of scope;** majorite tested and not adopted (S.10).
9. **Phase H on Nishi et al. (2014),** with its effective volume. Phase H dehydration temperature is a dataset choice against Ohtani et al. (2014) (S.8).
10. **Basalt on mbe to 9.5-10 GPa** (lawsonite limit), keeping the sector layout (S.3).
11. **`:ideal` melt water volume** (S.16).
12. **Post-perovskite water divided by the Townsend D** (4.4). The code previously multiplied, giving Al-free post-perovskite 75-525 ppm instead of 2-12 ppm on the 30 ppm reference. The fix is from 2026-09-27; its reference follows `Al_pv`.

   Effect on the production maps (only 4583 lower-mantle nodes that hold post-perovskite changed):
   - harzburgite ppv median 0.013 → 0.010 wt%;
   - the p = 0.3 mix 0.016 → 0.032 wt% (its ppv carries ~6 wt% Al2O3, so the Al-bearing ratio and the peak of the bridgmanite law apply).

13. **DHMS cut at 60 GPa,** the upper limit of MgSiO4H2 phase H (Ohtani et al. 2014; Tsuchiya & Umemoto 2019). The one-layer release of cold-slab water there is the dissociation into bridgmanite + ice VII (S.8).

14. **StagH2O at the ten compositions of the fO2 tables** (2026-09-27), in every sector, with the fO2 block order and a header that names them (S.1). The three-composition read was more than ×2 off at 24-31 % of the UM and LM compositions between blocks (S.11). No water axis: capacity at saturation does not depend on bulk water (S.1).
15. **UM database switch at `pmbe` = 0.65,** where ume and mbe agree best; the mixes from 0.65 follow mbe into the TZ to 9.5-10 GPa, as basalt does (S.3).

**Open:**

1. **Nodes that lose capacity with the DHMS term** (14 in the mixes, S.11). Proposed: keep the larger of C_sb24 and the rebuilt C in (10.5). NAM water is not in the K&O energies, so a trace DHMS cannot tell whether it beats a water-rich wadsleyite or Al-bearing bridgmanite. The larger value also removes the drop at the edge of the DHMS field.
2. **Production byte test.** The fO2, metal, step and density tables must still be compared with a fresh full run of the reference generator. That run needs `melt_H2O_model = :ideal`, and the Python splice of the shipped files means bytes are compared against a fresh reference run, not against the shipped files.
3. **Porting.** The generator goes into the StagPP and watercarb Julia copies at green light, and the old `src/sf_H2O.jl` (which writes the former layout) is retired. The Fortran reader of this layout is in watercarb since 2026-09-27 (`load_s_fO2_data`, `thermochem_int('solubility')`): header `nP nT nEM` + mixes, no water axis, linear in basalt fraction between neighbouring blocks, sector bounds checked against the fO2 table. `solubility_Xmid` is removed from the namelist (2026-09-27), and the par files that set it were cleaned.

## S.19 Constants and sources

**Constants:**
- R = 8.31446 J/mol/K;
- molar masses as in `MM` (MgO 40.30, SiO2 60.08, FeO 71.85, Fe 55.85, H2O 18.015);
- T₀ = 298.15 K;
- NAM constants of S.4; BDG anchors of S.5; activity anchors of S.7; phase H of S.8;
- hand-over windows, limits and grids of S.1.

**Sources by use.** Volume and page numbers or DOIs were re-checked for Litasov (2003), Okamoto & Maruyama (1999), Townsend (2016) and K&O; the others are as cited in the code and earlier notes.


- **Thermodynamic engine:**
  - Riel, Kaus, Green et al. (2022) MAGEMin, G³ 23, https://doi.org/10.1029/2022GC010427;
  - Stixrude & Lithgow-Bertelloni (2024) sb24;
  - Evans & Frost (2021) J. Petrol. (ultramafic set);
  - Green et al. (2016) J. Metamorph. Geol. 34:845 (metabasite set).
- **Water:**
  - Pitzer & Sterner (1994) J. Chem. Phys. 101:3111, https://doi.org/10.1063/1.467624;
  - Brodholt & Wood (1993) J. Geophys. Res. 98:519;
  - Robie, Hemingway & Fisher (1978) USGS Bull. 1452;
  - Gerya et al. (2004) Phys. Chem. Miner. 31:429.
- **NAM laws:**
  - Dong, Fischer, Stixrude & Lithgow-Bertelloni (2021) AGU Advances 2, e2020AV000323;
  - Mierdel et al. (2007) Science 315:364, https://doi.org/10.1126/science.1135422;
  - Rauch & Keppler (2002) Contrib. Mineral. Petrol. 143:525;
  - Bromiley et al. (2004) Am. Mineral. 89:941;
  - Bromiley & Keppler (2004) Contrib. Mineral. Petrol. 147:189;
  - Lu & Keppler (1997) Contrib. Mineral. Petrol. 129:35;
  - Bolfan-Casanova, Keppler & Rubie (2000) EPSL 182:209, https://doi.org/10.1016/S0012-821X(00)00244-2;
  - Umemoto et al. (2006) (CaCl2-type transition);
  - Townsend, Tsuchiya, Bina & Jacobsen (2016) EPSL 454:20, https://doi.org/10.1016/j.epsl.2016.08.009;
  - Guimond et al. (2023) MNRAS (compilation of calibration ranges);
  - "Chen et al. 2021" (stishovite, as labelled in `sf_H2O.jl`; not independently located).
- **Bridgmanite:**
  - Litasov et al. (2003) EPSL 211:189, https://doi.org/10.1016/S0012-821X(03)00200-0;
  - Fu et al. (2019) GRL;
  - Liu et al. (2021) EPSL 570.
- **DHMS:**
  - Komabayashi & Omori (2006) PEPI 156:89, https://doi.org/10.1016/j.pepi.2006.02.002;
  - Komabayashi, Omori & Maruyama (2005) PEPI 153:191, https://doi.org/10.1016/j.pepi.2005.07.001;
  - Fabrichnaya (2004) (dry reactions, via K&O);
  - Ohtani et al. (2003) (phase D + pv + pc).
- **Phase H:**
  - Nishi et al. (2014) Nat. Geosci. 7:224, https://doi.org/10.1038/ngeo2074;
  - Ohtani et al. (2014) GRL, https://doi.org/10.1002/2014GL061690;
  - Bindi et al. (2014) Am. Mineral. (phase H cell);
  - Tsuchiya (2013) GRL (ab initio phase H);
  - Tsuchiya & Umemoto (2019) GRL 46, https://doi.org/10.1029/2019GL083472 (phase H → bridgmanite + ice VII).
- **Basalt hand-over:** Okamoto & Maruyama (1999) Am. Mineral. 84:362, https://doi.org/10.2138/am-1999-0320.
- **Melt density:**
  - Jing & Karato (2011) GCA 75:6780, (2012) GCA 85:357;
  - Karki, Ghosh & Karato (2021) Sci. Rep. 11:10588;
  - Sakamaki, Suzuki & Ohtani (2006) Nature 439:192;
  - Matsukage, Jing & Karato (2005) Nature 438:488.
