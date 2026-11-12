# TOex numerics

How the excess-oxygen solve is computed, what it costs and how reliably it converges. The physics and
the equations are in `TOex_math.md` (sections S.n below); the comparison of the model with published
data is in `TOex_bench.md`. Everything here concerns the Fortran implementation in watercarb
(`volatiles_mod.f90`, `thermochem_mod.f90`) unless the Julia reference (`src/TOex.jl`) is named.

The numbers were produced on 2026-09-19 from 40-step StagYY continuations of the local run's
1.64 Gyr frame (512 × 64 cells, 8 ranks, core cooling on, ten-composition tables), each build run
twice, and from the Julia drivers `md_figs/scripts/bench_numerics.jl` (sections N1, N2) and
`md_figs/scripts/bench.jl` (B8, B10, B14, B15). The log extracts the figures are made from are in
`md_figs/scripts/numerics_logs/`.

| # | What is measured | Result |
|---|---|---|
| 1 | Cost of the melting and solubility section per time step | 1.75 s → 0.59 s (−66 %); whole step 2.29 s → 1.10 s |
| 2 | Convergence in the run | 0 to 2 cells of 32 768 unconverged per solve; 6.4 outer iterations per cell; oxygen closure ≤ 4·10⁻⁸ of the tolerance scale |
| 3 | Run-to-run reproducibility | horizontal-mean fO₂ differs by 0.04 to 0.17 log units between two runs of the same binary after 40 steps |
| 4 | Table read against a live minimisation | 98 % of 4000 random points within 0.3 log units, median 0.03 |
| 5 | Julia solver sweep, 1296 cases | 88.9 % converged; the rest is one infeasible corner; carbon closure 3·10⁻¹⁶ |
| 6 | Fortran against Julia | 53 melt-only cases agree to 5·10⁻⁷ log units |
| 7 | Hydrogen floor and closure test (2026-11-08) | worst C or H closure 47.7 → 5·10⁻⁹ in a 12-step continuation; failure counts unchanged |

## 1. Structure of the solve

One cell is one scalar equation: the total oxygen demand of all sinks at a trial log fO₂ equals the
cell's budget (S.3). Evaluating that demand is itself a nest of one-dimensional problems:

```
outer root in log fO₂              bisection, then Illinois, inside a warm-start window   (S.12)
 ├ solid iron demand                direct read of the tabulated oxybarometer              (S.4, section 2.4)
 ├ native iron at the root          read off the same columns, once per solve               (section 2.4)
 ├ melt iron demand                 Illinois inversion of the melt oxybarometer            (S.5)
 └ C-O-H state                      non-ideality loop on the fugacity coefficients         (S.11)
     └ hydrogen root                Illinois on the fluid hydrogen (or the scale λ)        (S.10)
         └ carbon-activity root     Illinois on aC with its cap                            (S.7)
             ├ speciation           closed-form quadratic                                  (S.9)
             └ fluid molar volume   Newton on the equation of state                        (section 2.2)
```

A call-stack sample of one rank shows where the time went before the 2026-09-19 changes and after:

| share of the TOex work | before | after |
|---|---|---|
| building the solid oxybarometer for the cell | 66 % | 11 % (two fO₂ and two native-iron columns, no inversion) |
| outer root with all demands | 15 % | 51 % |
| – of which the C-O-H state (hydrogen, carbon activity, speciation, fluid volume) | 11 % | 15 % |
| everything else | 19 % | 38 %: MPI reductions (waiting for the slowest rank) 18 %, gather, scatter and counters 20 %; the native-iron read of the converged state, 13 % while it went through `thermochem_int`, no longer shows |

## 2. The tools

### 2.1 Bracketed roots: bisection, regula falsi, Illinois

Every root whose function has kinks (caps, phase switches, a table) is solved with a bracket that is
never given up. Regula falsi replaces one endpoint by the zero of the secant; on a convex function
the same endpoint survives every step and convergence degrades to linear. The Illinois variant halves
the stored function value of an endpoint that has survived twice, which restores superlinear
convergence (order about 1.44) at one function evaluation per step ([Dowell & Jarratt 1971](https://doi.org/10.1007/BF01934364));
S.26 has the figure. It is used for the outer fO₂ root after a few bisections, the melt-iron
inversion, the hydrogen root and the carbon-activity root. A derivative obtained at the cost of a
second function evaluation would not beat it per evaluation (Newton's order 2 over two evaluations
is 1.41 per evaluation), which is why these roots were not converted to Newton.

### 2.2 Newton where the function is smooth and the derivative is free

Two inner problems are smooth with a closed-form derivative, and both now use Newton with a
bracketed fallback:

- **Fluid molar volume.** The [Zhang & Duan (2009)](https://doi.org/10.1016/j.gca.2009.01.021)
  equation of state is a polynomial-plus-exponential in density, so dP/dρ is analytic and shares
  every term with P. The iteration starts from the last root found in the cell, runs only if the old
  warm bracket (0.7 to 1.4 of that volume) straddles the root, stops on a relative step of 10⁻⁸
  (the error is then its square), and is accepted by the same residual test as before
  (|F| ≤ 10⁻¹⁰ P). A step that leaves the bracket, a non-positive dP/dρ or a failed acceptance hands
  the problem to the Illinois solver and its ladder of brackets. Residual evaluations per cell fell
  from 1245 to 294 to 308, and root failures per output interval from about 2250 to about 1750:
  the bracketed solver stopped on a volume tolerance and on stiff, liquid-like roots that could miss
  the residual test, which Newton does not. The Julia reference (`zd_findroot`) has no warm volume and
  runs the same Newton inside each bracket of its ladder: same roots as its old 200-step bisection to
  10⁻¹³ relative on 20 000 random states, 2 µs against 26 µs per root, benches B4 and B10 unchanged.
- **Column inversion near the silica-saturation step** (section 2.5): inside one segment between
  table nodes R = R₀ sinh(u) with u linear in fO₂, so dR/df = √(R₀² + R²) du/df needs no extra
  transcendental.

### 2.3 Warm starts

The outer root is searched first within ±0.25 log units of the cell's previous fO₂, then ±2, then on
the full bracket [−120, +40] (S.12). In the test runs 98.6 % of cells are solved in the narrow
window, 1.3 % in the wide one and 0.1 % on the full bracket. Where the residual has several roots
the window also acts as the selection rule: the root nearest last step's is reported. The fluid
volume root is warm-started from the last volume found in the same cell (hit rate 99.9 %).

### 2.4 The tabulated solid oxybarometer and how it is read

`StagfO2.dat` holds log fO₂ on 100 P × 100 T × 36 R nodes per pressure sector for ten compositions
(basalt fraction p = 0, 0.05, 0.1, 0.2, 0.5, 0.65, 0.75, 0.85, 0.95, 1); `StagMetal.dat` holds the
native iron on the same grid. The R axis is signed (ferric above the FeO reference, metal below) and
uniform in u = asinh(R/0.002), which resolves both the dilute end and the metal branch with 36
nodes. Reads are bilinear in (P, T) and linear in u; below the first node a column sits on its
Fe–FeO plateau, above the last it continues linearly in logit(R).

**Composition blend.** A cell at basalt fraction p uses the two compositions bracketing p. They are
not averaged in log fO₂: the two rocks share one fO₂, each holds the ferric iron its own column
gives at that fO₂, and the bulk ratio is the sum weighted by mass fraction and FeO content,

```
R_bulk(f) = [a₁ R₁(f) + a₂ R₂(f)] / (a₁ + a₂),    aᵢ = wᵢ · FeOᵢ
```

The solver's unknown is f, so the solid demand is this expression evaluated directly: each column
is inverted by a binary search and a linear step in u, and the loading follows from R_bulk. Until
2026-09-19 the reader first inverted the expression to tabulate f at the 36 bulk R nodes (36 root
finds per cell, 60 halvings each) and the solver then inverted that column again at every trial f.
Removing the double inversion takes out what the one-sweep column had left of that cost (section
3.1), and it is also more faithful: the 36-node column lost the kinks of the two block columns that fall between bulk nodes
(figure N2c: median 2·10⁻⁴ log units, 90th percentile 0.03, 99th 0.33, worst 1.1). The run's
horizontal means did not move beyond their repeat scatter. The native iron of the converged state
is read the same way: the reader hands over the two blocks' iron columns with their mass weights,
and each block's iron is taken at the ratio that block holds at the root fO₂. Through
`thermochem_int` the same number cost a second (P, T) bracket and a 60-step bisection per cell; the
two agree to 3·10⁻¹⁵ in weight fraction over a full run.

### 2.5 The step-aware read

Along p at fixed (P, T, R) the fO₂ is smooth except where olivine-family phases or ferropericlase
run out (p = pa) and a silica phase appears (p = pb): a jump of 1.7 to 3 log units within 0.05 in p,
at p ≈ 0.77 to 0.80 below 25 GPa and 0.66 to 0.72 in the lower mantle. A blend across the bracket
that contains the jump smears it over the whole bracket. `StagfO2step.dat` tabulates pa and pb on
40 × 16 × 10 nodes per sector (bisection in p on phase presence, inside the table generator
`solve_sH2O_fO2`, switch `fO2_step`). Where a cell's
bracket straddles [pa, pb] and the eight surrounding corners agree on one clean step (spread ≤ 0.06,
window ≤ 0.10), each side is read, or linearly extrapolated, from its own two nodes; inside the
window the read follows the olivine side below 25 GPa and ramps in the lower mantle. Anywhere else,
and whenever the file is absent, the read is the plain blend. Cells whose bracket can reach the step
keep a 36-node column (the step rule is applied per R node); its blended nodes are now found in one
sweep: the bulk ferric oxygen is evaluated once at the merged nodes of the two columns, which
brackets every target, and Newton (section 2.2) finishes inside that segment. The column differs
from the old per-node bisection by at most 1.8·10⁻¹¹ log units over a full run. A column assembled
from nodes that disagree on the step can lose its monotony in R; it is restored strictly (each node
at least 10⁻⁶ above the previous), because a tie makes the inversion non-unique. Curved
interpolation across p (shape-preserving cubic, spline) was tested against a fine scan and changes
the worst case by about 10 % and the mean by nothing; more nodes and the step rule are what help.

### 2.6 Column cache

The fO₂ and native-iron columns depend on (P, T, p) only, and are asked for twice per step for about half the cells
(`toex_solve_cells`, then `toex_react_cell` during water redistribution). A direct-mapped cache of
16 384 entries keyed on the exact three reals returns the second request without touching the
table; one third of all requests hit it. About 19 MB per rank.

### 2.7 Tolerances

The outer residual is accepted at

```
atol = max( rel · max(|TOₑₓ|, 10⁻⁶ capacity),  10⁻¹² (Mc + Mh + capacity) )
```

so a cell with a vanishing budget is judged against a floor set by what its sinks could hold, and
round-off in the carbon and hydrogen terms cannot masquerade as a residual (S.3). Inner tolerances
are chosen two to four orders tighter than the level above them: fugacity coefficients 10⁻⁸,
fluid volume 10⁻¹³ relative (Illinois) or a Newton step of 10⁻⁸, column Newton step 10⁻¹¹ log units.

### 2.8 Smooth switches

Every change of regime that a bracketed root has to cross is a ramp, not a switch, because a jump
in a sink puts a step in the residual and a step can hide the root: the graphite–carbonate host is a
logistic of width 0.05 log units around EMOG(D) (S.6); melt CO₂ and melt H₂ approach their ceilings
through a soft cap (S.8, S.10); the carbonate solidus ledge is a smoothstep in pressure and a tanh
in log concentration (S.16).

### 2.9 Failure handling

No sign change on the full bracket, a collapsed bracket with a large residual, or a non-finite state
is reported as a failure and never as a success. The cell keeps its previous partition, which is
exactly conservative under the incremental write-back, the failure is counted with its cause, and
the run continues (S.12).

The unified engine also requires the carbon and hydrogen balances of the final state to close to
`toex_closure_tol` = 10⁻⁶ (S.11). A cell whose solid holds more water than the budget even at the
lowest tabulated activity (the hydrogen floor, S.15) is closed onto the solid when it has no melt and
failed when it has melt. Each step the log prints the three counts:
`TOex: H floor (no melt, solid holds all H) N   H floor with melt (failed) N   closure failures N`.

## 3. Benchmarks

### 3.1 Cost

| build | melting and solubility, s per step | whole step, s | fluid EOS evaluations per cell | oxybarometer reads per cell | fluid EOS root failures per interval |
|---|---|---|---|---|---|
| per-node bisection (until 2026-09-19) | 1.754 (repeat 1.742) | 2.29 | 1245 | 205 | 2287 |
| one-sweep column and cache | 1.001 | 1.70 | 1243 | 205 | 2226 |
| plus Newton on the column and the fluid EOS | 0.756 | 1.27 | 294 | 205 | 1656 |
| plus direct blend | 0.615 (repeat 0.618) | 1.13 | 308 | 31 | 1782 |
| plus native iron read off the columns (shipped) | 0.588 (repeat 0.588) | 1.10 | 308 | 31 | 1699 |

The PETSc direct solve of the same steps takes 0.06 to 0.09 s, so the section is still the largest
item of a time step (53 %, from 77 %). Fluid reaction steps (285 to 298 per cell) and
carbon-activity solves (139 to 145) are unchanged: the chemistry does the same work, the table
handling around it does less. The ten-composition tables and the step-aware read cost nothing
measurable against the four-composition tables (1.848 against 1.861 s per step on the old build).

![N1](md_figs/N01_cost_per_build.png)

*Figure N1. (a) Wall time per step of the melting and solubility section and of the PETSc direct solve for the five builds. (b) Work per cell per solve.*

### 3.2 Convergence and reproducibility in the run

| build | unconverged cells per solve, mean (max) | outer iterations, mean (worst cell) | worst oxygen closure | carbon shares, worst \|Σ − 1\| |
|---|---|---|---|---|
| old | 0.43 (1); repeat 0.18 (2) | 6.32 (211; 217) | 1.1·10⁻⁸; 3.1·10⁻⁸ | 1.3·10⁻¹³ |
| shipped | 0.39 (1); repeat 0.53 (2) | 6.42 (206; 123) | 3.6·10⁻⁸; 3.7·10⁻⁸ | 5.1·10⁻¹⁴ |

The unconverged cells are melt-free, cold and shallow (325 K, 0.4 GPa, a budget of 0.002 wt% O in
the logged example); they keep their previous partition. Two runs of the same binary do not reproduce each
other bit for bit (parallel direct solver, tracer ordering): after 40 steps the horizontal-mean fO₂
profiles differ by 0.04 to 0.17 log units and the ferric ratio by about 10⁻³, the native-iron
profile not at all. Old against the direct-blend build differ by 0.05 to 0.07 and that build against
the shipped one by 0.08 to 0.09, inside that scatter, so profile
comparisons cannot separate the builds and the exactness statements of section 2 rest on in-code
checks instead. On the 0.13 Gyr frame of the same run the step from four to ten compositions moved
the profile by 0.13 (repeat scatter 0.09) and left 9 cells per solve unconverged against 3.

![N2](md_figs/N02_convergence_and_blend.png)

*Figure N2. (a) Cells left unconverged per solve and (b) outer iterations, old and shipped builds with their repeats. (c) Cumulative distribution of the difference between the 36-node column read and the exact blend at 1500 random points between R nodes.*

### 3.3 The table read against live minimisations

4000 random mantle-like points (P 1 to 130 GPa, T = 1600 + 7 P ± 200 K, Fe³⁺/ΣFe log-uniform 0.005
to 0.3, p uniform 0 to 1; `md_figs/scripts/randtest/`), |read − live sb24 minimisation| in log units:

| table and reader | median | P90 | P99 | share above 0.3 | share above 1 | same, p > 0.5 only |
|---|---|---|---|---|---|---|
| four compositions, blend | 0.100 | 0.68 | 1.95 | 22.1 % | 5.8 % | 39.9 % / 10.9 % |
| seven compositions, blend | 0.047 | 0.21 | 1.06 | 6.8 % | 1.1 % | 9.7 % / 1.6 % |
| ten compositions, blend | 0.031 | 0.17 | 0.98 | 5.2 % | 1.0 % | 8.2 % / 1.4 % |
| ten compositions, step-aware (shipped) | 0.026 | 0.14 | 0.58 | 2.0 % | 0.4 % | 1.8 % / 0.1 % |

For Fe³⁺/ΣFe ≤ 0.1, where the model lives, 0.9 % of points are above 0.3, and none between p = 0.2
and 0.65. At the nodes the table is faithful to its source: worst composition-node difference 0.32
log units at 130 GPa, where the R nodes are sparsest; at seven off-node fractions from 0.1 to 0.88
the worst read is 0.19; a pyrolite cell (p = 0.2) reads within 0.02 of the bulk minimisation on the
xenolith path and within 0.11 to 130 GPa. One corner remains: p < 0.1 at 12 to 17 GPa with
Fe³⁺/ΣFe above about 0.08, where an accessory ferric phase (akimotoite, calcium ferrite) switches on
in nearly pure harzburgite and fO₂ drops by up to 5 log units between p = 0.02 and 0.04 (worst read
3.9). It is a second step of the same kind; the same machinery would take it. What a mixed cell
should read is a modelling choice as well: a cell at p = 0.8 is a mechanical mixture of basalt and
harzburgite tracers, and the table reproduces the homogenised bulk.

![B15](md_figs/B15_step_aware_read.png)

*Figure N3. Table reads against live minimisations across the silica-saturation step, basalt fraction 0.5 to 1 at four conditions: live sb24 (black), the oxygen-weighted blend (orange) and the step-aware read (green dashed); dotted lines are the table nodes.*

### 3.4 Julia solver sweep

`partition_Oₑₓ` with the live solid branch (25 R nodes from MAGEMin per case) over P {1, 3, 10, 25}
GPa × T {1300, 1600, 1900, 2200} K × ϕ {10⁻³, 10⁻², 0.1} × TOex {10⁻⁴, 10⁻³, 10⁻²} × TC {10⁻⁶,
10⁻², 0.1} × H {0, 10⁻⁴, 10⁻²} (mass fractions):

| | value |
|---|---|
| converged | 1152 of 1296 (88.9 %) |
| not converged | 144, all at TOex = 1 wt% with TC = 10⁻⁶: 9 per (P, T) |
| oxygen closure, converged cases | worst \|Σ sinks − TOex\| = 0.99 atol |
| carbon closure, converged cases | worst 3·10⁻¹⁶ |
| wall time | 63 s, 0.05 s per case including the MAGEMin nodes |

The 144 failures are one infeasible corner: 1 wt% excess oxygen in a carbon-free cell exceeds the
ferric capacity of solid plus melt (about 0.9 wt% O for pyrolite), so no fO₂ satisfies the budget
and the solver reports it as such. With carbon present the carbonate sink absorbs the same budget
and every case converges.

### 3.5 Fortran against Julia

53 melt-only cases (2026-09-13) agree to 5·10⁻⁷ log units on the IW and FMQ offsets and to the
printed digits on the partition, so the scientific benchmarks of `TOex_bench.md`, which run the Julia
code, describe what the simulation computes.

### 3.6 Hydrogen floor and the closure test (2026-11-08)

The newcode2 run logged a worst joint residual of 3.3 and later 47.7. That is the C and H closure
`res_ideality` (S.11), printed for rank 0's worst cell only. The cause was the hydrogen floor of
S.15 in cold shallow solid. The fix (eqs. 15.2, 15.3 and the closure test) was checked in Julia
against the previous code on the same cells, then in Fortran:

| test set | cells | result |
|---|---|---|
| frame-600 floor cells | 7 | closure 47.7 → 0, solid share of H 48.7 → 1, fO₂ unchanged |
| cold cells with melt (frame 600) and the frame-530 dump | 81 | identical |
| floor cell with a melt trace, ϕ = 10⁻⁹, 10⁻⁶, 10⁻³ | 3 | failed, flagged floor with melt |
| random near-surface cells | 400 | 391 identical, 9 floor cells closed, closure above 10⁻⁶: 9 → 0 |
| regression grid, all sectors | 756 | 624 identical, 66 floor closed, 66 floor with melt failed, no other difference |
| legacy engine | 48 | identical |
| Fortran harness on the 7 cells | 7 | as Julia, a_solid = 10⁻³ |

The grid's three cells that failed before (20 GPa, 350 K, p = 1, no fO₂ root) still fail. In a
12-step continuation of frame 600 (`mpirun -np 16`) the worst closure stayed below 5·10⁻⁹ in every
step (3.3 and 47.7 before), the floor fired in one cell per step from step 6, never with melt,
no cell failed the closure test, the failed counts (29 to 35 per step) matched the old build step
for step and the water totals did not change.

## 4. Not done

- The hydrogen and carbon-activity roots restart from their full brackets at every outer iteration
  (139 to 145 carbon-activity solves per cell); reusing the previous roots as a bracket is the next
  saving of size.
- Cells whose inputs have not changed since the last solve are solved again.
- Near the step the reader interpolates all ten composition columns in (P, T); four are needed.
- The second step, in depleted harzburgite (section 3.3), is not tabulated.

## Sources

- Dowell, M. & Jarratt, P. (1971). A modified regula falsi method for computing the root of an equation. *BIT* 11, 168–174. [doi:10.1007/BF01934364](https://doi.org/10.1007/BF01934364)
- Riel, N., Kaus, B. J. P., Green, E. C. R. & Berlie, N. (2022). MAGEMin, an efficient Gibbs energy minimizer: application to igneous systems. *Geochemistry, Geophysics, Geosystems* 23, e2022GC010427. [doi:10.1029/2022GC010427](https://doi.org/10.1029/2022GC010427)
- Stixrude, L. & Lithgow-Bertelloni, C. (2024). Thermodynamics of mantle minerals – III: the role of iron. *Geophysical Journal International* 237, 1699–1733. [doi:10.1093/gji/ggae126](https://doi.org/10.1093/gji/ggae126)
- Zhang, C. & Duan, Z. (2009). A model for C–O–H fluid in the Earth's mantle. *Geochimica et Cosmochimica Acta* 73, 2089–2102. [doi:10.1016/j.gca.2009.01.021](https://doi.org/10.1016/j.gca.2009.01.021)
