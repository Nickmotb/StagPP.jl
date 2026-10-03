# TOex benchmarks against published data

Only quantities the model computes itself are benchmarked here: the shared oxygen fugacity of a
cell, the ferric ratio of its solid and melt, where its carbon sits (carbonate, graphite/diamond,
melt CO₂, fluid), the fluid speciation, and the native iron. Each is run through the same code the
simulation uses (`partition_Oₑₓ`, `coh_state`, `solid_fO2_interpolant`) and compared with a
published measurement or an independent model. How the solve is computed, what it costs, how the
tables are read and how well the Fortran follows the Julia is in `TOex_numerics.md`; checks of the
input laws against their own sources (water EOS, JANAF constants, the reference buffers, the melt
solubility and oxybarometer fits, the solidus brackets) are re-runnable from
`md_figs/scripts/bench_inputs.jl` and recorded in `TOex_bench.txt`.

The thermodynamic dataset is sb24 ([Stixrude & Lithgow-Bertelloni 2024](https://doi.org/10.1093/gji/ggae126)) through MAGEMin ([Riel et al.
2022](https://doi.org/10.1029/2022GC010427)), in the fixed build `~/Desktop/MAGEMin_v2.0.3_sb24`; the simulation reads it through the
ten-composition tables of 2026-09-19. Run data are the 1.64 Gyr frame of the local run
(`~/Desktop/localop/1.64 Gyr`, core cooling on). The driver is `md_figs/scripts/bench.jl`; each
section writes its numbers to `md_figs/scripts/bench_results_B<n>.txt` and its figure to
`md_figs/B<nn>_*.png`; `ONLY="B11 B12"` re-runs a subset and `RUN=<stem>` points sections 6 and 7 at
another frame. Verdicts: PASS inside the stated tolerance; OFFSET a systematic difference that is
understood and carried as a known bias; FAIL a disagreement not understood or not yet corrected.

| # | Model output tested | Against | Offset (model − published) | Verdict |
|---|---|---|---|---|
| 1 | Graphite-saturated C-O-H fluid speciation at IW | [Matveev et al. (1997)](https://doi.org/10.1016/S0016-7037%2897%2900142-7), GFluid; Sokol et al. (2004); [Hirschmann et al. (2012)](https://doi.org/10.1016/j.epsl.2012.06.031) | fO₂ of the match 0.12 log units; species mean \|Δx\| 0.017 | PASS |
| 2 | Cratonic lithosphere: fO₂ and carbonate → diamond depth from the full solve | [Stagno et al. (2013)](https://doi.org/10.1038/nature11679), [Woodland & Koch (2003)](https://doi.org/10.1016/S0012-821X%2803%2900379-0) | crossover −35 to 0 km; ΔFMQ +1.3 to +1.9 (mean +1.5) | PASS on the crossover, OFFSET on fO₂ |
| 3 | MORB: ferric ratio and fO₂ of the melt from a pyrolite source | [Cottrell & Kelley (2011)](https://doi.org/10.1016/j.epsl.2011.03.014), [Zhang et al. (2018)](https://doi.org/10.1016/j.chemgeo.2018.01.006), [Berry et al. (2018)](https://doi.org/10.1016/j.epsl.2017.11.032) | Fe³⁺/ΣFe +0.02 on the mean of the three; ΔFMQ −0.1; D(Fe³⁺) +2 | PASS; OFFSET on the partition |
| 4 | Upper-mantle self-reduction at fixed Fe³⁺/ΣFe, bulk minimisation and as the run reads it | [Frost & McCammon (2008)](https://doi.org/10.1146/annurev.earth.36.031207.124322), [Stagno et al. (2013)](https://doi.org/10.1038/nature11679) | slope 0.0; ΔFMQ +0.6 to +1.4 (mean +1.1), the same for the run's table read | OFFSET: the ferric garnet component of sb24 |
| 5 | Native-iron onset depth and amount | [Rohrbach et al. (2007)](https://doi.org/10.1038/nature06183), [Frost et al. (2004)](https://doi.org/10.1038/nature02413) | onset +50 to +150 km; amount −0.35 wt% | OFFSET |
| 6 | Lower-mantle redox at fixed Fe³⁺/ΣFe | HeFESTo ([SLB 2024](https://doi.org/10.1093/gji/ggae126)), [Frost & McCammon (2008)](https://doi.org/10.1146/annurev.earth.36.031207.124322) | ΔIW −0.7 vs HeFESTo, 0.0 vs F&M; metal −0.4 to −1 wt% | PASS to 120 GPa at Fe³⁺/ΣFe ≤ 0.05 |
| 7 | The local run at 1.64 Gyr against the data | xenolith, MORB-source and lower-mantle fO₂ | upper mantle 0 to +1.5 log units with Fe³⁺/ΣFe −0.025; lower mantle 0 above 1000 km, +2 to +5 below 2000 km | fO₂ PASS to 1500 km; FAIL on the inventory: ferric iron moved from the top to the base |

## 1. The C-O-H fluid

The fluid branch solves speciation and carbon activity for a given fO₂, budget and host (S.9,
S.11); its output is tested against the one direct measurement of a graphite-saturated fluid at
mantle pressure. [Matveev et al. (1997)](https://doi.org/10.1016/S0016-7037%2897%2900142-7) measured XH₂O at 1273.15 K and 2.4 GPa, buffered at IW. The
model, run through the production ideality loop (`coh_state`), reproduces the measured XH₂O at
log fO₂ = −13.975, 0.12 log units below the [Frost (1991)](https://doi.org/10.1515/9781501508684-004) IW the experiment was run on and 0.28
below the sb24 IW.

| quantity | model | published | offset |
|---|---|---|---|
| log fO₂ at which XH₂O = 0.1172 | −13.975 | IW ([Frost 1991](https://doi.org/10.1515/9781501508684-004)) −13.85; IW (sb24) −13.70 | −0.12; −0.28 |
| XCH₄ | 0.846 | 0.814 measured; 0.815 GFluid | +0.032 |
| XH₂ | 0.036 | 0.054 measured; 0.036 GFluid | −0.018; 0.000 |
| XCO₂ | 0.00004 | – | – |
| mean \|Δx\| over the measured species | | | 0.017 |

XCH₄ is high by 0.03 because C₂H₆ is excluded by design (GFluid, Zhang & Duan 2010, carries 0.032 of
it). Against the Sokol et al. (2004) conditions reproduced with GFluid (1693 K, 5.7 GPa, matched XO)
XH₂O is 0.8366 against 0.8478, the same C₂H₆ offset, and the O-H fluid at IW, IW−1 and IW−2 gives
X_H₂ 52.6 / 77.9 / 91.8 % against the 48 / 74 / 90 % of [Hirschmann et al. (2012)](https://doi.org/10.1016/j.epsl.2012.06.031), Fig. 7.

The redox sequence at 1500 K, 3 GPa with 0.5 wt% C and 0.1 wt% H₂O shows the emergent structure
the solve relies on: a methane-dominated fluid 4 log units below EMOG(D), graphite and carbonate
at 40.7 % each exactly on the buffer, carbonate a carbon-limited sink above it.

![B4](md_figs/B04_coh_speciation.png)

*Figure B1. The redox sequence across EMOG(D): (a) fluid speciation and carbon activity, (b) where the carbon sits.*

## 2. Cratonic lithosphere: redox and carbon speciation from the full solve

The closest thing to a field test of the whole partition: cratonic garnet peridotites give fO₂
against depth ([Woodland & Koch 2003](https://doi.org/10.1016/S0012-821X%2803%2900379-0); [Stagno et al. 2013](https://doi.org/10.1038/nature11679)), and the same xenolith suites carry the
carbonate-to-diamond transition that [Stagno et al. (2013)](https://doi.org/10.1038/nature11679) place at about 150 to 200 km. A pyrolite
cell (p = 0.2, no melt, no water) is run through `partition_Oₑₓ` along a 40 mW/m² geotherm
(T = 673 + 125 P). Everything in the tables is an output; the inputs are the source, and the source
is the weak point of this test, so it is stated first.

**How the source is set.** The xenoliths do not give the two numbers the solve needs. Their
Mössbauer ferric ratios are the state *after* the carbon has equilibrated (Stagno's garnets carry
more Fe³⁺ the deeper they come from, because the reduction of carbonate put it there), and their
carbon is not measured on the same samples. So the source is not sample-by-sample: the silicate
ferric ratio R is taken from the whole-rock range of spinel peridotites that have not seen the
diamond field (0.02 to 0.04, [Canil & O'Neill 1996](https://doi.org/10.1093/petrology/37.3.609)) and the carbon from the range quoted for
cratonic peridotite (30 to 100 ppm). Both levers feed the same excess-oxygen budget: R through
half an oxygen per Fe³⁺ (270 ppm O at R = 0.03), the carbon through two oxygens per carbonate
carbon (80 ppm O at 30 ppm C, 270 at 100 ppm). Whether the carbon starts as carbonate or as
graphite is therefore a third input of the same size, and the table brackets it. That is why the
carbonate curves in figure B2(b) fan out: the spread between them is the uncertainty of the
source, not of the solve.

| source | crossover to majority diamond (model) | published | offset | ΔFMQ at 3 / 4.5 / 6.5 GPa (model) | xenolith band centre | offset |
|---|---|---|---|---|---|---|
| R 0.03, 30 ppm C as carbonate | 5.5 GPa, 141 km | 150–200 km ([Stagno et al. 2013](https://doi.org/10.1038/nature11679)) | −9 to −59 km | +0.28 / −1.06 / −2.15 | −1.0 / −2.5 / −4.0 | +1.3 / +1.4 / +1.9 |
| R 0.03, 100 ppm C as carbonate | 6.0 GPa, 153 km | same | 0 to −47 km | +0.28 / −0.93 / −1.66 | same | +1.3 / +1.6 / +2.3 |
| R 0.02, 30 ppm C as carbonate | 4.5 GPa, 115 km | same | −35 to −85 km | −0.41 / −1.48 / −2.66 | same | +0.6 / +1.0 / +1.3 |
| R 0.03, 30 ppm C as graphite | 4.5 GPa, 115 km | same | −35 to −85 km | −0.32 / −1.43 / −2.59 | same | +0.7 / +1.1 / +1.4 |
| mean | | | −20 to −70 km | | | +1.5 (range +0.6 to +2.3) |

Three things come out that are not put in:

- **The carbonate-to-diamond transition lands where the xenoliths put it**, 115 to 153 km against
  150 to 200 km, for every source. The transition is gradual (30 ppm of carbonate goes from 90 % to
  10 % over 50 km) because the host is ramped over `wEMOG` and the ferric budget pays for it as it goes.
- **Carbon buys its own reduction.** As carbonate reduces to diamond it hands two oxygens per carbon
  to the silicate, so the solid's ferric ratio rises with depth (0.030 to 0.039 at 30 ppm, 0.030 to
  0.058 at 100 ppm) while fO₂ falls: the [Stagno et al. (2013)](https://doi.org/10.1038/nature11679) observation reproduced by mass balance.
- **The fO₂ curve is 1.3 to 1.9 log units above the xenoliths** for the R = 0.03 sources, 0.6 to 1.4
  for the two more reduced ones, and outside the ±1 log unit tolerance band of figure B2(a) below
  5 GPa. The offset is the ferric garnet component of sb24 (section 4).

**Would calibrating the ramp help?** No. The ramp (`wEMOG`, S.6) sets how sharply carbon switches
host around EMOG(D); it moves the width of the crossover in figure B2(b), not the fO₂ at a given
ferric ratio, which is what is off in B2(a). Sharpening it would only make the 50 km transition
steeper. The crossover depth already agrees because it is fixed by where the geotherm's fO₂ crosses
EMOG(D), and both are shifted together by the garnet offset, so the error cancels there and shows
only in the absolute fO₂. The fix is a ferric garnet component with the right energy (skiagite, or
Fe³⁺ in majorite; section 4), which is a dataset change, not a parameter.

![B11](md_figs/B11_cratonic_carbon.png)

*Figure B2. (a) ΔFMQ along the cratonic geotherm for four sources, the EMOG(D) buffer on the same path, the garnet-peridotite band (read from [Stagno et al. 2013](https://doi.org/10.1038/nature11679) and [Woodland & Koch 2003](https://doi.org/10.1016/S0012-821X%2803%2900379-0)) and a ±1 log unit tolerance around it. (b) The fraction of the carbon held as carbonate against depth; grey: where [Stagno et al. (2013)](https://doi.org/10.1038/nature11679) place the transition.*

## 3. MORB: the ferric ratio of the melt from a pyrolite source

MORB glasses are the best-measured ferric ratio on Earth, and the source is constrained
independently: spinel peridotite whole rocks carry Fe³⁺/ΣFe 0.03 to 0.04 ([Canil & O'Neill 1996](https://doi.org/10.1093/petrology/37.3.609)).
The composition is the same pyrolite as every other section (p = 0.2: 0.2 basalt + 0.8 harzburgite
in the two-endmember mixing of S.4). A cell at 1.5 GPa and 1573 K with 100 ppm H₂O and 30 ppm C is melted to ϕ = 0.05,
0.10, 0.15 with the source ferric ratio set through the excess-oxygen budget. The melt's ferric
ratio, its fO₂ and the ferric partition between melt and residue are outputs.

| quantity (source 0.03–0.04, ϕ = 0.10, C as carbonate) | model | published | offset |
|---|---|---|---|
| melt Fe³⁺/ΣFe | 0.133–0.166 | 0.16 ± 0.01 ([Cottrell & Kelley 2011](https://doi.org/10.1016/j.epsl.2011.03.014)); 0.14 ± 0.01 ([Zhang et al. 2018](https://doi.org/10.1016/j.chemgeo.2018.01.006)); 0.10 ± 0.02 ([Berry et al. 2018](https://doi.org/10.1016/j.epsl.2017.11.032)) | −0.01 / +0.01 / +0.05; +0.02 on the mean of the three |
| melt ΔFMQ | −0.37 to +0.21 | 0.0 ± 0.2 ([Cottrell & Kelley 2011](https://doi.org/10.1016/j.epsl.2011.03.014)) | −0.1 |
| D(Fe³⁺) melt/residue | 5.9–6.6 | 3–5 ([Sorbadere et al. 2018](https://doi.org/10.1016/j.gca.2018.07.019); [Davis & Cottrell 2018](https://doi.org/10.2138/am-2018-6280)) | +2 |
| residue Fe³⁺/ΣFe | 0.020–0.028 | abyssal peridotite 0.02–0.03 | 0.00 |
| same source, C as graphite: melt Fe³⁺/ΣFe, ΔFMQ | 0.101–0.137, −1.08 to −0.29 | as above | −0.03; −0.7 |

The full sweep (ϕ 0.05 to 0.15, R 0.02 to 0.05, both carbon states) is in `bench_results_B12.txt`.
For a source at 0.03 to 0.04 and 10 % melting the melt comes out at 0.13 to 0.17 and ΔFMQ −0.4 to
+0.2, which spans the three published MORB values and the measured fO₂. The ferric partition is
5.9 to 6.6, above the 3 to 5 of the experiments: the melt oxybarometer ([Hirschmann 2022](https://doi.org/10.1016/j.gca.2022.04.005), basalt
composition) takes more ferric iron per log unit than the sb24 solid gives up, so the same source
produces a melt about 0.02 richer in Fe³⁺ than the experiments would. Carried as OFFSET; the sign
is such that the run's melts are marginally too oxidised for their source.

The dashed curve shows the same source with its carbon as graphite instead of carbonate: the melt
then dissolves the carbon as CO₂ and pays the two oxygens per carbon out of the ferric budget, so
the melt is 0.7 log units more reduced and 0.03 lower in ferric ratio. Thirty ppm of carbon is
enough because carbonate oxygen and ferric oxygen are the same size at these concentrations (80 and
270 ppm O). With 100 ppm of graphite the melt drops to ΔFMQ −3. The initial speciation of the
source carbon is therefore part of the source definition, and in the run it is whatever the
previous step left in the cell.

![B12](md_figs/B12_morb_ferric.png)

*Figure B3. (a) Melt ferric ratio against source ferric ratio for three melt fractions with the carbon as carbonate, and for ϕ = 0.10 with the carbon as graphite (dashed), against the three published MORB values; green: the source range. (b) The melt's ΔFMQ against the MORB range.*

## 4. Upper-mantle self-reduction, and what the run reads

Pyrolite (p = 0.2) with the ferric ratio held at 0.03 and T = 1550 + 12 P, the solid fO₂ from the
live sb24 oxybarometer (`solid_fO2_interpolant`) and, next to it, what a p = 0.2 cell of the
simulation reads from its table at the same points. The published picture ([Frost & McCammon 2008](https://doi.org/10.1146/annurev.earth.36.031207.124322)) is
a decrease of 0.6 to 1.0 log units per GPa from ΔFMQ ≈ 0 at 1 GPa, with cratonic garnet peridotites
([Stagno et al. 2013](https://doi.org/10.1038/nature11679)) at ΔFMQ −1 near 3 GPa, −2.5 near 4.5 and −4 near 6.5.

| quantity | model, bulk minimisation | model, as the run reads it | published | offset (bulk / run) |
|---|---|---|---|---|
| ΔFMQ at 3 GPa | −0.38 | −0.37 | −1.0 | +0.6 / +0.6 |
| ΔFMQ at 4.5 GPa | −1.36 | −1.35 | −2.5 | +1.1 / +1.2 |
| ΔFMQ at 6.5 GPa | −2.56 | −2.55 | −4.0 | +1.4 / +1.5 |
| slope 1–8 GPa (log units per GPa) | −0.60 | −0.60 | −0.6 to −1.0 | 0.0 / 0.0 |
| mean ΔFMQ offset | | | | +1.1 / +1.1 |

**The bulk minimisation** has the right slope and sits 0.6 to 1.4 log units above the xenolith
points, and the run reads the same curve: its ten-composition table with the oxygen-weighted,
step-aware read is within 0.02 log units of the bulk minimisation on this path (the earlier
two-composition table read 1.7 to 2.6 above it; how the table is read and how well is in
`TOex_numerics.md`, section 3.3). The offset is a property of the dataset. sb24's only ferric garnet
endmember is andradite, but calcium is not what limits it in a fertile rock: the model garnet at
6.5 GPa keeps 7 % grossular next to 2.8 % andradite and carries Fe³⁺/ΣFe ≈ 0.12, about what cratonic
garnets of that depth carry, and still reads 1.4 log units high. The [Holland et al. (2018)](https://doi.org/10.1093/petrology/egy048) igneous
set, whose garnet is also andradite-only, gives the same or more (ΔFMQ −0.31 / −0.84 / −1.96 at
3 / 4.5 / 6.5 GPa). Closing 1.4 log units at 1600 K takes about −23 kJ per mole of andradite, several
times a plausible uncertainty on that endmember, which points to a missing ferric garnet component
(skiagite, ferric majorite) whose energy is set by the skiagite–fayalite–ferrosilite equilibrium the
xenolith oxybarometer uses, rather than to a retuning. In harzburgite (0.9 wt% CaO) calcium does
limit the andradite, and the rock reads several log units more oxidised than pyrolite at the same
ferric ratio. Carried as OFFSET; untested until the endmember is added to the dataset.

![B13](md_figs/B13_metal_onset_blend.png)

*Figure B4. (a) Native iron at fixed net ferric ratio along T = 1600 + 7 P against the [Rohrbach et al. (2007)](https://doi.org/10.1038/nature06183) onset and the [Frost et al. (2004)](https://doi.org/10.1038/nature02413) lower-mantle amount (section 5). (b) Upper-mantle self-reduction on the xenolith path T = 1550 + 12 P at constant Fe³⁺/ΣFe = 0.03: bulk pyrolite minimisation and the table read of a p = 0.2 cell, against the [Frost & McCammon (2008)](https://doi.org/10.1146/annurev.earth.36.031207.124322) band and the [Stagno et al. (2013)](https://doi.org/10.1038/nature11679) points.*

## 5. Native iron: where it appears and how much

Along T = 1600 + 7 P with the net ferric ratio fixed, the metal the oxybarometer reports (the
disproportionation product, S.4) is an output. [Rohrbach et al. (2007)](https://doi.org/10.1038/nature06183) found Fe–Ni metal saturation
from about 250 km (8 GPa) in fertile peridotite; [Frost et al. (2004)](https://doi.org/10.1038/nature02413) estimate about 1 wt% metal in
the lower mantle.

| quantity | model (net R 0.02 / 0.03 / 0.04) | published | offset |
|---|---|---|---|
| first metal (> 0.001 wt%) | 12 / 16 / 16 GPa; 303 / 401 / 401 km | ~250 km, 8 GPa ([Rohrbach et al. 2007](https://doi.org/10.1038/nature06183)) | +50 to +150 km |
| metal at 26–30 GPa | 0.67 / 0.64 / 0.61 wt% | ~1 wt% ([Frost et al. 2004](https://doi.org/10.1038/nature02413)) | −0.33 to −0.39 wt% |
| ΔIW at 26–30 GPa | −1.24 / −1.23 / −1.23 | −1 to −1.5 ([Frost & McCammon 2008](https://doi.org/10.1146/annurev.earth.36.031207.124322)) | 0.0 |

Metal appears 50 to 150 km deeper than in the experiments, in traces (0.01 to 0.07 wt%) that come and
go through the transition zone, and reaches 0.6 to 0.7 wt% only where bridgmanite takes over. Two
reasons: the fO₂ on this path is still 0.1 to 1.2 log units above IW at 8 to 10 GPa (the garnet
offset of section 4 again), and sb24's majorite carries less ferric iron than the experimental
garnets, so it disproportionates less. The lower-mantle amount is 60 to 70 % of the published
estimate. Carried as OFFSET; both errors pull the same way (the model is late and light on metal).

## 6. The lower mantle, and the disagreement

The question this section answers: is carbonate near the CMB of the runs a property of the solid
oxybarometer, or of the oxygen the runs carry there? With the unfixed sb24 port MAGEMin burned its
metal by 60 to 80 GPa and climbed to ΔIW +2.3 at the CMB where HeFESTo stays at −0.5 with metal; the
fixed build removes that.

**At fixed net Fe³⁺/ΣFe, sb24 behaves like HeFESTo.** Pyrolite (p = 0.2) along T = 1600 + 7 P:

| P (GPa) | T (K) | ΔIW, R = 0.02 (model) | metal (model) | HeFESTo, [SLB 2024](https://doi.org/10.1093/gji/ggae126) (R 0.02) | [Frost & McCammon 2008](https://doi.org/10.1146/annurev.earth.36.031207.124322) | offset vs HeFESTo / F&M | R = 0.05 | R = 0.10 | R = 0.15 |
|---|---|---|---|---|---|---|---|---|---|
| 25 | 1775 | −1.20 | 0.59 % | −0.5 to 0, ~1 % metal | −1 to −1.5, ~1 % | −0.7 / 0.0; metal −0.4 % | −1.19 (0.50 %) | −1.17 (0.36 %) | −1.15 (0.21 %) |
| 60 | 2020 | −1.18 | 0.39 % | same | same | −0.7 / 0.0; −0.6 % | −1.17 (0.30 %) | −1.16 (0.15 %) | −1.11 (0.01 %) |
| 80 | 2160 | −1.22 | 0.22 % | same | same | −0.7 / 0.0; −0.8 % | −1.21 (0.13 %) | −1.04 (0) | −0.16 (0) |
| 100 | 2300 | −1.32 | 0.09 % | same | same | −0.8 / 0.0; −0.9 % | −1.24 (0) | +0.08 (0) | +0.96 (0) |
| 120 | 2440 | −1.47 | 0.01 % | same | same | −1.0 / 0.0; −1.0 % | −0.14 (0) | +1.24 (0) | +2.12 (0) |
| 130 | 2510 | −1.23 | 0 | same | same | −0.7 / 0.0; −1.0 % | +0.39 (0) | +1.63 (0) | +2.38 (0) |

At the ferric ratio of a pyrolite (0.02) the model holds ΔIW = −1.2 to −1.5 with metal present to
120 GPa: the metal-buffered lower mantle of [Frost & McCammon (2008)](https://doi.org/10.1146/annurev.earth.36.031207.124322) and of [SLB (2024)](https://doi.org/10.1093/gji/ggae126) Fig. 13. The
metal content is 0.2 to 0.6 wt%, half of the ~1 wt% quoted, and it fades toward the CMB where
HeFESTo keeps it; those are the residual differences. The oxidised branch appears only where the
rock carries more net ferric iron than disproportionation supplies: above ~100 GPa at R = 0.05,
above ~70 GPa at R = 0.15. That is the physical behaviour, since net ferric beyond what metal
balances has to be paid by fO₂.

**The run's deep mantle is inventory-driven.** The local run (core cooling on) at 1.64 Gyr, ΔIW
recomputed with the same buffer along its own horizontal-mean temperature profile:

| depth (km) | P (GPa) | T (K) | Fe³⁺/ΣFe | O_ex (wt% O) | ΔIW (model) | published | offset |
|---|---|---|---|---|---|---|---|
| 689 | 28 | 2377 | 0.144 | 0.12 | −1.45 | −1 to −1.5 with metal | 0.0 |
| 992 | 41 | 2449 | 0.157 | 0.14 | −1.20 | same | 0.0 |
| 1506 | 65 | 2540 | 0.152 | 0.13 | −0.62 | same | +0.4 |
| 2020 | 90 | 2346 | 0.197 | 0.19 | +1.42 | same | +2.4 |
| 2477 | 114 | 2394 | 0.215 | 0.26 | +2.85 | same | +3.9 |
| 2788 | 130 | 2563 | 0.218 | 0.38 | +3.84 | same | +4.8 |
| 2888 | 135 | 2955 | 0.202 | 0.39 | +2.75 | same | +3.8 |

The top 800 km of the lower mantle sits at ΔIW −1.2 to −1.5 with a reported Fe³⁺/ΣFe of 0.14 to 0.16:
the metal-bearing branch, on the data. From 1500 km down the rock carries a net O_ex of 0.13 rising
to 0.39 wt% O at the base, and at that budget the oxybarometer has no metal left to pay with, so ΔIW
climbs to +3 to +4. The budget, not the solver, sets this: at R = 0.15 to 0.20 the fixed-ratio
curves above leave the metal branch between 60 and 80 GPa, which is where the run does. That oxygen
is advected: subducted crust arrives with inherited ferric iron and carbonate and piles on the CMB
(basal O_ex was 0.009 wt% at the start and reached 0.25 wt% with the first downwelling). No
oxybarometer can reduce carbonate in rock that holds 0.4 wt% excess oxygen; the model is answering
the question it was asked.

So the disagreement splits in two. The dataset part is closed by the fixed build. The remaining part
is the oxygen budget of subducted material: how much ferric iron and carbonate the ingassing leg puts
into the crust, whether the surface boxes let it accumulate, and whether hot piles should hold it.
Those are model inputs, not thermodynamics, and section 7 shows the other end of the same transfer.

![B7](md_figs/B07_lower_mantle_redox.png)

*Figure B5. (a) ΔIW at fixed net ferric ratio against the HeFESTo band ([SLB 2024](https://doi.org/10.1093/gji/ggae126) Fig. 13) and the [Frost & McCammon (2008)](https://doi.org/10.1146/annurev.earth.36.031207.124322) bridgmanite range. (b) Native iron at the same points. (c) The local run at 1.64 Gyr below 660 km: ΔIW, ten times the reported ferric ratio, and the O_ex the rock carries.*

## 7. The run's upper mantle against the data

The same run above 660 km, horizontal means, against what sections 2 to 4 say the data require:

| depth (km) | P (GPa) | T (K) | Fe³⁺/ΣFe | O_ex (wt% O) | ΔFMQ (model) | ΔIW (model) | published | offset |
|---|---|---|---|---|---|---|---|---|
| 88 | 3 | 1816 | 0.008 | 0.014 | −0.92 | +2.68 | xenoliths ΔFMQ −0.5 to −1; MORB source 0 ± 0.5; source Fe³⁺/ΣFe 0.03 to 0.04 | 0.0 on fO₂; −0.025 on the ferric ratio |
| 192 | 8 | 2023 | 0.011 | 0.011 | −1.96 | +1.99 | xenoliths ΔFMQ −3 to −4 | +1.5 |
| 312 | 12 | 2079 | 0.015 | 0.015 | −4.06 | +0.28 | metal onset region, ΔIW ≈ 0 to +0.5 | 0.0 |
| 401 | 16 | 2147 | 0.027 | 0.021 | −5.78 | −1.28 | ΔIW −1 to 0 | −0.3 |

The fO₂ of the upper mantle is where the data put it at 90 km and from 300 km down, and 1.5 log units
high at 190 km, which is the dataset offset of section 4 and no more. What is off is the inventory
behind it: a silicate ferric ratio of 0.008 to 0.015 above 300 km against the 0.03 to 0.04 of the
MORB source. The run starts at Fe³⁺/ΣFe = 0.01 with 0.02 wt% carbon as graphite, so its source is
three to four times poorer in ferric iron than a pyrolite from the outset, and it stays there:
melting extracts ferric iron preferentially (section 3, D ≈ 6), everything erupted goes to the
surface boxes, and the upper-mantle O_ex, which fell from 0.0093 to 0.0052 wt% O in the first 10 Myr,
is back at 0.011 to 0.015 without the ferric ratio following. A melt of this source is more reduced
and poorer in Fe³⁺ than the MORB values of section 3, which were obtained from a source at 0.03 to
0.04. The oxygen that a pyrolite-like upper mantle would hold is of the size of what the basal
1400 km carry in section 6: the degassing and ingassing legs move it from the top of the mantle to
the bottom. Verdict FAIL for the run's inventory, not for the solver; which input is responsible
(the initial ferric ratio and carbon, the extraction, or the return leg) is for the parameter sweep
to separate. The frame was computed with the four-composition table; the ten-composition read moves
the horizontal means by about 0.1 log units (`TOex_numerics.md`, section 3.2).

## Sources

- Berry, A. J., Stewart, G. A., O'Neill, H. St C., Mallmann, G. & Mosselmans, J. F. W. (2018). A re-assessment of the oxidation state of iron in MORB glasses. *Earth and Planetary Science Letters* 483, 114–123. [doi:10.1016/j.epsl.2017.11.032](https://doi.org/10.1016/j.epsl.2017.11.032)
- Canil, D. & O'Neill, H. St C. (1996). Distribution of ferric iron in some upper-mantle assemblages. *Journal of Petrology* 37, 609–635. [doi:10.1093/petrology/37.3.609](https://doi.org/10.1093/petrology/37.3.609)
- Cottrell, E. & Kelley, K. A. (2011). The oxidation state of Fe in MORB glasses and the oxygen fugacity of the upper mantle. *Earth and Planetary Science Letters* 305, 270–282. [doi:10.1016/j.epsl.2011.03.014](https://doi.org/10.1016/j.epsl.2011.03.014)
- Davis, F. A. & Cottrell, E. (2018). Experimental investigation of basalt and peridotite oxybarometers: implications for spinel thermodynamic models and Fe³⁺ compatibility during generation of upper mantle melts. *American Mineralogist* 103, 1056–1067. [doi:10.2138/am-2018-6280](https://doi.org/10.2138/am-2018-6280)
- Frost, B. R. (1991). Introduction to oxygen fugacity and its petrologic importance. *Reviews in Mineralogy* 25, 1–9. [doi:10.1515/9781501508684-004](https://doi.org/10.1515/9781501508684-004)
- Frost, D. J. & McCammon, C. A. (2008). The redox state of Earth's mantle. *Annual Review of Earth and Planetary Sciences* 36, 389–420. [doi:10.1146/annurev.earth.36.031207.124322](https://doi.org/10.1146/annurev.earth.36.031207.124322)
- Frost, D. J., Liebske, C., Langenhorst, F., McCammon, C. A., Trønnes, R. G. & Rubie, D. C. (2004). Experimental evidence for the existence of iron-rich metal in the Earth's lower mantle. *Nature* 428, 409–412. [doi:10.1038/nature02413](https://doi.org/10.1038/nature02413)
- Hirschmann, M. M. (2022). Magma oceans, iron and chromium redox, and the origin of comparatively oxidized planetary mantles. *Geochimica et Cosmochimica Acta* 328, 221–241. [doi:10.1016/j.gca.2022.04.005](https://doi.org/10.1016/j.gca.2022.04.005)
- Hirschmann, M. M., Withers, A. C., Ardia, P. & Foley, N. T. (2012). Solubility of molecular hydrogen in silicate melts and consequences for volatile evolution of terrestrial planets. *Earth and Planetary Science Letters* 345–348, 38–48. [doi:10.1016/j.epsl.2012.06.031](https://doi.org/10.1016/j.epsl.2012.06.031)
- Holland, T. J. B., Green, E. C. R. & Powell, R. (2018). Melting of peridotites through to granites: a simple thermodynamic model in the system KNCFMASHTOCr. *Journal of Petrology* 59, 881–900. [doi:10.1093/petrology/egy048](https://doi.org/10.1093/petrology/egy048)
- Matveev, S., Ballhaus, C., Fricke, K., Truckenbrodt, J. & Ziegenbein, D. (1997). Volatiles in the Earth's mantle: I. Synthesis of CHO fluids at 1273 K and 2.4 GPa. *Geochimica et Cosmochimica Acta* 61, 3081–3088. [doi:10.1016/S0016-7037(97)00142-7](https://doi.org/10.1016/S0016-7037%2897%2900142-7)
- Riel, N., Kaus, B. J. P., Green, E. C. R. & Berlie, N. (2022). MAGEMin, an efficient Gibbs energy minimizer: application to igneous systems. *Geochemistry, Geophysics, Geosystems* 23, e2022GC010427. [doi:10.1029/2022GC010427](https://doi.org/10.1029/2022GC010427)
- Rohrbach, A., Ballhaus, C., Golla-Schindler, U., Ulmer, P., Kamenetsky, V. S. & Kuzmin, D. V. (2007). Metal saturation in the upper mantle. *Nature* 449, 456–458. [doi:10.1038/nature06183](https://doi.org/10.1038/nature06183)
- Rohrbach, A. & Schmidt, M. W. (2011). Redox freezing and melting in the Earth's deep mantle resulting from carbon–iron redox coupling. *Nature* 472, 209–212. [doi:10.1038/nature09899](https://doi.org/10.1038/nature09899) (figure B6 of the driver only)
- Sokol, A. G. et al. (2004), as reproduced with GFluid; the exact reference is not recorded in the notes and has no link here.
- Sorbadere, F., Laurenz, V., Frost, D. J., Wenz, M., Rosenthal, A., McCammon, C. & Rivard, C. (2018). The behaviour of ferric iron during partial melting of peridotite. *Geochimica et Cosmochimica Acta* 239, 235–254. [doi:10.1016/j.gca.2018.07.019](https://doi.org/10.1016/j.gca.2018.07.019)
- Stagno, V., Ojwang, D. O., McCammon, C. A. & Frost, D. J. (2013). The oxidation state of the mantle and the extraction of carbon from Earth's interior. *Nature* 493, 84–88. [doi:10.1038/nature11679](https://doi.org/10.1038/nature11679)
- Stixrude, L. & Lithgow-Bertelloni, C. (2024). Thermodynamics of mantle minerals – III: the role of iron. *Geophysical Journal International* 237, 1699–1733. [doi:10.1093/gji/ggae126](https://doi.org/10.1093/gji/ggae126)
- Woodland, A. B. & Koch, M. (2003). Variation in oxygen fugacity with depth in the upper mantle beneath the Kaapvaal craton, Southern Africa. *Earth and Planetary Science Letters* 214, 295–310. [doi:10.1016/S0012-821X(03)00379-0](https://doi.org/10.1016/S0012-821X%2803%2900379-0)
- Zhang, C. & Duan, Z. (2010). GFluid: an Excel spreadsheet for investigating C–O–H fluid composition under high temperatures and pressures. *Computers & Geosciences* 36, 569–572 (DOI not verified).
- Zhang, H. L., Cottrell, E., Solheid, P. A., Kelley, K. A. & Hirschmann, M. M. (2018). Determination of Fe³⁺/ΣFe of XANES basaltic glass standards by Mössbauer spectroscopy and its application to the oxidation state of iron in MORB. *Chemical Geology* 479, 166–175. [doi:10.1016/j.chemgeo.2018.01.006](https://doi.org/10.1016/j.chemgeo.2018.01.006)
