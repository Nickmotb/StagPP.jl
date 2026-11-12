# Benchmark slides of TOex_solver_flow.tex: (B) the model against the Stagno & Frost (2010) experiments run by run,
#  (C) the model's fO2 against natural peridotites and MORB glasses compiled by Cottrell et al. (2021, Table 3.1).
# Run: include after the Makie-free TOexLite loader, or julia --project=<StagPP.jl> bench_slides.jl

isdefined(Main, :TOexLite) || @eval using StagPP
const M = isdefined(Main, :TOexLite) ? Main.TOexLite : StagPP
using CairoMakie, MAGEMin_C, Printf, Statistics, DelimitedFiles

const OUT = joinpath(@__DIR__, "figs", "bench"); mkpath(OUT)

# ---------------- B: Stagno & Frost (2010), Tables 1-2 (data/stagno_frost_2010.csv) ----------------
raw = readdlm(joinpath(@__DIR__, "data", "stagno_frost_2010.csv"), ',', Any; header=true)[1]
runs = [(run=string(r[1]), comp=string(r[2]), TC=Float64(r[4]), P=Float64(r[5]), X=Float64(r[6]), f=Float64(r[7]),
         ph=split(string(r[8]))) for r in eachrow(raw)]
carb(r)  = any(p -> p in ("Mst", "Dol"), r.ph)                      # carbonate mineral observed
melt(r)  = "Liq" in r.ph                                            # melt observed
meltq(r) = "(Liq)" in r.ph                                          # melt inferred, hard to analyse
TK(r) = r.TC + 273.15
f_melt(r) = 5.44 - 21380/TK(r) + 0.078*(1e4r.P - 1)/TK(r) + log10(r.X)   # model melt-carbon law at aC = 1 (Eq. 14 of SF10)
f_emog(r) = M.logfO2_EMOG(TK(r), 1e4r.P)                                  # model EMOG(D), HP98

sub = [r for r in runs if carb(r) && !melt(r)]                      # subsolidus graphite/diamond + carbonate
liq = [r for r in runs if melt(r) || (meltq(r) && !carb(r))]        # melt-bearing
st(o, p) = (d = p .- o; (mean(d), sqrt(mean(d.^2)), length(d)))
ms, rs, ns = st([r.f for r in sub], f_emog.(sub)); ml, rl, nl = st([r.f for r in liq], f_melt.(liq))
@printf("B-1 subsolidus vs EMOG(D): n=%d mean %+.2f rms %.2f | melt runs vs melt law: n=%d mean %+.2f rms %.2f\n", ns, ms, rs, nl, ml, rl)

f = Figure(size=(1500, 760), fontsize=16)
a = Axis(f[1, 1], xlabel="measured log₁₀ fO₂ (Ir–Fe sensor)", ylabel="model log₁₀ fO₂",
         title="(a) fO₂ of graphite/diamond + carbonate, run by run", titlealign=:left)
lims = (-9.5, -2.0); lines!(a, [lims...], [lims...], color=:gray50, linestyle=:dash); limits!(a, lims, lims)
band!(a, [lims...], [lims[1]-0.5, lims[2]-0.5], [lims[1]+0.5, lims[2]+0.5], color=(:gray, 0.12))
scatter!(a, [r.f for r in liq], f_melt.(liq), color="#2B7BD6", markersize=12, strokewidth=0.5, strokecolor=:white,
         label=@sprintf("melt runs, melt-carbon law at the measured X_CO₂ (n=%d, mean %+.2f, rms %.2f)", nl, ml, rl))
scatter!(a, [r.f for r in sub], f_emog.(sub), color="#E8692C", marker=:diamond, markersize=15, strokewidth=0.5, strokecolor=:white,
         label=@sprintf("subsolidus carbonate runs, EMOG(D) from HP98 (n=%d, mean %+.2f, rms %.2f)", ns, ms, rs))
Legend(f[2, 1], a, framevisible=false, labelsize=13)

# (b) phases: carbonated solidus of the model against the observed melt
Tper(P) = P <= 3 ? 973.0 : P <= 5 ? 973 + 71.5(P - 3) : P <= 7 ? 1116 + 40.5(P - 5) : P <= 10 ? 1197 + 25.67(P - 7) : 1274 + 23.8(P - 10)
b = Axis(f[1, 2], xlabel="T (°C)", ylabel="P (GPa)", yreversed=true, title="(b) Melt and carbonate observed against the model's carbonated solidus", titlealign=:left)
Pl = range(2.0, 12.0, length=100); lines!(b, Tper.(Pl), collect(Pl), color=:black, linewidth=2.5, label="model carbonated-peridotite solidus")
jit(r) = r.P + 0.12*(r.comp == "EMOG1" ? 1 : -1)                     # separate the two capsules of one experiment
for r in runs
    obs = melt(r) ? :melt : meltq(r) ? :q : :none
    col = carb(r) ? "#E8692C" : "#2B7BD6"
    scatter!(b, [r.TC], [jit(r)], color=obs == :none ? (:white, 1.0) : obs == :q ? (col, 0.45) : col,
             strokecolor=col, strokewidth=2, marker=carb(r) ? :diamond : :circle, markersize=14)
end
agree = count(r -> (melt(r) || meltq(r)) == (r.TC > Tper(r.P)), runs)
@printf("B-2 melt observed vs model solidus: %d of %d runs agree\n", agree, length(runs))
elems = [MarkerElement(marker=:circle, color="#2B7BD6", strokecolor="#2B7BD6", strokewidth=2, markersize=14),
         MarkerElement(marker=:circle, color=(:white, 1.0), strokecolor="#2B7BD6", strokewidth=2, markersize=14),
         MarkerElement(marker=:circle, color=("#2B7BD6", 0.45), strokecolor="#2B7BD6", strokewidth=2, markersize=14),
         MarkerElement(marker=:diamond, color="#E8692C", strokecolor="#E8692C", strokewidth=2, markersize=14),
         LineElement(color=:black, linewidth=2.5)]
Legend(f[2, 2], elems, ["melt observed", "no melt", "melt inferred (bracketed)", "carbonate mineral observed",
       @sprintf("model solidus: melt predicted to its right (%d of %d runs agree)", agree, length(runs))],
       nbanks=2, framevisible=false, labelsize=13)
xlims!(b, 1050, 1650); ylims!(b, 11.8, 2.2)
Label(f[0, 1:2], "Benchmark B: Stagno & Frost (2010), graphite/diamond-saturated carbonated peridotite, 2.5–11 GPa, measured fO₂", font=:bold, fontsize=18)
rowsize!(f.layout, 1, Relative(0.72)); save(joinpath(OUT, "bench_B.png"), f, px_per_unit=2); println("wrote bench_B.png"); flush(stdout)

# ---------------- C: model fO₂ against natural rocks (Cottrell et al. 2021, Table 3.1) ----------------
data = Initialize_MAGEMin("sb24", verbose=false); CACHE = Dict{Any,Any}()
part(P, T, ϕ, TO) = M.partition_Oₑₓ(P, T, 0.2, ϕ, TO, 3e-5, 1e-4; verbose=false, data=data, out=true, mag_cache=CACHE)
TOs = 10 .^ range(-5, log10(1.5e-3), length=40)
# x: the source Fe³⁺/ΣFe, i.e. the solid's ratio for the same oxygen budget without melt (all Fe in the solid)
function curve(P, T, ϕ)
    R0 = [part(P, T, 0.0, b) for b in TOs]; R = [part(P, T, ϕ, b) for b in TOs]
    ok = [a.converged && b.converged for (a, b) in zip(R0, R)]
    x = [o.FeR_net for o in R0[ok]]; y = [o.FMQ for o in R[ok]]; o = sortperm(x)
    return x[o], y[o]
end
SRC = (0.03, 0.04)                                                   # fertile-mantle Fe³⁺/ΣFe, Canil & O'Neill (1996) range used in TOex_bench.md

g = Figure(size=(1100, 760), fontsize=16)
c2 = Axis(g[1, 1], xlabel="Fe³⁺/ΣFe of the source (same budget, no melt)", ylabel="ΔFMQ of the partially molten cell",
          title="MORB source melting at 1.5 GPa, 1350 °C: model (pyrolite, 30 ppm C, 100 ppm H₂O) against MORB glasses", titlealign=:left)
vspan!(c2, SRC..., color=(:gray, 0.15))
for (lab, μ, σ, col) in (("MORB glass −0.17 ± 0.15 (Zhang et al. 2018 calibration, n=166)", -0.17, 0.15, "#2B7BD6"),
                         ("MORB glass +0.19 ± 0.36 (Berry et al. 2018 calibration, n=42)", 0.19, 0.36, "#E8692C"))
    hspan!(c2, μ - σ, μ + σ, color=(col, 0.15)); hlines!(c2, [μ], color=col, linewidth=2, label=lab)
end
for (ϕ, ls) in ((0.05, :dot), (0.10, :solid), (0.20, :dash))
    x, y = curve(1.5, 1623.15, ϕ); lines!(c2, x, y, color=:black, linewidth=2.2, linestyle=ls, label="model, $(round(Int, 100ϕ)) % melt")
end
xlims!(c2, 0, 0.08); ylims!(c2, -4, 3)
Legend(g[2, 1], c2, nbanks=2, framevisible=false, labelsize=13)
Label(g[3, 1], "grey: fertile-mantle Fe³⁺/ΣFe 0.03–0.04 (Canil & O'Neill 1996). MORB values: Cottrell et al. (2021) Table 3.1, relative to QFM at 1 atm;\nthe model's FMQ is computed from sb24.",
      fontsize=13, color=:gray40)
rowsize!(g.layout, 1, Relative(0.70)); save(joinpath(OUT, "bench_C.png"), g, px_per_unit=2); println("wrote bench_C.png")
Finalize_MAGEMin(data); println("done")
