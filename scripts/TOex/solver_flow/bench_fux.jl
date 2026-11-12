# Benchmark slide (page 28 of TOex_solver_flow.tex) against Fux et al. (2024, EPSL 645, 118955): closed-system fO₂ (ΔIW)
# and metal along their 5–16 GPa path for bulk Fe³⁺/ΣFe 0.036 at 5 GPa (their Fig. 8, extracted from the PDF).
# Run: include after the Makie-free TOexLite loader, or julia --project=<StagPP.jl> bench_fux.jl

isdefined(Main, :TOexLite) || @eval using StagPP
const M = isdefined(Main, :TOexLite) ? Main.TOexLite : StagPP
using CairoMakie, MAGEMin_C, Printf, DelimitedFiles, Interpolations

const OUT = joinpath(@__DIR__, "figs", "bench"); mkpath(OUT)
const PB = 0.2                                              # pyrolite
Tf   = linear_interpolation([5.0, 9, 12, 16], [1420.0, 1450, 1490, 1590] .+ 273.15)   # their run temperatures (K)
wFeO = PB*0.0718 + (1 - PB)*0.0835                          # bulk FeO mass fraction
TOr(r) = wFeO/71.844*(r/2)*15.999                           # excess O (mass fraction) of bulk Fe³⁺/ΣFe = r
const TCb, Hb = 1.2e-4, 1e-4                                # 120 ppm C (their BSE 100–150), 100 ppm H₂O

data = Initialize_MAGEMin("sb24", verbose=false); CACHE = Dict{Any,Any}()
Ps = collect(5.0:0.5:16.0)
R  = [M.partition_Oₑₓ(P, Tf(P), PB, 0.0, TOr(0.036), TCb, Hb; verbose=false, data=data, out=true, mag_cache=CACHE) for P in Ps]
Finalize_MAGEMin(data)
ok  = [o.converged for o in R]
iw  = [o.IW for o in R]
fe0 = [o.FeM_solid*wFeO/71.844*55.845*1e6 for o in R]       # metal Fe⁰, ppm by mass of rock
for (P, o, m) in zip(Ps, R, fe0)
    @printf("P %4.1f dIW %+.2f FeR_solid %.4f FeM_solid %.4f Fe0 %5.0f ppm gph %.2f\n", P, o.IW, o.FeR_solid, o.FeM_solid, m, o.fC_gph)
end

fx  = readdlm(joinpath(@__DIR__, "data", "fux_2024_fig8.csv"), ',', Any; header=true, comments=true)[1]
ic  = fx[:, 1] .== "curve"; ib = fx[:, 1] .== "band"
fP, fD = Float64.(fx[ic, 2]), Float64.(fx[ic, 3])
fit = linear_interpolation(fP, fD, extrapolation_bc=Flat())
d   = [iw[i] - fit(Ps[i]) for i in eachindex(Ps) if ok[i]]
@printf("TOex - Fux: mean %+.2f, max |d| %.2f; at 5/8/12/16 GPa: %s\n", sum(d)/length(d), maximum(abs.(d)),
        join([@sprintf("%+.2f", iw[findfirst(==(p), Ps)] - fit(p)) for p in (5.0, 8.0, 12.0, 16.0)], " "))
zIW = Ps[findfirst(<(0), iw)]; @printf("TOex crosses IW at %.1f GPa\n", zIW)

fig = Figure(size=(1500, 760), fontsize=16)
a = Axis(fig[1, 1], xlabel="P (GPa)", ylabel="ΔIW", title="(a) fO₂ of a closed BSE-like mantle, Fe³⁺/ΣFe 0.036 at 5 GPa",
         titlealign=:left, xgridvisible=false, ygridvisible=false)
poly!(a, Point2f.(Float64.(fx[ib, 2]), Float64.(fx[ib, 3])), color=(:gray, 0.18))
lines!(a, fP, fD, color=:black, linewidth=1.3, linestyle=:dash)
scatter!(a, fP, fD, color=:white, strokecolor=:black, strokewidth=1, markersize=8)
lines!(a, Ps[ok], iw[ok], color="#E8692C", linewidth=4)
hlines!(a, [0.0], color=:gray60, linewidth=1)
xlims!(a, 4.5, 16.5); ylims!(a, -2, 3.5)
Legend(fig[2, 1], [LineElement(color="#E8692C", linewidth=4),
                   [LineElement(color=:black, linestyle=:dash, linewidth=1.3), MarkerElement(marker=:circle, color=:white, strokecolor=:black, strokewidth=1, markersize=8)],
                   PolyElement(color=(:gray, 0.18))],
       [@sprintf("TOex (sb24, pyrolite, 120 ppm C): reaches IW at %.1f GPa", zIW), "Fux et al. (2024) Fig. 8, alloy redox sensor",
        "their ±0.5 uncertainty"], nbanks=2, framevisible=false, labelsize=13, tellwidth=false)

b = Axis(fig[1, 2], xlabel="P (GPa)", ylabel="reduced phase (ppm by mass)", title="(b) Metal made by Fe²⁺ disproportionation",
         titlealign=:left, xgridvisible=false, ygridvisible=false)
lines!(b, Ps[ok], fe0[ok], color="#E8692C", linewidth=4)
rangebars!(b, [5.0], [600.0], [750.0], color=:black, whiskerwidth=10, linewidth=1.5)
scatter!(b, [5.0, 16.0], [675.0, 3800.0], color=:white, strokecolor=:black, strokewidth=1.5, markersize=12)
xlims!(b, 4.5, 16.5); ylims!(b, 0, 4200)
Legend(fig[2, 2], [LineElement(color="#E8692C", linewidth=4), MarkerElement(marker=:circle, color=:white, strokecolor=:black, strokewidth=1.5, markersize=12)],
       ["TOex: Fe⁰ metal (sb24 has no Ni or S)", "Fux et al.: whole Fe-Ni-S-C(-O) alloy melt, 250 ppm S"],
       nbanks=1, framevisible=false, labelsize=13, tellwidth=false)
rowsize!(fig.layout, 1, Relative(0.80))
save(joinpath(OUT, "bench_Fux.png"), fig, px_per_unit=2); println("wrote bench_Fux.png")
