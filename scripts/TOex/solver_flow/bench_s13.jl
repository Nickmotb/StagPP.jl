# Benchmark slide (page 27 of TOex_solver_flow.tex) against Stagno et al. (2013, Nature 493, 84): (a) garnet
# oxybarometry run by run from their Supplementary Tables S2 and S5, (b) closed-system decompression as their Fig. 3.
# Run: include after the Makie-free TOexLite loader, or julia --project=<StagPP.jl> bench_s13.jl

isdefined(Main, :TOexLite) || @eval using StagPP
const M = isdefined(Main, :TOexLite) ? Main.TOexLite : StagPP
using CairoMakie, MAGEMin_C, Printf, Statistics, DelimitedFiles, Interpolations

const OUT = joinpath(@__DIR__, "figs", "bench"); mkpath(OUT)
const PB = 0.2                                             # model bulk: pyrolite (20 % basalt end-member)

# ---------------- (a) garnet Fe³⁺/ΣFe → fO₂ ----------------
raw  = readdlm(joinpath(@__DIR__, "data", "stagno_2013.csv"), ',', Any; header=true)[1]
runs = [(run=string(r[1]), P=Float64(r[2]), TC=Float64(r[3]), flux=string(r[4]), g=Float64(r[5]), ge=Float64(r[6]),
         f=Float64(r[7]), fe=Float64(r[8]), f3=Float64(r[9])) for r in eachrow(raw)]
(; X, Xox) = M.endmember_composition(PB)
Xv   = Vector{String}(Xox); Xd = (X=zeros(length(Xv)), Xox=copy(Xv), mm=M.get_Xoxmm(Xv))
Xred = M.oxidized_Xox(Xv, "Fe_O")
rs   = collect(range(0.0, 0.30, length=31))                # bulk Fe³⁺/ΣFe grid of the solid map
PT   = unique([(r.P, r.TC) for r in runs])
pts  = [(P, T, r) for (P, T) in PT for r in rs]            # one multi_point call over every (P, T, r)
Xl   = [M.oxidize_bulk(collect(X), r, Xd; wt_in=true, wt_out=true, frac=true, FeFormat="Fe_O", vector=true) for (_, _, r) in pts]
data = Initialize_MAGEMin("sb24", verbose=false)
println("MAGEMin: $(length(pts)) points"); flush(stdout)
mo   = multi_point_minimization([10p[1] for p in pts], [p[2] for p in pts], data, X=Xl, Xoxides=Xred, name_solvus=true, sys_in="wt", progressbar=false)
# garnet Fe³⁺/ΣFe from sb24 end-members: andradite Ca₃Fe³⁺₂Si₃O₁₂, almandine Fe²⁺₃Al₂Si₃O₁₂
function gtR(o)
    k = findfirst(s -> "andr" in s.emNames && "alm" in s.emNames, o.SS_vec); isnothing(k) && return NaN
    s = o.SS_vec[k]; a = s.emFrac[findfirst(==("andr"), s.emNames)]; l = s.emFrac[findfirst(==("alm"), s.emNames)]
    return 2a/(2a + 3l)
end
gmap = Dict{Tuple{Float64,Float64},Any}()
for (P, T) in PT
    idx = findall(p -> p[1] == P && p[2] == T, pts)
    g = gtR.(mo[idx]); f = [o.fO2 for o in mo[idx]]; q = [o.fO2 - o.dQFM for o in mo[idx]]
    ok = isfinite.(g) .& isfinite.(f); g, f = g[ok], f[ok]; o = sortperm(g)
    gmap[(P, T)] = (g=g[o], f=f[o], qfm=median(q[isfinite.(q)]))
    @printf("P %.0f GPa T %.0f C: garnet Fe3+/SFe %.3f-%.3f, monotone %s\n", P, T, minimum(g), maximum(g), issorted(f[o]))
end
model_f(r, g) = (m = gmap[(r.P, r.TC)]; linear_interpolation(m.g, m.f, extrapolation_bc=Line())(g))
fm  = [model_f(r, r.g) for r in runs]
fmh = [model_f(r, r.g + r.ge) for r in runs]; fml = [model_f(r, max(r.g - r.ge, 1e-3)) for r in runs]
st(d) = (mean(d), sqrt(mean(d.^2)))
mM, rM = st(fm .- [r.f for r in runs]); m3, r3 = st([r.f3 - r.f for r in runs])
@printf("(a) model - measured: mean %+.2f rms %.2f | S13 eq.(3) - measured: mean %+.2f rms %.2f\n", mM, rM, m3, r3)
for (r, x) in zip(runs, fm); @printf("   %-9s %d GPa %4.0f C gt %.3f  meas %.2f  model %.2f  eq3 %.2f\n", r.run, r.P, r.TC, r.g, r.f, x, r.f3); end

# ---------------- (b) closed-system decompression, their Fig. 3 ----------------
CACHE = Dict{Any,Any}()
zP(z) = 0.0316z                                            # GPa from depth (km), PREM upper mantle
zT(z) = 1320 + 273.15 + 0.4z                               # 1320 °C adiabat, 0.4 K/km (K)
wFeO  = (PB*0.0718 + (1 - PB)*0.0835)                       # bulk FeO mass fraction
TOr(r) = wFeO/71.844*(r/2)*15.999                          # excess O (mass fraction) of bulk Fe³⁺/ΣFe = r
const TCb, Hb = 3e-5, 1e-4                                 # 30 ppm C, 100 ppm H₂O
part(z, TO) = M.partition_Oₑₓ(zP(z), zT(z), PB, 0.0, TO, TCb, Hb; verbose=false, data=data, out=true, mag_cache=CACHE)
zs  = collect(250.0:-5.0:60.0)
R0s = [0.02, 0.03, 0.04, 0.05, 0.06, 0.07]
paths = Dict{Float64,Any}()
for r0 in R0s
    t = @elapsed R = [part(z, TOr(r0)) for z in zs]
    paths[r0] = R
    @printf("(b) r0 %.2f: %d/%d converged in %.1f s; FeR_net at 250 km %.4f\n", r0, count(o -> o.converged, R), length(R), t, R[1].FeR_net); flush(stdout)
end
fE  = [M.logfO2_EMOG(zT(z), 1e4zP(z)) for z in zs]
qfm = [paths[0.04][i].fO2 - paths[0.04][i].FMQ for i in eachindex(zs)]
# onset and end of graphite oxidation for 4 %: graphite share first below 99 % and below 1 % (some C sits in the fluid)
p4  = paths[0.04]; fg = [o.fC_gph for o in p4]
i1  = findfirst(<(0.99), fg); i2 = findfirst(<(0.01), fg)
zon = isnothing(i1) ? NaN : zs[i1]; zend = isnothing(i2) ? NaN : zs[i2]
dR  = isnothing(i2) ? NaN : p4[i1].FeR_net - p4[i2].FeR_net
lin = findall(i -> zs[i] > (isnan(zon) ? 0 : zon) + 10, eachindex(zs))
slope = length(lin) > 2 ? -(p4[lin[end]].FMQ - p4[lin[1]].FMQ)/(zs[lin[end]] - zs[lin[1]]) : NaN
@printf("(b) 4%% carbon at %.0f km: graphite %.2f carbonate %.2f fluid %.2f\n", zend, p4[i2].fC_gph, p4[i2].fC_carb, p4[i2].fC_fluid)
@printf("(b) 4%%: carbonate from %.0f km, all C oxidised at %.0f km, bulk Fe3+/SFe drop %.4f, dFMQ/dz %.4f per km (1/%.0f)\n",
        zon, zend, dR, slope, 1/slope)

# ---------------- figure ----------------
fig = Figure(size=(1500, 760), fontsize=16)
a = Axis(fig[1, 1], xlabel="measured log₁₀ fO₂ (Ir–Fe sensor)", ylabel="predicted log₁₀ fO₂ at the measured garnet Fe³⁺/ΣFe",
         title="(a) Garnet oxybarometry, run by run", titlealign=:left, xgridvisible=false, ygridvisible=false)
lims = (-8.0, -2.5); lines!(a, [lims...], [lims...], color=:gray50, linestyle=:dash); limits!(a, lims, lims)
band!(a, [lims...], [lims[1]-0.5, lims[2]-0.5], [lims[1]+0.5, lims[2]+0.5], color=(:gray, 0.12))
cols = Dict(3.0 => "#2B7BD6", 6.0 => "#E8692C", 7.0 => "#7E57C2")
scatter!(a, [r.f for r in runs], [r.f3 for r in runs], color=:white, strokecolor=:gray45, strokewidth=1.5, markersize=11)
for r in runs, (i, x) in ((findfirst(==(r), runs), 0),)
    rangebars!(a, [r.f], [fml[i]], [fmh[i]], color=(cols[r.P], 0.5), whiskerwidth=6)
    scatter!(a, [r.f], [fm[i]], color=cols[r.P], marker=r.flux == "H2O" ? :utriangle : :circle, markersize=14,
             strokecolor=:white, strokewidth=0.5)
end
dP(P) = (d = [fm[i] - runs[i].f for i in eachindex(runs) if runs[i].P == P]; @sprintf("%d GPa: %+.2f ± %.2f (n=%d)", P, mean(d), std(d), length(d)))
elems = [MarkerElement(marker=:circle, color=cols[3.0], markersize=14), MarkerElement(marker=:circle, color=cols[6.0], markersize=14),
         MarkerElement(marker=:circle, color=cols[7.0], markersize=14), MarkerElement(marker=:utriangle, color=:gray30, markersize=14),
         MarkerElement(marker=:circle, color=:white, strokecolor=:gray45, strokewidth=1.5, markersize=11)]
Legend(fig[2, 1], elems, [dP(3.0), dP(6.0), dP(7.0), "H₂O flux (circles: carbonate)",
       @sprintf("their eq. (3): %+.2f, rms %.2f", m3, r3)], nbanks=2, framevisible=false, labelsize=13, tellwidth=false)
text!(a, -7.85, -2.8, text=@sprintf("model (sb24, pyrolite): mean %+.2f, rms %.2f, n = %d", mM, rM, length(runs)),
      align=(:left, :top), fontsize=14, font=:bold)

b = Axis(fig[1, 2], xlabel="depth (km)", ylabel="ΔFMQ", title="(b) Upwelling mantle, 30 ppm C: TOex (solid), Stagno et al. (dashed)",
         titlealign=:left, xgridvisible=false, ygridvisible=false)
lines!(b, zs, fE .- qfm, color="#1F9E89", linewidth=3, label="model EMOG(D)")
for r0 in R0s
    R = paths[r0]; ok = [o.converged for o in R]
    lines!(b, zs[ok], [o.FMQ for o in R[ok]], color=r0 == 0.04 ? "#E8692C" : :gray45, linewidth=r0 == 0.04 ? 4 : 2.5)
    k = findfirst(ok); text!(b, zs[k] + 2, R[k].FMQ, text=@sprintf("%d %%", round(Int, 100r0)), align=(:left, :center), fontsize=12, color=:gray30)
end
# Stagno et al. (2013) Fig. 3a, extracted from the PDF's vector paths: thin dashed, under the model curves' colours
sf = readdlm(joinpath(@__DIR__, "data", "stagno_2013_fig3a.csv"), ',', Any; header=true, comments=true)[1]
for c in unique(sf[:, 1])
    c == "FeNi" && continue
    if c != "buffer"                                       # TOex − Stagno at the same Fe³⁺/ΣFe
        R = paths[parse(Int, chop(c))/100]; sit = linear_interpolation(Float64.(sf[sf[:, 1] .== c, 2]), Float64.(sf[sf[:, 1] .== c, 3]))
        @printf("(b) %s TOex - Stagno: 250 km %+.2f, 150 km %+.2f, 70 km %+.2f\n", c, (R[findfirst(==(zz), zs)].FMQ - sit(zz) for zz in (250.0, 150.0, 70.0))...)
    end
    i = sf[:, 1] .== c
    col = c == "buffer" ? :black : c == "4%" ? "#E8692C" : :gray45
    lines!(b, Float64.(sf[i, 2]), Float64.(sf[i, 3]), color=col, linewidth=1.3, linestyle=:dash)
end
xlims!(b, 50, 260); ylims!(b, -5, 1.5)
Legend(fig[2, 2], [LineElement(color="#E8692C", linewidth=4), LineElement(color=:gray45, linewidth=2.5),
                   LineElement(color="#1F9E89", linewidth=3),
                   LineElement(color=:gray45, linestyle=:dash, linewidth=1.3), LineElement(color=:black, linestyle=:dash, linewidth=1.3)],
       [@sprintf("model, 4 %% (ΔFMQ +1/%.0f per km)", 1/slope), "2, 3, 5, 6, 7 %", "model EMOG(D)",
        "Stagno et al. (2013) Fig. 3a, same Fe³⁺/ΣFe", "Stagno et al.: graphite + carbonate melt"], nbanks=2, framevisible=false, labelsize=13, tellwidth=false)
rowsize!(fig.layout, 1, Relative(0.78))
save(joinpath(OUT, "bench_S13.png"), fig, px_per_unit=2); println("wrote bench_S13.png")
Finalize_MAGEMin(data); println("done")
