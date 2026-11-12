# Appendix figures of TOex_solver_flow.tex (figs/behaviour): capacity maps (A1, A2), melt water (A4) and the
#  oxygen, water, carbon and temperature sweeps at 3 GPa (B1-B5), all from the current partition_Oₑₓ.
# Run: julia --project=<StagPP.jl> behaviour_figs.jl   (or include after a Makie-free TOexLite loader)

isdefined(Main, :TOexLite) || @eval using StagPP
const M = isdefined(Main, :TOexLite) ? Main.TOexLite : StagPP
using CairoMakie, MAGEMin_C, Printf

const OUT = joinpath(@__DIR__, "figs", "behaviour"); mkpath(OUT)
const P, TAD, PH = 3.0, 1636.0, 0.0                        # 3 GPa, adiabat, harzburgite
const CG, CC, CM, CF, CS = "#4D4D4D", "#DAA520", "#D64B3C", "#2B7BD6", "#1BB37C"   # graphite, carbonate, melt, fluid, solid
const SP = (("H₂O", :XH2O, "#2B7BD6"), ("CO₂", :XCO2, "#E8692C"), ("CH₄", :XCH4, "#1BA77A"), ("CO", :XCO, "#8A5BD6"), ("H₂", :XH2, "#C8417F"))

data = Initialize_MAGEMin("sb24", verbose=false); CACHE = Dict{Any,Any}()   # one multi_point per (P, T, bulk)
solve(T, ϕ, TO, C, H) = M.partition_Oₑₓ(P, T, PH, ϕ, TO, C, H; verbose=false, data=data, out=true, mag_cache=CACHE)

# ---------- labels inside stacked bands: one per run of the same text where the band is thick enough ----------
function band_labels!(ax, x, lo, hi, txt; minw=0.2, col=:white, xscale=identity)
    n = length(x); th = hi .- lo; i = 1
    u = xscale === log10 ? log10.(x) : collect(float.(x))              # position along the axis
    u0, u1 = minimum(u) + 0.13*(maximum(u) - minimum(u)), maximum(u) - 0.13*(maximum(u) - minimum(u))
    while i <= n
        if th[i] < minw; i += 1; continue; end
        j = i; t = txt(i)
        while j < n && th[j+1] >= minw && txt(j+1) == t; j += 1; end
        a, b = max(u[i], u0), min(u[j], u1)                             # the run, kept off the edges
        if b - a >= 0.10*(maximum(u) - minimum(u))
            k = argmin(abs.(u .- (a + b)/2)); k = clamp(k, i, j)
            text!(ax, x[k], (lo[k] + hi[k])/2, text=t, color=col, font=:bold, fontsize=14, align=(:center, :center))
        end
        i = j + 1
    end
end
dom(o, keys, names) = names[argmax([getfield(o, k) for k in keys])]

# ---------- one row of four panels ----------
function row!(f, r, xs, R, xlab; xcol1, ylab1="ΔFMQ (log units)", xscale=identity, xfluid=nothing, last=false, rowlab="")
    xb = xfluid === nothing ? xs : xfluid
    a1 = Axis(f[r, 1], ylabel=ylab1, xscale=xscale, xlabel=last ? xlab[1] : "", title=r == 1 ? xlab[3] : "", titlealign=:left)
    y1 = [o.FMQ for o in R]
    lines!(a1, xcol1, y1, color=:black, linewidth=2.2, label="ΔFMQ")
    lo, hi = extrema(y1); hi - lo < 0.5 && ylims!(a1, (lo + hi)/2 - 0.25, (lo + hi)/2 + 0.25)   # near-flat fO₂: no 10⁻⁵ ticks
    Label(f[r, 0], replace(rowlab, ", " => "\n"), rotation=π/2, font=:bold, tellheight=false)     # own column, clear of the y label
    a2 = Axis(f[r, 2], yscale=log10, ylabel="mole fraction in fluid", xscale=xscale, xlabel=last ? xlab[2] : "",
              title=r == 1 ? "fluid composition" : "", titlealign=:left)
    absent = [o.fluid_absent for o in R]
    for k in eachindex(xb)                                         # grey band where no fluid exists
        absent[k] || continue
        x0 = k > 1 ? (xb[k-1] + xb[k])/2 : xb[k]; x1 = k < length(xb) ? (xb[k] + xb[k+1])/2 : xb[k]
        xscale === log10 && (x0 = k > 1 ? sqrt(xb[k-1]*xb[k]) : xb[k]; x1 = k < length(xb) ? sqrt(xb[k]*xb[k+1]) : xb[k])
        vspan!(a2, x0, x1, color=(:gray, 0.15))
    end
    for (nm, key, c) in SP
        lines!(a2, xb, max.([getfield(o, key) for o in R], 1e-6), color=c, linewidth=2, label=nm)
    end
    ylims!(a2, 1e-4, 1.6)
    # carbon
    a3 = Axis(f[r, 3], ylabel="share of carbon", xscale=xscale, xlabel=last ? xlab[2] : "", title=r == 1 ? "carbon" : "", titlealign=:left)
    g = [o.fC_gph for o in R]; cb = [o.fC_carb for o in R]; m = [o.fC_melt for o in R]; fl = [o.fC_fluid for o in R]
    s = max.(g .+ cb .+ m .+ fl, 1e-30); g, cb, m, fl = g ./ s, cb ./ s, m ./ s, fl ./ s
    l1 = g; l2 = l1 .+ cb; l3 = l2 .+ m; l4 = l3 .+ fl
    band!(a3, xb, zeros(length(xb)), l1, color=CG, label="graphite/diamond"); band!(a3, xb, l1, l2, color=CC, label="carbonate")
    band!(a3, xb, l2, l3, color=CM, label="melt"); band!(a3, xb, l3, l4, color=CF, label="fluid")
    band_labels!(a3, xb, zeros(length(xb)), l1, _ -> "graphite"; xscale=xscale)
    band_labels!(a3, xb, l1, l2, _ -> "carbonate"; col=:black, xscale=xscale)
    band_labels!(a3, xb, l2, l3, _ -> "C in melt"; xscale=xscale)
    band_labels!(a3, xb, l3, l4, i -> "C in fluid\nas " * dom(R[i], (:XCO2, :XCO, :XCH4), ("CO₂", "CO", "CH₄")); xscale=xscale)
    ylims!(a3, 0, 1)
    # hydrogen
    a4 = Axis(f[r, 4], ylabel="share of hydrogen", xscale=xscale, xlabel=last ? xlab[2] : "", title=r == 1 ? "hydrogen" : "", titlealign=:left)
    hs = [o.fH_solid for o in R]; hm = [o.fH_melt for o in R]; hf = [o.fH_fluid for o in R]
    s = max.(hs .+ hm .+ hf, 1e-30); hs, hm, hf = hs ./ s, hm ./ s, hf ./ s
    k1 = hs; k2 = k1 .+ hm; k3 = k2 .+ hf
    band!(a4, xb, zeros(length(xb)), k1, color=CS, label="solid"); band!(a4, xb, k1, k2, color=CM, label="melt")
    band!(a4, xb, k2, k3, color=CF, label="fluid")
    band_labels!(a4, xb, zeros(length(xb)), k1, _ -> "solid"; xscale=xscale)
    band_labels!(a4, xb, k1, k2, _ -> "melt"; xscale=xscale)
    band_labels!(a4, xb, k2, k3, i -> "fluid\n(" * dom(R[i], (:XH2O, :XH2, :XCH4), ("H₂O", "H₂", "CH₄")) * ")"; xscale=xscale)
    ylims!(a4, 0, 1)
    return a1, a2, a3, a4
end

function finish!(f, axs, title, name)
    Label(f[0, 1:4], title, font=:bold, fontsize=20)
    Legend(f[4, 1], axs[1], framevisible=false, tellwidth=false)                 # each legend under its own column
    Legend(f[4, 2], axs[2], nbanks=3, framevisible=false, tellwidth=false)
    Legend(f[4, 3], axs[3], nbanks=2, framevisible=false, tellwidth=false)
    Legend(f[4, 4], axs[4], nbanks=1, orientation=:horizontal, framevisible=false, tellwidth=false)
    Label(f[5, 2], "grey band: no fluid; species give the composition of the incipient fluid", fontsize=12, color=:gray40, tellwidth=false)
    rowgap!(f.layout, 4, 14); rowgap!(f.layout, 5, 10)   # layout rows: title, 3 panel rows, legends, note
    save(joinpath(OUT, name), f, px_per_unit=2); println("wrote ", name); flush(stdout)
end
sweep(f, xs) = (R = [f(x) for x in xs]; ok = [o.converged for o in R]; @printf("  %d/%d converged\n", count(ok), length(ok)); flush(stdout); (xs[ok], R[ok]))

ROWS = (("reduced, TOex -0.05 wt%", -5e-4), ("graphite-buffered, TOex 0.02 wt%", 2e-4), ("oxidized, TOex 0.3 wt%", 3e-3))

# ---------- B1: excess oxygen swept, three melt fractions ----------
println("B1"); flush(stdout)
f = Figure(size=(1900, 1300), fontsize=14); axs = nothing
for (r, (ϕ, rl)) in enumerate(((0.0, "solid + fluid"), (0.01, "solid + 1 % melt"), (0.1, "solid + 10 % melt")))
    xs, R = sweep(TO -> solve(TAD, ϕ, 1e-2TO, 5e-4, 1e-3), collect(range(-0.6, 0.6, length=81)))
    o = sortperm([q.FMQ for q in R]); xs, R = xs[o], R[o]
    a = row!(f, r, xs, R, ("excess oxygen TOex (wt%)", "ΔFMQ", "fO₂ against the oxygen budget"); xcol1=xs, xfluid=[q.FMQ for q in R],
             last=r == 3, rowlab=rl)
    r == 1 && (global axs = a)
end
finish!(f, axs, "Harzburgite at 3 GPa, $(Int(TAD)) K (adiabat): 0.05 wt% C, 0.1 wt% H₂O, excess oxygen swept", "B1_fo2_3GPa.png")

# ---------- B2, B3, B5: water, carbon and temperature swept, three oxygen budgets ----------
function budget_fig(name, title, xlabel, xs0, solver; xscale=log10)
    println(name); flush(stdout)
    f = Figure(size=(1900, 1300), fontsize=14); axs = nothing
    for (r, (rl, TO)) in enumerate(ROWS)
        xs, R = sweep(x -> solver(x, TO), xs0)
        a = row!(f, r, xs, R, (xlabel, xlabel, "fO₂"); xcol1=xs, xscale=xscale, last=r == 3, rowlab=rl)
        r == 1 && (axs = a)
    end
    finish!(f, axs, title, name)
end
budget_fig("B2_H_solid.png", "Water budget swept: harzburgite at 3 GPa, $(Int(TAD)) K, solid + fluid, 0.05 wt% C",
           "water budget (wt% H₂O)", 10 .^ range(-3, log10(3), length=50), (x, TO) -> solve(TAD, 0.0, TO, 5e-4, 1e-2x))
budget_fig("B3_C_solid.png", "Carbon budget swept: harzburgite at 3 GPa, $(Int(TAD)) K, solid + fluid, 0.1 wt% H₂O",
           "carbon budget (wt% C)", 10 .^ range(-4, 0, length=50), (x, TO) -> solve(TAD, 0.0, TO, 1e-2x, 1e-3))
budget_fig("B5_T_solid.png", "Temperature swept: harzburgite at 3 GPa, solid + fluid, 0.05 wt% C, 0.1 wt% H₂O",
           "T (K)", collect(range(900, 2000, length=45)), (x, TO) -> solve(x, 0.0, TO, 5e-4, 1e-3); xscale=identity)
Finalize_MAGEMin(data)

# ---------- A1, A2: solid capacity maps from StagH2Oa ----------
tb = M.load_StagH2Oa(); AWS = (1.0, 0.1, 0.01, 0.001); SN = ("Upper mantle", "Transition zone", "Lower mantle")
for (c, comp, name) in ((1, "Harzburgite", "A1_harz.png"), (length(tb.pX), "Basalt", "A2_basalt.png"))
    f = Figure(size=(1500, 1500), fontsize=14); hm = nothing
    for (r, a) in enumerate(AWS), q in 1:3
        ka = argmin(abs.(log.(tb.aw) .- log(a)))
        P0, P1, T0, T1 = tb.sectors[q]
        Ps = range(P0, P1, length=tb.nP); Ts = range(T0, T1, length=tb.nT)
        Z = log10.(max.(tb.S[:, :, q, c, ka], 1e-3))
        ax = Axis(f[r, q], title=(r == 1 ? "$(SN[q]),  " : "") * "a = $(a)", titlealign=:left,
                  xlabel=r == 4 ? "P (GPa)" : "", ylabel=q == 1 ? "T (K)" : "")
        hm = heatmap!(ax, Ps, Ts, Z, colormap=:Blues, colorrange=(-3, 1.3))
        contour!(ax, Ps, Ts, Z, levels=-2:1, color=(:white, 0.8), linewidth=0.8)
    end
    Colorbar(f[1:4, 4], hm, label="log₁₀ storage capacity (wt% H₂O)")
    Label(f[0, 1:3], "$(comp): solid storage capacity at water activity a (StagH2Oa, white contours every decade)", font=:bold, fontsize=18)
    save(joinpath(OUT, name), f, px_per_unit=2); println("wrote ", name)
end

# ---------- A4: melt water solubility (Sun & Yao), clamped at 2 GPa and 933 K ----------
f = Figure(size=(1500, 560), fontsize=14); hm = nothing
Ps = range(0.05, 7, length=180); Ts = range(700, 2200, length=160)
for (k, a) in enumerate((1.0, 0.1, 0.01))
    Z = [log10(sum(M.melt_h2o_wt(p, t, a))) for p in Ps, t in Ts]
    ax = Axis(f[1, k], title="melt,  a = $(a)", titlealign=:left, xlabel="P (GPa)", ylabel=k == 1 ? "T (K)" : "")
    global hm = heatmap!(ax, Ps, Ts, Z, colormap=:Blues, colorrange=(-0.5, 1.7))
    contour!(ax, Ps, Ts, 10 .^ Z, levels=[0.5, 1, 2, 5, 10, 20, 40], color=:black, linewidth=1, labels=true)
    vlines!(ax, [2.0], color=:gray30, linestyle=:dash); hlines!(ax, [933.15], color=:gray30, linestyle=:dash)
end
Colorbar(f[1, 4], hm, label="log₁₀ melt water solubility (wt% H₂O)")
Label(f[2, 1:3], "contours in wt% H₂O; dashed: 2 GPa pressure cap and 933 K temperature floor of the law", color=:gray40)
save(joinpath(OUT, "A4_melt.png"), f, px_per_unit=2); println("wrote A4_melt.png")
println("done")
