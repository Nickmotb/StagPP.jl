# Methods figures of the water paper: hydrous solidus depression (Katz 2003 + TZ extension) and hydrous density decrease.
# Run from the StagPP.jl root: julia --project=. scripts/Paper1/methods_figs.jl   (OUT=<dir> to redirect the output)
using CairoMakie, Printf

OUT = get(ENV, "OUT", "/Users/nickmotb/Desktop/Drafts/Water/figures/methods")
TAB = "/Users/nickmotb/Desktop/watercarb/Perple_X_MAGEMin/MAGEMin/Tables"
ink = colorant"#0b0b0b"; muted = colorant"#52514e"; gshade = colorant"#f3f2ee"
seq = cgrad([colorant"#f3f7fd", colorant"#cde2fb", colorant"#9ec5f4", colorant"#6da7ec", colorant"#3987e5", colorant"#256abf", colorant"#184f95", colorant"#0d366b"])
thm = Theme(fontsize=12, textcolor=ink, Axis=(xgridvisible=false, ygridvisible=false, xlabelsize=13.5, ylabelsize=13.5, xticklabelsize=11.5, yticklabelsize=11.5,
            xticksmirrored=true, yticksmirrored=true, xminorticksvisible=true, yminorticksvisible=true, xtickalign=1, ytickalign=1, xminortickalign=1, yminortickalign=1,
            spinewidth=0.9), Legend=(framevisible=false, labelsize=11, patchsize=(18, 10), rowgap=1))
include(joinpath(@__DIR__, "katz_fit.jl"))   # calibration data (lo02, lo02e, kh97, fr17), model and fitted θ = (Pc, w, ginf, Psat)
PK = 8.0                                     # Katz et al. (2003) calibration limit

# Solidus depression (K) of bulk water w (wt%) at P (GPa) at the solidus (F = 0), with the fitted extension
katz_dT(P, w) = (Ph = min(P, fitres.θ[4]); 43 * min(w / 0.01, 12 * Ph^0.6 + Ph)^0.75 * gstep(P, fitres.θ[1:3]...))

# ---------------- Figure 1: solidus depression ----------------
# draw the brackets and markers of all three sets on axis a; returns the marker plots for the legend
function katz_data!(a; ms=1.0)
    xs(d) = [x[1] == 10.0 && d === lo02 ? 9.6 : x[1] for x in d]       # Litasov & Ohtani 10 GPa bar beside Kawamoto & Holloway
    for (d, c) in ((lo02, seq[1.0]), (kh97, seq[0.84]), (fr17, seq[0.7]))
        rangebars!(a, xs(d), getindex.(d, 2), getindex.(d, 3), color=c, whiskerwidth=6ms, linewidth=0.9)
    end
    mid(d) = [(x[2] + x[3]) / 2 for x in d]
    sl = scatter!(a, first.(lo02e), last.(lo02e), marker=:utriangle, markersize=10ms, color=seq[1.0], strokecolor=:white, strokewidth=0.8)
    sk = scatter!(a, first.(kh97), mid(kh97), marker=:diamond, markersize=10ms, color=seq[0.84], strokecolor=:white, strokewidth=0.8)
    sf = scatter!(a, first.(fr17), mid(fr17), marker=:circle, markersize=8.5ms, color=seq[0.7], strokecolor=:white, strokewidth=0.8)
    bar = [LineElement(points=Point2f[(0.5, 0), (0.5, 1)], color=muted, linewidth=0.9),
           LineElement(points=Point2f[(0.3, 0), (0.7, 0)], color=muted, linewidth=0.9), LineElement(points=Point2f[(0.3, 1), (0.7, 1)], color=muted, linewidth=0.9)]
    return [sl, sk, sf, bar, PolyElement(color=(seq[1.0], 0.16))]
end
labs = ["Litasov & Ohtani (2002), 2 wt%", "Kawamoto & Holloway (1997), sat.", "Freitas et al. (2017), > 0.3 wt%", "experimental bounds", "fit, 1σ (saturated)"]

f1 = with_theme(thm) do
    f = Figure(size=(620, 370), figure_padding=(6, 12, 6, 8))
    a = Axis(f[1, 1], xlabel="Pressure (GPa)", ylabel="Solidus depression, ΔT (K)", xticks=0:5:25, yticks=0:200:1000)
    xlims!(a, 0, 25.6); ylims!(a, 0, 1120)
    vspan!(a, 0, PK, color=gshade)
    vlines!(a, PK, color=muted, linestyle=(:dash, :dense), linewidth=0.9)
    text!(a, PK / 2, 1105, text="Katz et al. (2003)", align=(:center, :top), color=muted, font=:italic, fontsize=11)
    text!(a, PK + 0.25, 1105, text="extension, g(P)", align=(:left, :top), color=muted, font=:italic, fontsize=11)
    # transition zone bracket
    lines!(a, [13.5, 13.5, 23.5, 23.5], [1050, 1065, 1065, 1050], color=muted, linewidth=0.8)
    text!(a, 18.5, 1070, text="transition zone", align=(:center, :bottom), color=muted, fontsize=10.5)
    band!(a, fitres.P, fitres.lo, fitres.hi, color=(seq[1.0], 0.16))   # 1-sigma of the saturated curve, under the lines
    P = range(0.01, 25.6, 2000)
    ws = (0.01, 0.05, 0.1, 0.2, 2.0); cs = (0.42, 0.56, 0.7, 0.84, 1.0)
    for (w, c) in zip(ws, cs)
        lines!(a, P, katz_dT.(P, w), color=seq[c], linewidth=2.2)
    end
    for (w, c) in zip(ws[1:4], cs[1:4])
        text!(a, 4.6, katz_dT(4.6, w) + 10, text=@sprintf("%g wt%%", w), align=(:center, :bottom), color=seq[c], fontsize=11, font=:bold)
    end
    text!(a, PK - 0.4, katz_dT(PK, 2.0) + 20, text="H₂O-saturated", align=(:right, :bottom), color=seq[1.0], fontsize=11, font=:bold)
    axislegend(a, katz_data!(a; ms=1.05), labs, position=(1.0, 0.86), margin=(0, 6, 0, 0), rowgap=0)
    f
end
mkpath(OUT)
save(joinpath(OUT, "katz_depression.png"), f1, px_per_unit=4); save(joinpath(OUT, "katz_depression.pdf"), f1)

# ---------------- Figure 2: hydrous density decrease ----------------
# StagSdens/StagMdens: header nP nT nH / axis ends, then slot (1 basalt, 2 harzburgite) x H x P rows of T values
tabs = map(("StagSdens.dat", "StagMdens.dat")) do fn
    L = readlines(joinpath(TAB, fn)); ax = parse.(Float64, split(L[3])); n = parse.(Int, split(L[1]))
    A = permutedims(reshape(parse.(Float64, split(join(L[5:end], " "))), n[2], n[1], n[3], 2), (2, 1, 3, 4))
    (P=range(ax[1], ax[2], n[1]), T=range(ax[3], ax[4], n[2]), H=range(ax[5], ax[6], n[3]), A=A)
end
# Bilinear read in (T, H) at every pressure node, as thermochem_int does between nodes
function dens_TH(t, s, T, w)
    x = (T - t.T[1]) / step(t.T) + 1; i = clamp(floor(Int, x), 1, length(t.T) - 1); fx = x - i
    y = (w - t.H[1]) / step(t.H) + 1; j = clamp(floor(Int, y), 1, length(t.H) - 1); fy = y - j
    A = t.A
    return [(1 - fx) * (1 - fy) * A[p, i, j, s] + fx * (1 - fy) * A[p, i+1, j, s] + (1 - fx) * fy * A[p, i, j+1, s] + fx * fy * A[p, i+1, j+1, s] for p in axes(A, 1)]
end

ws2 = (1.0, 2.0, 4.0, 8.0); cs2 = (0.45, 0.62, 0.8, 1.0); Tref, Tlo, Thi = 2000.0, 1500.0, 2500.0
f2 = with_theme(thm) do
    f = Figure(size=(620, 540), figure_padding=(6, 8, 6, 6))
    axs = Matrix{Axis}(undef, 2, 2)
    for (r, (s, comp)) in enumerate(((2, "Harzburgite"), (1, "Basalt"))), (c, ph) in enumerate(("Solid", "Melt"))
        t = tabs[c]
        a = axs[r, c] = Axis(f[r+2, c], xticks=0:10:50, yticks=0:5:25,
                             xlabel=(r == 2 ? "Pressure (GPa)" : ""),
                             xticklabelsvisible=(r == 2), yticklabelsvisible=(c == 1))
        xlims!(a, 0, 50); ylims!(a, 0, 25)
        vspan!(a, 13.5, 23.5, color=gshade)
        r == 1 && text!(a, 18.5, 24.3, text="TZ", align=(:center, :top), color=muted, fontsize=10.5)
        for (w, cc) in zip(ws2, cs2)
            band = reduce(hcat, [dens_TH(t, s, T, w) for T in range(Tlo, Thi, 21)])
            band!(a, t.P, vec(minimum(band, dims=2)), vec(maximum(band, dims=2)), color=(seq[cc], 0.28))
            lines!(a, t.P, dens_TH(t, s, Tref, w), color=seq[cc], linewidth=1.9)
        end
        text!(a, 0.965, 0.955, text=string("abcd"[2(r-1)+c]), space=:relative, align=(:right, :top), font=:bold, fontsize=13)
        r == 1 && Label(f[2, c], ph, font=:bold, fontsize=12.5, tellwidth=false, padding=(0, 0, 0, 0))
        c == 2 && Label(f[r+2, 3], comp, font=:bold, fontsize=12.5, rotation=-π/2, tellheight=false, padding=(4, 0, 0, 0))
    end
    linkaxes!(axs...)
    Label(f[3:4, 0], "Density decrease relative to dry (%)", rotation=π/2, fontsize=13.5, tellheight=false)
    els = [[PolyElement(color=(seq[cc], 0.28)), LineElement(color=seq[cc], linewidth=1.9)] for cc in cs2]
    Legend(f[1, 1:2], els, [@sprintf("%g wt%% H₂O", w) for w in ws2], orientation=:horizontal, colgap=14, tellwidth=false)
    Label(f[1, 1:2, Bottom()], @sprintf("lines %d K, bands %d–%d K", Tref, Tlo, Thi), color=muted, fontsize=10.5, font=:italic, padding=(0, 0, 0, 2))
    colgap!(f.layout, 1, 2); colgap!(f.layout, 2, 16); colgap!(f.layout, 3, 4); rowgap!(f.layout, 1, 10); rowgap!(f.layout, 2, 2); rowgap!(f.layout, 3, 18)
    f
end
save(joinpath(OUT, "hydrous_density.png"), f2, px_per_unit=4); save(joinpath(OUT, "hydrous_density.pdf"), f2)

# ---------------- Figure 3: Katz extension (a), its calibration (b), solid (c) and melt (d) density decrease, harzburgite to basalt ----------------
f3 = with_theme(thm) do
    f = Figure(size=(730, 700), figure_padding=(4, 8, 4, 6))
    gt = f[1, 1] = GridLayout(); gc = gt[1, 2] = GridLayout(); gd = f[2, 1] = GridLayout()
    a = Axis(gt[1, 1], xlabel="Pressure (GPa)", ylabel="Solidus depression, ΔT (K)", xticks=0:5:25, yticks=0:200:1000)
    xlims!(a, 0, 25.6); ylims!(a, 0, 1120)
    vspan!(a, 0, PK, color=gshade)
    vlines!(a, PK, color=muted, linestyle=(:dash, :dense), linewidth=0.9)
    text!(a, PK - 0.25, 470, text="Katz et al. (2003)", rotation=π/2, align=(:left, :bottom), color=muted, font=:italic, fontsize=10)
    band!(a, fitres.P, fitres.lo, fitres.hi, color=(seq[1.0], 0.16))   # 1-sigma of the saturated curve, under the lines
    P = range(0.01, 25.6, 2000)
    ws = (0.01, 0.05, 0.1, 0.2, 2.0); cs = (0.42, 0.56, 0.7, 0.84, 1.0)
    for (w, c) in zip(ws, cs)
        lines!(a, P, katz_dT.(P, w), color=seq[c], linewidth=2.0)
    end
    for (w, c) in zip(ws[1:4], cs[1:4])
        text!(a, 4.4, katz_dT(4.4, w) + 10, text=@sprintf("%g wt%%", w), align=(:center, :bottom), color=seq[c], fontsize=10.5, font=:bold)
    end
    text!(a, PK - 0.4, katz_dT(PK, 2.0) + 20, text="H₂O-saturated", align=(:right, :bottom), color=seq[1.0], fontsize=10.5, font=:bold)
    text!(a, 24.6, 790, text="Smaller depression: the rising H₂O\nstorage capacity of the solid narrows\nthe solid–melt Gibbs energy gap",
          align=(:right, :top), color=muted, font=:italic, fontsize=9.5, justification=:right, lineheight=1.1)
    axislegend(a, katz_data!(a), labs, position=(1.0, 1.0), margin=(0, 4, 0, 4), labelsize=9.5, rowgap=-1, patchsize=(12, 10))

    # b: the logistic step with Psat fitted or held at Katz's 8 GPa, against the same bounds (legend: Akaike weights of the pair),
    #  then best value and 1σ range of each parameter of the adopted fit
    m = Axis(gc[1, 1], xlabel="Pressure (GPa)", ylabel="Saturated ΔT (K)", xticks=10:5:25, yticks=0:250:1000, xlabelsize=11.5, ylabelsize=11.5,
             xticklabelsize=10, yticklabelsize=10, xlabelpadding=1, xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5))
    xlims!(m, 7, 25.6); ylims!(m, 0, 1120)
    lgrey = colorant"#cfcdc7"
    rangebars!(m, [k <= length(lo02) && x[1] == 10.0 ? 9.6 : x[1] for (k, x) in enumerate(brk)], getindex.(brk, 2), getindex.(brk, 3), color=lgrey, whiskerwidth=4, linewidth=0.9)
    scatter!(m, first.(lo02e), last.(lo02e), marker=:utriangle, markersize=7, color=lgrey)
    Pm = range(7, 25.6, 800); pair = fitres.cmp[[1, 3]]                    # logistic, Psat fitted / Psat = 8 GPa
    w2 = (e = exp.(-([c[5] for c in pair] .- minimum(c[5] for c in pair)) ./ 2); e ./ sum(e))
    sty = ((color=seq[1.0], linewidth=2.2), (color=seq[1.0], linewidth=1.6, linestyle=(:dash, :dense)))
    for n in (2, 1)
        lines!(m, Pm, model.(Pm, Ref(pair[n][6])); sty[n]...)
    end
    axislegend(m, [LineElement(; sty[n]...) for n in 1:2], [rich("P", subscript("sat"), @sprintf(" fitted (adopted)   %.2f", w2[1])), rich("P", subscript("sat"), @sprintf(" = 8 GPa   %.2f", w2[2]))],
               position=:rt, margin=(0, 4, 0, 17), labelsize=9.5, rowgap=-1, patchsize=(16, 8), padding=(2, 2, 2, 2), backgroundcolor=:transparent)
    text!(m, 0.975, 0.975, text="Akaike weights", space=:relative, align=(:right, :top), font=:bold, fontsize=9.5)
    gs = gc[2, 1] = GridLayout()
    pn = (("drop centre", "P", "c", " = %.2f GPa"), ("drop width", "δ", "P", " = %.2f GPa"), ("deep residual", "g", "∞", " = %.2f"), ("saturation limit", "P", "sat", " = %.0f GPa"))
    xt = (11:1:16, 0:0.5:2.5, 0.1:0.1:0.4, 8:1:13); xl = ((11, 16), (0, 2.5), (0.05, 0.45), (8, 13.4))   # rows span the searched range
    for i in 1:4
        s = Axis(gs[i, 2], xticks=xt[i], xticklabelsize=9.5, xticksmirrored=false, xminorticks=IntervalsBetween(2), xticklabelpad=1,
                 yticksvisible=false, yticklabelsvisible=false, yminorticksvisible=false, yticksmirrored=false,
                 topspinevisible=false, leftspinevisible=false, rightspinevisible=false)
        xlims!(s, xl[i]...); ylims!(s, -1, 1)
        lo, hi = fitres.rng[i]; th = fitres.θ[i]
        lines!(s, [lo, hi], [0, 0], color=(seq[1.0], 0.4), linewidth=6)
        i == 4 && scatter!(s, [hi + 0.1], [0], marker=:rtriangle, markersize=9, color=(seq[1.0], 0.4))   # bound open above the searched range
        scatter!(s, [th], [0], color=seq[1.0], markersize=8.5, strokecolor=:white, strokewidth=0.8)
        Label(gs[i, 1], rich(rich(pn[i][1], color=muted, font=:italic), "  ", rich(pn[i][2], rich(pn[i][3], offset=(0, -0.3), fontsize=7.5)),
              Printf.format(Printf.Format(pn[i][4]), th)), halign=:right, fontsize=9.5, tellheight=false)
    end
    colgap!(gs, 1, 6); rowgap!(gs, 3); colsize!(gs, 1, Auto(0.8))
    rowgap!(gc, 1, 40); rowsize!(gc, 2, Relative(0.47))
    colsize!(gt, 1, Relative(0.565)); colgap!(gt, 1, 16)

    # c, d: one ribbon per water content, harzburgite (solid line) to basalt (dashed) at Tref; log axis spaces 1-2-4-8 wt% evenly
    ab = Axis[]
    for (r, ph) in enumerate(("Solid", "Melt"))
        t = tabs[r]
        d = Axis(gd[1, r], xticks=0:10:50, yscale=log10, yticks=[1, 2, 5, 10, 20], yminorticks=IntervalsBetween(5), xlabel="Pressure (GPa)",
                 ylabel=(r == 1 ? "Density decrease vs dry (%)" : ""), yticklabelsvisible=(r == 1))
        xlims!(d, 0, 50); ylims!(d, 0.55, 40); push!(ab, d)
        vspan!(d, 13.5, 23.5, color=gshade)
        for (w, cc) in zip(ws2, cs2)
            hz = dens_TH(t, 2, Tref, w); bs = dens_TH(t, 1, Tref, w)
            band!(d, t.P, min.(hz, bs), max.(hz, bs), color=(seq[cc], 0.3))
            lines!(d, t.P, hz, color=seq[cc], linewidth=1.6)
            lines!(d, t.P, bs, color=seq[cc], linewidth=1.3, linestyle=(:dash, :dense))
            k = findfirst(>=(46.0), t.P)
            text!(d, 49.2, 1.13 * max(hz[k], bs[k]), text=@sprintf("%g wt%%", w), align=(:right, :bottom), color=seq[cc], fontsize=9)
        end
        text!(d, 0.97, 0.955, text=@sprintf("%s, %d K", ph, Tref), space=:relative, align=(:right, :top), font=:bold, fontsize=11)
    end
    linkaxes!(ab...)
    axislegend(ab[1], [LineElement(color=ink, linewidth=1.6), LineElement(color=ink, linewidth=1.3, linestyle=(:dash, :dense))], ["harzburgite", "basalt"],
               position=(0.42, 0.0), margin=(0, 0, 1, 0), labelsize=9.5, patchsize=(16, 8), rowgap=-2, orientation=:horizontal, colgap=8)
    colgap!(gd, 1, 18); rowsize!(f.layout, 1, Relative(0.6)); rowgap!(f.layout, 1, 12)

    Makie.update_state_before_display!(f)
    # a: the arrow runs parallel to the saturated curve, 18 units off it along the normal (pixel space, so it keeps the drawn slope)
    W, H = a.scene.viewport[].widths; q = [Point2f(P * W / 25.6, katz_dT(P, 2.0) * H / 1120) for P in range(12.7, 14.8, 200)]
    tg = [(v = q[min(i + 1, end)] - q[max(i - 1, 1)]; v / hypot(v...)) for i in eachindex(q)]
    q = [q[i] + 18 * Point2f(-tg[i][2], tg[i][1]) for i in eachindex(q)]
    tip = q[end]; u = tg[end]; base = tip - 8 * u; nrm = Point2f(-u[2], u[1])
    lines!(a, push!(filter(x -> hypot((tip - x)...) > 8, q), base), space=:pixel, color=muted, linewidth=1.1)
    poly!(a, Point2f[tip, base + 3.5 * nrm, base - 3.5 * nrm], space=:pixel, color=muted)
    # b: key to the parameter rows, written at 45° above the dot and the 1σ bar start of the first row (figure space, so not clipped)
    s1 = contents(gs[1, 2])[1]; o = s1.scene.viewport[].origin; W1, H1 = s1.scene.viewport[].widths
    for (x, txt, al) in ((fitres.θ[1], "best value", 1.0), (fitres.rng[1][1], "1σ range", 0.45))
        text!(f.scene, Point2f(o[1] + (x - 11) / 5 * W1 - 2, o[2] + H1 / 2 + 7), text=txt, rotation=π/4, align=(:left, :bottom), color=(seq[1.0], al), fontsize=9.5)
    end
    # panel letters: equal 19-unit white boxes 5 units in from each top-left corner, ink centred (x-height letters sit lower)
    for (ax, s) in zip((a, m, ab...), ("a", "b", "c", "d"))
        h = ax.scene.viewport[].widths[2]; c = (5 + 9.5, h - 5 - 9.5)
        poly!(ax, Rect2f(c[1] - 9.5, c[2] - 9.5, 19, 19), space=:pixel, color=:white, strokecolor=muted, strokewidth=0.6)
        text!(ax, c[1], c[2] - (s in ("b", "d") ? 4.5 : 3.3), text=s, space=:pixel, align=(:center, :baseline), font=:bold, fontsize=12.5)
    end
    f
end
save(joinpath(OUT, "katz_density.png"), f3, px_per_unit=4); save(joinpath(OUT, "katz_density.pdf"), f3)
println("saved to ", OUT)
