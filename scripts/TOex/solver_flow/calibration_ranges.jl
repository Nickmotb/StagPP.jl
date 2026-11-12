# Calibration ranges of every model law in P-T (figs/calibration.png for TOex_solver_flow.tex).
# Boxes are the source data ranges; dashed edges mark a bound not yet checked against the paper.
using CairoMakie

const OUT = joinpath(@__DIR__, "figs")
const PLO = 0.01                                                           # log axis floor (GPa)
sectors = [(PLO, 7, 500, 2000), (7, 25, 700, 4000), (25, 135, 700, 4000)]   # model tables: P (GPa), T (K)

rect(P1, P2, T1, T2) = Point2f[(T1, P1), (T2, P1), (T2, P2), (T1, P2)]
function box!(ax, P1, P2, T1, T2, col; dash=(false, false, false), alpha=0.12, lw=2)
    P1 = max(P1, PLO)
    poly!(ax, rect(P1, P2, T1, T2), color=(col, alpha), strokewidth=0)
    ls(d) = d ? :dash : :solid
    lines!(ax, [T1, T2], [P1, P1], color=col, linewidth=lw)
    lines!(ax, [T1, T1], [P1, P2], color=col, linewidth=lw)
    lines!(ax, [T1, T2], [P2, P2], color=col, linewidth=lw, linestyle=ls(dash[1]))   # high-P edge
    lines!(ax, [T2, T2], [P1, P2], color=col, linewidth=lw, linestyle=ls(dash[2]))   # high-T edge
end
lab!(ax, T, P, s, col; al=(:left, :center), fs=13) = text!(ax, T, P, text=s, color=col, fontsize=fs, align=al, font=:bold)
function frame(pos, title)
    ax = Axis(pos, title=title, xlabel="Temperature (K)", ylabel="Pressure (GPa)", yreversed=true, yscale=log10,
              yticks=([0.01, 0.1, 1, 10, 100], ["0.01", "0.1", "1", "10", "100"]), titlesize=17)
    for (P1, P2, T1, T2) in sectors
        poly!(ax, rect(P1, P2, T1, T2), color=(:gray, 0.10), strokecolor=(:gray, 0.6), strokewidth=1)
    end
    limits!(ax, 250, 4250, 150, PLO)
    return ax
end

f = Figure(size=(1900, 980), fontsize=15)

# (a) melt, fluid, redox and phase equilibria
a = frame(f[1, 1], "(a) Melt, fluid, redox and phase-equilibrium laws")
box!(a, PLO, 10, 273, 2573, "#56B4E9")                                     # ZD09 calibrated
lines!(a, [273, 4200], [50, 50], color="#2b7fb0", linewidth=2, linestyle=:dash)
box!(a, PLO, 10, 273, 2000, "#E69F00"; alpha=0.10)                         # PS94
box!(a, 5, 60, 773, 1873, "#7E57C2")                                       # DHMS kappa map
box!(a, 2.5, 11, 1373, 1873, "#CC79A7"; alpha=0.25)                        # Stagno & Frost
box!(a, 0.7, 3, 1623, 1723, "#0072B2"; alpha=0.45)                         # H12, ~1400 C
box!(a, PLO, 2, 933, 1700, "#009E73"; dash=(false, true, false), alpha=0.25)   # Sun & Yao
box!(a, PLO, 150, 1523, 1773, "#D55E00"; dash=(true, false, false), alpha=0.08, lw=1.5)  # H22 T band
lab!(a, 2600, 0.03, "ZD09 fluid EoS: calibrated\n(used above 773 K)", "#2b7fb0")
lab!(a, 2600, 20, "ZD09 extrapolated", "#2b7fb0")
lab!(a, 2600, 55, "above 50 GPa: ideal pure-water fluid", "#2b7fb0")
lab!(a, 290, 0.02, "PS94 pure H₂O, CO₂\n(hybrid, below 773 K)", "#B07800"; al=(:left, :top))
lab!(a, 790, 30, "DHMS κ map\n(K&O 2006,\nNishi et al. 2014)", "#7E57C2")
lab!(a, 1890, 6, "Melt CO₂, Stagno & Frost (2010)", "#A0507F")
lab!(a, 1740, 1.4, "Melt H₂, Hirschmann\net al. (2012)", "#0072B2")
lab!(a, 950, 0.5, "Melt H₂O\nSun & Yao (2026)", "#00775A")
lab!(a, 1648, 120, "Melt Fe³⁺, Hirschmann (2022):\nT of the experiments, P to check", "#D55E00"; al=(:center, :center), fs=12)
text!(a, 3300, 140, text="sb24 (MAGEMin): whole mantle", fontsize=13, align=(:center, :bottom), color=:gray30)

# (b) NAM water laws
b = frame(f[1, 2], "(b) Water in nominally anhydrous minerals")
nams = [  # label, P1, P2, T1, T2, label position (T, P, align)
    ("ol", 0.01, 14.2, 1200, 2400, (2380, 0.02, (:right, :top))),
    ("opx En", 0.01, 14, 500, 1500, (520, 0.02, (:left, :top))),
    ("opx Al", 0.01, 14, 1073, 1500, (1090, 0.2, (:left, :top))),
    ("di, jd", 0.01, 10, 873, 973, (880, 2.0, (:left, :center))),
    ("gt, maj", 0.01, 10, 700, 1500, (720, 0.6, (:left, :center))),
    ("coe, qtz", 0.01, 10, 800, 1073, (810, 5.0, (:left, :center))),
    ("wad", 14, 20, 1473, 2173, (2190, 16.5, (:left, :center))),
    ("ring (ak = ring/21)", 19, 23, 1723, 2000, (2190, 21, (:left, :center))),
    ("hpcpx", 8, 13.4, 1373, 1673, (1690, 10.5, (:left, :center))),
    ("st", 0.01, 60, 800, 1500, (1150, 40, (:center, :center))),
    ("CaCl₂ silica", 25, 100, 2900, 3500, (3200, 50, (:center, :center))),
    ("α-PbO₂ silica", 128, 144, 3650, 4100, (3640, 120, (:right, :center))),
    ("ppv", 90, 135, 1000, 4000, (1100, 110, (:left, :center))),
]
pal = ["#E69F00", "#56B4E9", "#009E73", "#D55E00", "#0072B2", "#CC79A7", "#7E57C2", "#B07800", "#00775A", "#A0507F", "#2b7fb0", "#666666", "#C1272D"]
for (k, (l, P1, P2, T1, T2, (tx, py, al))) in enumerate(nams)
    box!(b, P1, P2, T1, T2, pal[k]; alpha=0.10, lw=1.8)
    lab!(b, tx, py, l, pal[k]; al=al, fs=12)
end
text!(b, 3300, 0.012, text="bridgmanite: Al₂O₃ law (Litasov et al. 2003)\nferropericlase, CF, Ca-pv: constants", fontsize=12, align=(:center, :top), color=:gray30)

Label(f[2, 1:2], "Grey: the model's table sectors (upper mantle 0–7 GPa, transition zone 7–25 GPa, lower mantle 25–135 GPa). Dashed edges: a bound not yet checked against the source.",
      fontsize=13, tellwidth=false)
save(joinpath(OUT, "calibration.png"), f, px_per_unit=2)
println("wrote ", joinpath(OUT, "calibration.png"))
