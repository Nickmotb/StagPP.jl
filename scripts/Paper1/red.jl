# Tectonic regime survey of the water-friction sweep
#
# Reads the ysX_effY_muwZ suite in /Volumes/Sinergia/EarthWater, where water weakens the
# brittle yield stress by ramping the friction coefficient from its dry value (0.6, Byerlee
# 1978) down to muw where the lid is hydrated. Three panels:
#
#   1. Mobility–plateness plane — one point per run, colored by surface heat flux. The
#      threshold-free view: stagnant and plutonic-squishy lids both sit at low mobility and
#      separate on plateness, since squishy localizes deformation into weak zones between
#      blocks while a stagnant lid deforms diffusely.
#   2. Regime heatmaps — the same runs classified, rows = friction set, columns = the swept
#      parameter (yield stress in one slice, eruption efficiency in the other).
#   3. Regime raster — one strip per run through time, mobile intervals highlighted. Keeps
#      the temporal information the two time-averaged panels above discard, which is what
#      separates a bursty episodic lid from a persistently deforming squishy one.
using StagPP, GLMakie, Printf, Statistics

# ------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------
simroot = "/Volumes/Sinergia/EarthWater"
load_sims = false

# Per-instant mobility threshold above which the lid counts as locally "mobile".
# 0.6 matches the `tectonic_threshold` convention already used by `regime(...;
# tectonic=true)` in src/api.jl.
const Mob_thresh = 0.6
# Time-fraction cutoffs turning that instantaneous threshold into a whole-run
# regime label, matching the (s)/(e)/(m) classification in `regime()`.
const stagnant_cut = 0.08
const mobile_cut = 0.60
# Plateness above which a low-mobility run counts as plutonic-squishy rather than stagnant.
# PROVISIONAL: the archived ductile-reduction runs sat at P = 0.14-0.28 with nothing above,
# so this is placed clear of that floor rather than calibrated. Re-set it once the friction
# runs show where the population actually splits — panel 1 is the threshold-free check.
const squishy_P_cut = 0.40
# Surface-area fraction accommodating 80% of the surface deformation for isoviscous
# convection — the reference Tackley (2000) normalizes plateness against, P = 1 - f80/f80_iso.
const f80_iso = 0.6
# Skip the initial transient before computing time-averaged diagnostics.
tstart = 0.5 # Gyr

const REGIMES = ["Stagnant", "Squishy", "Episodic", "Mobile"]
const REGIME_COLORS = [:grey60, :mediumpurple3, :goldenrod2, :seagreen]

# ------------------------------------------------------------------------
# Discover + load simulations
# ------------------------------------------------------------------------
# ys is MPa; eff is the raw directory integer / 100 (eff20 -> 0.2 erupted, 0.8 intruded);
# muw is the wet friction coefficient x 100 (muw05 -> 0.05).
namepat = r"^ys(\d+)_eff(\d+)_muw(\d+)$"

function parse_name(name)
    m = match(namepat, name)
    isnothing(m) && return nothing
    ys  = parse(Float64, m.captures[1])        # yield stress [MPa]
    eff = parse(Float64, m.captures[2]) / 100  # eruption efficiency [-]
    muw = parse(Float64, m.captures[3]) / 100  # wet friction coefficient [-]
    return (ys=ys, eff=eff, muw=muw)
end

simnames = sort(filter(n -> !isnothing(parse_name(n)), readdir(simroot)))

if load_sims
    sims = NamedTuple[]
    for name in simnames
        p = parse_name(name)
        sroot = joinpath(simroot, name, "+op", "EW")
        if !isfile(sroot * "_plates_analyse.dat")
            @warn "Skipping $name: no plates_analyse.dat found"
            continue
        end
        db = load_sim(sroot, name) # rprof=false trips a load_sim bug (haskey(nothing,"r")) in backend.jl
        push!(sims, (name=name, ys=p.ys, eff=p.eff, muw=p.muw, db=db))
    end
end

# ------------------------------------------------------------------------
# Regime diagnostics
# ------------------------------------------------------------------------
# First index at or after tstart, falling back to the whole series for runs that
# have not yet reached it.
function window_start(t)
    i0 = something(findfirst(>=(tstart), t), 1)
    return i0 >= length(t) ? 1 : i0
end

function regime_diagnostics(db)
    idxT, _, idxP = data_encoding(db)
    # --- Plates block: mobility and plateness
    t   = db.platesdata[:, idxP["time"]]
    mob = db.platesdata[:, idxP["Mob"]]
    f80 = db.platesdata[:, idxP["f80"]]
    i0  = window_start(t)
    tw, mw, fw = @view(t[i0:end]), @view(mob[i0:end]), @view(f80[i0:end])
    mob05 = 0.5(mw[1:end-1] .+ mw[2:end])
    fmob  = sum(diff(tw)[mob05.>=Mob_thresh]) / (tw[end] - tw[1])
    P     = 1 - t_avg(fw, tw)/f80_iso
    # --- Time block: surface heat flux, on its own (denser) sampling
    tt = db.timedata[:, idxT["time"]]
    j0 = window_start(tt)
    Ftop = 1e3t_avg(@view(db.timedata[j0:end, idxT["F_top"]]), @view(tt[j0:end])) # W/m² -> mW/m²
    # Mobility alone cannot separate stagnant from squishy, so plateness breaks the tie
    # among the low-mobility runs.
    code = fmob >= mobile_cut   ? 4 :
           fmob >= stagnant_cut ? 3 :
           P    >= squishy_P_cut ? 2 : 1
    # The median tracks the persistent baseline, the mean is dominated by overturn spikes:
    # a squishy lid sits above the stagnant floor in both, an episodic one only in the mean.
    return (t=t, mob=mob, mean_Mob=t_avg(mw, tw), med_Mob=median(mw),
            frac_mobile=fmob, plateness=P, Ftop=Ftop, code=code, regime=REGIMES[code])
end

diag = Dict(s.name => regime_diagnostics(s.db) for s in sims)

# Console summary
@printf("%-20s %5s %5s %5s %10s %10s %9s %7s %8s %10s\n",
        "name", "ys", "eff", "muw", "meanMob", "medMob", "frac_mob", "P", "Ftop", "regime")
for s in sims
    d = diag[s.name]
    @printf("%-20s %5.0f %5.2f %5.2f %10.2e %10.2e %8.1f%% %7.3f %8.1f %10s\n",
            s.name, s.ys, s.eff, s.muw, d.mean_Mob, d.med_Mob, 100d.frac_mobile, d.plateness, d.Ftop, d.regime)
end

# ------------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------------
muvals = sort(unique([s.muw for s in sims]))
markershapes = [:circle, :rect, :utriangle, :diamond, :pentagon]
mumarker = Dict(mu => markershapes[i] for (i, mu) in enumerate(muvals))

# The suite is two slices through (ys, eff, muw) sharing one corner: yield stress is swept at
# the most-sampled eruption efficiency, and eruption efficiency at the most-sampled yield stress.
most_common(v) = (u = unique(v); u[argmax([count(==(x), v) for x in u])])
effA = most_common([s.eff for s in sims])
ysB  = most_common([s.ys for s in sims])

fig = Figure(size=(1900, 1250))
Label(fig[1, 1:4], "Tectonic regimes of the water-friction sweep", fontsize=19, tellwidth=false)

# --- Panel 1: mobility–plateness plane (no thresholds applied to the data itself)
axsc = Axis(fig[2, 1], xscale=log10, title="Mobility–plateness plane",
            xlabel="Mob  (whisker: median → time-mean)", ylabel="Plateness  P = 1 - f80/$(f80_iso)")
vlines!(axsc, [Mob_thresh], color=:gray40, linestyle=:dash, linewidth=1)
hlines!(axsc, [squishy_P_cut], color=:gray40, linestyle=:dash, linewidth=1)
Fvals = [diag[s.name].Ftop for s in sims]
crange = (minimum(Fvals), maximum(Fvals))
for s in sims
    d = diag[s.name]
    lines!(axsc, [d.med_Mob, d.mean_Mob], [d.plateness, d.plateness], color=(:gray55, 0.6), linewidth=1)
    scatter!(axsc, [d.mean_Mob], [d.plateness], color=[d.Ftop], colormap=:inferno, colorrange=crange,
             marker=mumarker[s.muw], markersize=14, strokewidth=0.6, strokecolor=:black)
end
text!(axsc, 0.02, 0.02, text="Stagnant", space=:relative, fontsize=11, color=:grey40)
text!(axsc, 0.02, 0.96, text="Plutonic-squishy", space=:relative, fontsize=11, color=:mediumpurple4)
text!(axsc, 0.98, 0.86, text="Mobile lid", space=:relative, align=(:right, :baseline), fontsize=11, color=:seagreen)
axislegend(axsc,
           [MarkerElement(marker=mumarker[mu], color=:gray40, strokecolor=:black, strokewidth=0.6, markersize=12) for mu in muvals],
           ["μ_wet = $(mu)" for mu in muvals],
           position=(0.02, 0.62), labelsize=10, framevisible=true, backgroundcolor=(:white, 0.85))
Colorbar(fig[2, 2], limits=crange, colormap=:inferno, label="Mean surface heat flux [mW/m²]",
         width=14, ticklabelsize=9, labelsize=10)

# --- Panel 2: regime heatmaps, one per slice
function regime_heatmap!(pos, subset, xkey, xlabel, title)
    xvals = sort(unique([getproperty(s, xkey) for s in subset]))
    codes = fill(NaN, length(xvals), length(muvals)) # heatmap! wants [x, y]
    ax = Axis(fig[pos...], title=title, xlabel=xlabel, ylabel="Wet friction μ_wet",
              xticks=(1:length(xvals), [x == round(x) ? string(Int(x)) : string(x) for x in xvals]),
              yticks=(1:length(muvals), string.(muvals)))
    for (i, mu) in enumerate(muvals), (j, xv) in enumerate(xvals)
        k = findfirst(s -> s.muw == mu && getproperty(s, xkey) == xv, subset)
        isnothing(k) && continue
        d = diag[subset[k].name]
        codes[j, i] = d.code
        text!(ax, j, i, text=@sprintf("%.0f%%\nP %.2f", 100d.frac_mobile, d.plateness),
              align=(:center, :center), fontsize=9, color=:black)
    end
    heatmap!(ax, 1:length(xvals), 1:length(muvals), codes,
             colormap=cgrad(REGIME_COLORS, categorical=true), colorrange=(0.5, 4.5))
    translate!(ax.scene.plots[end], 0, 0, -1) # keep the cell text on top
    return ax
end

regime_heatmap!((2, 3), filter(s -> s.eff == effA, sims), :ys,
                "Ductile yield stress [MPa]", "Regime vs yield stress  (eff = $(effA))")
regime_heatmap!((2, 4), filter(s -> s.ys == ysB, sims), :eff,
                "Eruption efficiency [-]", "Regime vs eruption efficiency  (ys = $(Int(ysB)) MPa)")
Legend(fig[3, 1:4],
       [PolyElement(color=c) for c in REGIME_COLORS], REGIMES,
       orientation=:horizontal, framevisible=false, labelsize=12, tellheight=true)

# --- Panel 3: regime through time, one strip per run
axr = Axis(fig[4, 1:4], title="Mobile intervals through time  (Mob ≥ $(Mob_thresh))",
           xlabel="Time [Gyr]", yticks=(1:length(sims), [s.name for s in sims]),
           yticklabelsize=8)
vlines!(axr, [tstart], color=:gray30, linestyle=:dash, linewidth=1)
for (i, s) in enumerate(sims)
    d = diag[s.name]
    lines!(axr, [d.t[1], d.t[end]], [i, i], color=(:grey91, 1.0), linewidth=7)
    segs = Point2f[]
    for j in 1:length(d.t)-1
        0.5(d.mob[j] + d.mob[j+1]) >= Mob_thresh || continue
        push!(segs, Point2f(d.t[j], i)); push!(segs, Point2f(d.t[j+1], i))
    end
    isempty(segs) || linesegments!(axr, segs, color=REGIME_COLORS[d.code], linewidth=7)
end
ylims!(axr, 0.3, length(sims) + 0.7)

rowsize!(fig.layout, 4, Relative(0.52))
display(fig)
savepath = joinpath(@__DIR__, "imgs", "red_mobility_regimes.png")
mkpath(dirname(savepath))
# save(savepath, fig)
# println("Saved figure to $savepath")
