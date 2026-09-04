# Water-weakening reduction-coefficient calibration
#
# Reads the ysX_effY_redZ simulation suite in /Volumes/Sinergia/EarthWater and
# plots lithosphere mobility ("Mob" = Vsurf_rms/Vrms, plates_analyse.dat) against
# time for every (yield stress, eruption efficiency) pair, one line per water
# yield-stress reduction fraction "red" (ys_wet_min = (1-red)*ys). A companion
# panel tracks the fraction of run time spent above the mobile-lid threshold
# as a function of "red" — the diagnostic to read off to pick the reduction
# fraction that lands each (ys, eff) pair in the regime Lourenço et al.
# predict for it (stagnant / episodic / mobile lid).
using StagPP, GLMakie, Printf

# ------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------
simroot = "/Volumes/Sinergia/EarthWater"
load_sims = true

# Per-instant mobility threshold above which the lid counts as locally "mobile".
# 0.6 matches the `tectonic_threshold` convention already used by `regime(...;
# tectonic=true)` in src/api.jl.
const Mob_thresh = 0.6
# Time-fraction cutoffs turning that instantaneous threshold into a whole-run
# regime label, matching the (s)/(e)/(m) classification in `regime()`.
const stagnant_cut = 0.08
const mobile_cut = 0.60
# Skip the initial transient before computing time-averaged diagnostics.
tstart = 0.5 # Gyr

# ------------------------------------------------------------------------
# Discover + load simulations
# ------------------------------------------------------------------------
# NOTE: ys is read directly in MPa, eff is the raw directory-name integer
# divided by 100 (eff20 -> 0.2). red encodes the two-decimal-digit fraction
# by which water weakening can reduce the yield stress: red20 means the
# wet minimum achievable yield stress is (1 - 0.20)*ys.
namepat = r"^ys(\d+)_eff(\d+)_red(\d+)$"

function parse_name(name)
    m = match(namepat, name)
    isnothing(m) && return nothing
    ys  = parse(Float64, m.captures[1])        # yield stress [MPa]
    eff = parse(Float64, m.captures[2]) / 100  # eruption efficiency [-]
    red = parse(Float64, m.captures[3]) / 100  # wet yield-stress reduction fraction: ys_wet_min = (1-red)*ys
    return (ys=ys, eff=eff, red=red)
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
        db = load_sim(sroot, name) # rprof=false trips a load_sim bug (haskey(nothing,"r")) in backend.jl, so load rprof too
        push!(sims, (name=name, ys=p.ys, eff=p.eff, red=p.red, db=db))
    end
end

# ------------------------------------------------------------------------
# Mobility regime diagnostics
# ------------------------------------------------------------------------
function mobility_diagnostics(db)
    idxP = data_encoding(db)[3]
    t   = db.platesdata[:, idxP["time"]]
    mob = db.platesdata[:, idxP["Mob"]]
    i0  = something(findfirst(>=(tstart), t), 1)
    (i0 >= length(t)) && (i0 = 1) # run doesn't reach tstart; fall back to the whole series
    tw, mw = @view(t[i0:end]), @view(mob[i0:end])
    mob05 = 0.5(mw[1:end-1] .+ mw[2:end])
    dt    = diff(tw)
    fmob  = sum(dt[mob05.>=Mob_thresh]) / (tw[end] - tw[1])
    mmean = t_avg(mw, tw)
    regime = fmob < stagnant_cut ? "Stagnant lid" : fmob < mobile_cut ? "Episodic lid" : "Mobile lid"
    return (t=t, mob=mob, mean_Mob=mmean, frac_mobile=fmob, regime=regime)
end

diag = Dict(s.name => mobility_diagnostics(s.db) for s in sims)

# Console summary
@printf("%-22s %6s %6s %6s %10s %10s %14s\n", "name", "ys", "eff", "red", "meanMob", "frac_mob", "regime")
for s in sims
    d = diag[s.name]
    @printf("%-22s %6.0f %6.2f %6.2f %10.2e %9.1f%% %14s\n", s.name, s.ys, s.eff, s.red, d.mean_Mob, 100d.frac_mobile, d.regime)
end

# ------------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------------
combos = sort(unique([(s.ys, s.eff) for s in sims]))
redvals = sort(unique([s.red for s in sims]))
redcolors = Dict(r => c for (r, c) in zip(redvals, cpalette(:viridis, length(redvals))))
combocolors = Dict(combo => c for (combo, c) in zip(combos, cpalette(:tab10, length(combos))))

ncols = 3
nrows = cld(length(combos), ncols)

fig = Figure(size=(2000, 750))
Label(fig[0, 1:ncols+1], "Mobility regime calibration against water yield-stress reduction", fontsize=18, tellwidth=false)

for (k, combo) in enumerate(combos)
    ys, eff = combo
    row, col = fldmod1(k, ncols)
    ax = Axis(fig[row, col], yscale=log10, xlabel="Time [Gyr]", ylabel="Mob",
              title="ys = $(Int(ys)) MPa, eff = $(eff)")
    hlines!(ax, [Mob_thresh], color=:gray40, linestyle=:dash, linewidth=1)
    for s in filter(s -> (s.ys, s.eff) == combo, sims)
        d = diag[s.name]
        lines!(ax, d.t, max.(d.mob, 1e-8), color=redcolors[s.red],
               label="red = $(s.red) ($(d.regime), $(round(Int, 100d.frac_mobile))%)")
    end
    axislegend(ax, position=:rb, labelsize=9, framevisible=true, backgroundcolor=(:white, 0.85))
end

# Summary: fraction of time spent mobile vs. reduction fraction, per (ys, eff) — placed to the
# right of the per-simulation grid, spanning its full height.
axsum = Axis(fig[1:nrows, ncols+1], xlabel="Water yield-stress reduction (red)", ylabel="Time spent mobile [%]",
             title="Regime shift with water weakening", xticks=redvals)
hspan!(axsum, 0, 100stagnant_cut, color=(:firebrick, 0.08))
hspan!(axsum, 100stagnant_cut, 100mobile_cut, color=(:goldenrod, 0.08))
hspan!(axsum, 100mobile_cut, 100, color=(:seagreen, 0.08))
text!(axsum, maximum(redvals), 100stagnant_cut/2, text="Stagnant lid", align=(:right, :center), fontsize=11, color=:firebrick)
text!(axsum, maximum(redvals), 50(stagnant_cut+mobile_cut), text="Episodic lid", align=(:right, :center), fontsize=11, color=:goldenrod4)
text!(axsum, maximum(redvals), 100(mobile_cut+1)/2, text="Mobile lid", align=(:right, :center), fontsize=11, color=:seagreen)
for combo in combos
    ys, eff = combo
    csims = sort(filter(s -> (s.ys, s.eff) == combo, sims), by=s -> s.red)
    rv = [s.red for s in csims]
    fv = [100diag[s.name].frac_mobile for s in csims]
    scatterlines!(axsum, rv, fv, color=combocolors[combo], label="ys=$(Int(ys)), eff=$(eff)")
end
ylims!(axsum, 0, 100)
axislegend(axsum, position=:lt, labelsize=10, framevisible=true, backgroundcolor=(:white, 0.85), nbanks=1)
colsize!(fig.layout, ncols+1, Relative(0.30))

display(fig)
savepath = joinpath(@__DIR__, "imgs", "red_mobility_regimes.png")
mkpath(dirname(savepath))
# save(savepath, fig)
# println("Saved figure to $savepath")
