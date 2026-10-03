# Lunar water post-processing: sweep survey + per-run time series
#
# Reads one sweep folder (e.g. /Volumes/Sinergia2/MoonWater/IDgood), loads every case in it and
# plots either:
#
#   MODE = :heatmap  one cell per run on the two swept parameters, for the four quantities the
#                    sample data constrain: erupted volume, lava H2O, bulk mantle H2O and the
#                    outer-layer H2O. The axes are discovered by diffing the pf files, so the
#                    same script works for any sweep, and the cells carry the number as text so
#                    a near-miss reads as such instead of being lost in the colour scale.
#                    A second figure scores each run against the samples: how far the time series
#                    sits from the measured ranges, in log10 units (0 = inside the range).
#   MODE = :normal   one five-panel figure reading the runs in SELECT against the samples:
#                    (a) erupted magma H2O against every sample estimate, (b) when the mare was
#                    emplaced, (c) magma H2O against the melting-zone source it came from,
#                    (d) the match target by target, (e) the CE5-age source fit against the late
#                    source, for every run of every sweep. (a), (b), (c) and (d) overlay every selected run.
#   MODE = :water    one three-panel figure on the hydration state of the runs in SELECT:
#                    (a) water and storage capacity against depth at five times, (c) what each
#                    shell holds against the bulk, (b) saturation of the capacity in depth and
#                    time for the first run, with the melting region and the shells of (c)
#                    outlined. (a) and (c) overlay every selected run.
#   MODE = :fit      one plane per sample stage (Apollo, CE6, CE5): log distance of the window lava
#                    from the centre of the stage's melt range on x, of the window source from the
#                    centre of its mantle range on y; the ranges drawn as the box. Every run of
#                    every sweep as a grey dot, SELECT in colour.
#   MODE = :sector   one four-panel figure reading Moon_sector.dat for the sample-site sectors of
#                    every run in SELECT: (a) lava H2O by sample sector against the sample boxes,
#                    (b) the same for the melt source, (c) cumulative erupted volume near vs far
#                    sectors, (d) window means, far/near ratio and last-eruption age per case.
using StagPP, GLMakie, DelimitedFiles, Statistics, Printf

LOCKED = false   # true: reuse data loaded earlier in this Julia session (run with include in a REPL)

# In-session cache, keyed per data source; survives repeated `include`, cleared by restarting Julia.
if !@isdefined(LW_CACHE)
    const LW_CACHE = Dict{Any,Any}()
end
LW_CACHE_HITS = 0; LW_CACHE_MISS = 0   # reset every run

# Reuses the stored value when LOCKED and the source files are unchanged; else recomputes and stores.
function cached(f, key, files...)
    mt = maximum(mtime, files)
    if LOCKED && haskey(LW_CACHE, key) && LW_CACHE[key][2] == mt
        global LW_CACHE_HITS += 1
        return LW_CACHE[key][1]
    end
    global LW_CACHE_MISS += 1
    v = f()
    LW_CACHE[key] = (v, mt)
    return v
end

# ------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------
sweep = "/Volumes/Sinergia2/MoonWater/IDgood2"
MODE  = :sector          # :heatmap | :normal | :water | :fit | :sector
# Case folder names to use, in order. Empty vector = every case in the sweep. Also restricts the
# heatmap, which is how you drop a subset that would share cells (e.g. capped vs uncapped runs).
SELECT = String["r09_D008"]
# Heatmap axes, as pf parameter names. Empty = whichever parameters differ between the cases.
AXES = String[]

tstart = 0.6              # Gyr; skips the initial transient everywhere
outer_default = 300e3     # m, outer-layer thickness if the run wrote no water_layer_depth
savein = "r09_D008"               # basename to save the figure, "" = display only
INTERP = false            # heatmap colour: false = one flat colour per cell, true = interpolated
STATS  = true             # second figure scoring the runs against the sample ranges
obs_halfwidth = 0.15      # Gyr; model window averaged around each sample's eruption age
obs_volume = 1e7          # km3 of mare basalt the volume score aims at
TRADE_SKIP = ["drydeep"]  # sweeps left out of the grey cloud in (e)
TRADE_EXTRA = ["/Volumes/Sinergia2/Testing"]   # folders outside the study also read into the cloud of (e)
MARE_END   = 1.2          # Ga, youngest mare units: eruption should fade out by here
VOL_BAND   = 0.30         # half-width of the band drawn around it in (b), as a fraction
FIT_RANGE  = :widest      # :fit stage ranges: :widest = min of mins..max of maxs, :central = median min..median max

SECTORS = ["CE5 site", "PKT core", "Apollo band", "CE6 site", "SPA", "rest"]   # names in file order
SECTOR_OF = Dict("Apollo" => 3, "CE6" => 4, "CE5" => 1)   # which sector each sample window is read from
NEAR_SECTORS = [1, 2, 3]; FAR_SECTORS = [4, 5]            # :sector mode: volume grouping for panel (c)

# erupta reaches StagPP in kg (mantle-mass fraction x Mmantle), so km3 needs 1/(rho*1e9).
const KM3_PER_KG = 1 / (3300 * 1e9)
const WT2UGG = 1e4        # wt% -> ug/g
const OBS_CSV = "/Users/nickmotb/Desktop/LunarWater/Table2.csv"
const CMAP = :vik100
const MISFIT_CMAP = cgrad(:lajolla, rev=true)   # one-sided: pale = in the data, dark = the worst run

# Sample windows quoted in the text, as model time (t = 4.5 - age).
const WINDOWS = ("Apollo" => (0.6, 1.4), "CE6" => (1.5, 1.9), "CE5" => (2.3, 2.7))

# ------------------------------------------------------------------------
# Sweep discovery
# ------------------------------------------------------------------------
cases(sw) = sort([d for d in readdir(sw) if isfile(joinpath(sw, d, "+op", "Moon_time.dat"))])

# Eruption time series straight from Moon_time.dat: the fast path for cases not loaded by StagPP.
function read_timedat(path)
    h = split(readline(path))
    A = readdlm(path, skipstart=1)
    ix = Dict(h[i] => i for i in eachindex(h))
    # overflowed columns lose the E of the exponent ("1.07+238"), so parse rather than convert
    num(x) = x isa Number ? Float64(x) :
             something(tryparse(Float64, String(x)),
                       tryparse(Float64, replace(String(x), r"(\d)([+-]\d+)$" => s"\1e\2")), NaN)
    col(k) = [num(A[i, ix[k]]) for i in axes(A, 1)]
    return (t = col("time") ./ 3.1557e16,                      # s -> Gyr
            dE = vcat(0.0, max.(diff(col("erupta")), 0.0)),
            lava = col("cH2Oerupt") .* WT2UGG)
end

# Same idea for Moon_sector.dat: n sectors, each with its own 8 columns "<name>_<k>".
function read_sectordat(path)
    h = split(readline(path))
    A = readdlm(path, skipstart=1)
    ix = Dict(h[i] => i for i in eachindex(h))
    num(x) = x isa Number ? Float64(x) :
             something(tryparse(Float64, String(x)),
                       tryparse(Float64, replace(String(x), r"(\d)([+-]\d+)$" => s"\1e\2")), NaN)
    col(k) = [num(A[i, ix[k]]) for i in axes(A, 1)]
    n = div(length(h) - 2, 8)
    E = [col("erupta_$k") for k in 1:n]
    return (t = col("time") ./ 3.1557e16, n = n,
            dE = [vcat(0.0, max.(diff(E[k]), 0.0)) for k in 1:n],
            lava = [col("cH2Oerupt_$k") .* WT2UGG for k in 1:n],
            src = [col("cH2Osrc_$k") .* WT2UGG for k in 1:n],
            E = E, crthick = [col("crthick_$k") ./ 1e3 for k in 1:n])
end

# Melt-weighted source water per rprof frame, as panel (c) defines it, straight from the file:
#  shells with melt at depth <= 420 km, weighted by fmelt x shell mass (rho r^2 stands in for dM).
function read_rprof_src(path)
    t = Float64[]; src = Float64[]
    ix = Dict{String,Int}(); rows = Vector{Vector{SubString{String}}}()
    flush!(tt) = begin
        isempty(rows) && return
        r = [parse(Float64, w[ix["r"]]) for w in rows]; zc = (maximum(r) .- r) ./ 1e3
        f = [parse(Float64, w[ix["fmeltmean"]]) for w in rows]
        W = [(parse(Float64, w[ix["boundH2O"]]) + parse(Float64, w[ix["freeH2O"]])) * WT2UGG for w in rows]
        m = [parse(Float64, w[ix["rhomean"]]) for w in rows] .* (r .+ 350e3) .^ 2
        k = (zc .<= 420) .& (f .> 0)
        push!(t, tt); push!(src, any(k) ? sum(W[k] .* f[k] .* m[k]) / sum(f[k] .* m[k]) : NaN)
        empty!(rows)
    end
    tcur = NaN
    for line in eachline(path)
        if isempty(ix)
            h = split(line); ix = Dict(h[i] => i for i in eachindex(h)); continue
        end
        if startswith(line, '*')
            flush!(tcur)
            tcur = parse(Float64, match(r"time\s*=\s*([\dE+.-]+)", line).captures[1]) / 3.1557e16
        elseif !isempty(strip(line))
            push!(rows, split(line))
        end
    end
    flush!(tcur)
    return (t = t, src = src)
end

# Every assignment in a pf, comments stripped. Used to find what the sweep varies.
function parse_pf(path)
    out = Dict{String,String}()
    for raw in eachline(path)
        line = split(raw, '!')[1]
        for m in eachmatch(r"([A-Za-z_]\w*)\s*(\([^)]*\))?\s*=\s*([^=]*?)(?=\s+[A-Za-z_]\w*\s*(\([^)]*\))?\s*=|$)", line)
            out[lowercase(m.captures[1])] = strip(m.captures[3])
        end
    end
    return out
end

# Parameters whose value differs between cases are the sweep axes; everything else is shared.
function sweep_axes(sw, cs)
    pfs = Dict(c => parse_pf(joinpath(sw, c, "pf")) for c in cs)
    keys_all = union((Set(keys(p)) for p in values(pfs))...)
    axes = Dict{String,Vector{String}}()
    for k in keys_all
        vals = [String(get(pfs[c], k, "(unset)")) for c in cs]
        length(unique(vals)) > 1 || continue
        axes[k] = vals
    end
    return axes, pfs
end

# Axis ticks keep the pf spelling; order numerically when every value is a number, since a
# parameter missing from some cases ("(unset)") makes the column non-numeric.
function ordered_unique(vals)
    u = unique(vals)
    k = [something(tryparse(Float64, v), NaN) for v in u]
    return all(isfinite, k) ? u[sortperm(k)] : sort(u)
end

# ------------------------------------------------------------------------
# Sample constraints
# ------------------------------------------------------------------------
# Table2 rows come in min/max pairs: [group; t; min; max] for melt and for mantle.
function observations()
    csv = readdlm(OBS_CSV, ','); replace!(csv, "" => NaN)
    n = Int(0.5(size(csv, 1) - 2))
    for tag in ("AP11", "AP12", "AP15", "AP17", "Ap7")
        replace!(csv, tag => "Apollo")
    end
    # Apollo suites and the two Chang'e missions keep their identity; every other entry is a
    # single "Rest" colour, so the legend stays readable.
    grps = string.(csv[2:2:end-1, 1])
    grps = [g in ("Apollo", "CE5", "CE6") ? g : "Rest" for g in grps]
    tmat = zeros(5, n)
    tmat[1, :] .= 4.5 .- csv[2:2:end-1, 2]
    tmat[2, :] .= csv[2:2:end-1, 4];   tmat[3, :] .= csv[3:2:end-1, 4]  # primary melt
    tmat[4, :] .= csv[2:2:end-1, 5];   tmat[5, :] .= csv[3:2:end-1, 5]  # mantle source
    return vcat(permutedims(grps), tmat)
end

# ------------------------------------------------------------------------
# Per-run diagnostics
# ------------------------------------------------------------------------
# Thickness of the initially water-enriched outer layer, as StagYY wrote it out.
function outer_layer_depth(db::DataBlock)
    pf = joinpath(db.metadata.outdir, db.metadata.sroot * "_parameters.dat")
    m = isfile(pf) ? match(r"WATER_LAYER_DEPTH\s*=\s*([0-9.eE+-]+)", read(pf, String)) : nothing
    return isnothing(m) ? outer_default : parse(Float64, m.captures[1])
end

wmean(v, w) = (s = sum(w); s > 0 ? sum(v .* w) / s : NaN)

function diagnose(db::DataBlock)
    iT, iR, _ = data_encoding(db.timeheader, db.rprofheader, db.platesheader)
    T = db.timedata
    t, er = T[:, iT["time"]], T[:, iT["erupta"]]
    lava, bulk = T[:, iT["cH2Oerupt"]] .* WT2UGG, T[:, iT["cH2O_mean"]] .* WT2UGG
    i0 = findlast(t .<= tstart)
    dE = vcat(0.0, diff(er))
    keep = (t .> tstart) .& (dE .> 0) .& (lava .> 1e-6)
    # Erupted volume weights the lava composition: a magma counts for as much as erupted of it.
    lava_vw = wmean(lava[keep], dE[keep])
    win = Dict(name => let k = (t .> a) .& (t .<= b) .& (dE .> 0) .& (lava .> 1e-6)
                   wmean(lava[k], dE[k])
               end for (name, (a, b)) in WINDOWS)
    r = db.rprofdata[:, 1, iR["r"]]
    top = findall(>=(maximum(r) - outer_layer_depth(db)), r)
    W, dM = db.rprofdata[top, :, iR["Water"]], db.rprofdata[top, :, iR["dM"]]
    outer = [wmean(W[:, j], dM[:, j]) for j in eachindex(db.rproftime)] .* WT2UGG
    jr = findall(>=(tstart), db.rproftime)
    return (erupt_km3 = (er[end] - er[i0]) * KM3_PER_KG,
            lava = lava_vw, windows = win,
            bulk = mean(bulk[t .>= tstart]), bulk_end = bulk[end],
            outer = mean(outer[jr]), outer_end = outer[end],
            tend = t[end],
            series = (t = t, dE = dE, lava = lava, bulk = bulk))   # kept for the sample scores
end

cs = isempty(SELECT) ? cases(sweep) : SELECT
isempty(cs) && error("no cases with +op/Moon_time.dat found in $sweep")
# :sector reads Moon_sector.dat directly; skip the StagPP load for cases that only have that file.
if MODE != :sector
    simfiles(c) = (joinpath(sweep, c, "+op", "Moon_time.dat"), joinpath(sweep, c, "+op", "Moon_rprof.dat"))
    sims = Dict(c => cached(() -> load_sim(joinpath(sweep, c, "+op", "Moon"), c, plates=false),
                            ("sim", sweep, c), simfiles(c)...) for c in cs)
    diag = Dict(c => cached(() -> diagnose(sims[c]), ("diag", sweep, c), simfiles(c)...) for c in cs)

    for c in cs
        d = diag[c]
        @printf("%-14s  V=%8.3e km3  lava=%7.1f  bulk=%6.2f  outer=%6.2f  ug/g   (t_end %.2f Gyr)\n",
                c, d.erupt_km3, d.lava, d.bulk, d.outer, d.tend)
    end
end

# ------------------------------------------------------------------------
# Heatmap mode
# ------------------------------------------------------------------------
# Text on top of the colour, so pick black or white from the luminance of the cell colour.
function label_color(v, lo, hi, cmap)
    c = (cmap isa Symbol ? cgrad(cmap) : cmap)[clamp((v - lo) / max(hi - lo, eps()), 0, 1)]
    return 0.299c.r + 0.587c.g + 0.114c.b > 0.55 ? :black : :white
end

function sweep_heatmap!(fig, pos, xs, ys, xname, yname, Z, title, fmt;
                        logscale=false, cmap=CMAP, crange=nothing, interpolate=INTERP)
    C = logscale ? log10.(Z) : Z
    lo, hi = isnothing(crange) ? extrema(filter(isfinite, C)) : crange
    ax = Axis(fig[pos...], title=title, xlabel=xname, ylabel=yname,
              xticks=(1:length(xs), string.(xs)), yticks=(1:length(ys), string.(ys)),
              titlesize=11, xlabelsize=10, ylabelsize=10, xticklabelsize=9, yticklabelsize=9)
    hm = heatmap!(ax, 1:length(xs), 1:length(ys), C, colormap=cmap, colorrange=(lo, hi),
                  interpolate=interpolate)
    translate!(hm, 0, 0, -1)   # cell labels stay on top of the colour
    for i in eachindex(xs), j in eachindex(ys)
        isfinite(Z[i, j]) || continue
        text!(ax, i, j, text=fmt(Z[i, j]), align=(:center, :center), fontsize=9,
              color=label_color(C[i, j], lo, hi, cmap))
    end
    Colorbar(fig[pos[1], pos[2] + 1], limits=(lo, hi), colormap=cmap, width=10,
             ticklabelsize=8, label=logscale ? "log10" : "")
    return ax
end

# ------------------------------------------------------------------------
# Sample scores: how far a run sits from the measured ranges, in log10 units
# ------------------------------------------------------------------------
# 0 while the model is inside [lo, hi]; otherwise the distance to the nearer edge, signed so a
# positive score is a model that is too wet.
function log_offset(m, lo, hi)
    (isfinite(m) && m > 0) || return NaN
    l, h = min(lo, hi), max(lo, hi)
    isfinite(l) && m < l && return log10(m / l)
    isfinite(h) && m > h && return log10(m / h)
    (isfinite(l) || isfinite(h)) ? 0.0 : NaN
end

# Lava is sampled the way the observations are: the volume-weighted mean of the magma erupted
# around that sample's age, so a run that erupts nothing there scores NaN rather than 0.
function lava_at(s, t0, hw)
    k = (s.t .> t0 - hw) .& (s.t .<= t0 + hw) .& (s.dE .> 0) .& (s.lava .> 1e-6)
    any(k) ? wmean(s.lava[k], s.dE[k]) : NaN
end

bulk_at(s, t0) = (i = findlast(s.t .<= t0); isnothing(i) ? NaN : s.bulk[i])

# Per-run scores against every sample in Table2: RMS distance (how far off) and mean signed
# distance (which way), plus the volume score against the mare estimate.
function scores(d, obs; hw=obs_halfwidth)
    tobs = Float64.(obs[2, :])
    js = 1:size(obs, 2)
    lav = [log_offset(lava_at(d.series, tobs[j], hw), obs[3, j], obs[4, j]) for j in js]
    man = [log_offset(bulk_at(d.series, tobs[j]),     obs[5, j], obs[6, j]) for j in js]
    rms(v) = (u = filter(isfinite, v); isempty(u) ? NaN : sqrt(mean(abs2, u)))
    bias(v) = (u = filter(isfinite, v); isempty(u) ? NaN : mean(u))
    vol = log10(d.erupt_km3 / obs_volume)
    return (lava_rms = rms(lav), lava_bias = bias(lav), mantle_rms = rms(man),
            vol = vol, hit = count(x -> isfinite(x) && x == 0, vcat(lav, man)),
            n = count(isfinite, vcat(lav, man)),
            total = sqrt(mean(abs2, filter(isfinite, vcat(lav, man, vol)))))
end

if MODE == :heatmap
    axs, _ = sweep_axes(sweep, cs)
    isempty(axs) && error("the pf files in $sweep are identical - nothing to lay a heatmap on")
    ks = isempty(AXES) ? sort(collect(keys(axs))) : AXES
    for k in ks
        haskey(axs, k) || error("AXES asks for $k but every case has the same value for it")
    end
    length(ks) > 2 && @warn "sweep varies $(length(ks)) parameters; using the first two - set AXES to choose" ks
    xk = ks[1]; yk = length(ks) > 1 ? ks[2] : nothing
    xv = axs[xk]; yv = isnothing(yk) ? fill("", length(cs)) : axs[yk]
    xs, ys = ordered_unique(xv), ordered_unique(yv)
    # Two cases on the same cell mean a third parameter varies as well: warn rather than let one
    # silently overwrite the other, and name them so SELECT or AXES can separate them.
    owner = Dict{Tuple{Int,Int},String}()
    for (n, c) in enumerate(cs)
        ij = (findfirst(==(xv[n]), xs), findfirst(==(yv[n]), ys))
        haskey(owner, ij) && @warn "cases share a cell ($xk=$(xv[n]), $(something(yk,""))=$(yv[n])); only the last is shown" first=owner[ij] second=c
        owner[ij] = c
    end
    grid(f) = begin
        Z = fill(NaN, length(xs), length(ys))
        for (ij, c) in owner
            Z[ij...] = f(diag[c])
        end
        Z
    end
    fig = Figure(size=(1500, 950))
    Label(fig[0, 1:4], "$(basename(sweep)): $(length(cs)) runs, diagnostics after t = $(tstart) Gyr",
          fontsize=15, font=:bold)
    yname = isnothing(yk) ? "" : yk
    sweep_heatmap!(fig, (1, 1), xs, ys, xk, yname, grid(d -> d.erupt_km3),
                   "Erupted volume [km³]  (mare ≈ 1e7)", v -> @sprintf("%.2e", v), logscale=true)
    sweep_heatmap!(fig, (1, 3), xs, ys, xk, yname, grid(d -> d.lava),
                   "Lava H₂O, volume-weighted [μg/g]", v -> @sprintf("%.0f", v))
    sweep_heatmap!(fig, (2, 1), xs, ys, xk, yname, grid(d -> d.bulk),
                   "Bulk mantle H₂O [μg/g]", v -> @sprintf("%.2f", v))
    sweep_heatmap!(fig, (2, 3), xs, ys, xk, yname, grid(d -> d.outer),
                   "Outer-layer H₂O [μg/g]", v -> @sprintf("%.2f", v))
    display(fig)
    (savein != "") && save(savein * ".png", fig)

    if STATS
        obs = observations()
        sc = Dict(c => scores(diag[c], obs) for c in cs)
        gridS(f) = begin
            Z = fill(NaN, length(xs), length(ys))
            for (ij, c) in owner; Z[ij...] = f(sc[c]); end
            Z
        end
        for c in cs
            s = sc[c]
            @printf("%-14s  lava RMS=%5.2f (bias %+5.2f)  mantle RMS=%5.2f  volume=%+5.2f  inside %d/%d  total=%5.2f  [log10]\n",
                    c, s.lava_rms, s.lava_bias, s.mantle_rms, s.vol, s.hit, s.n, s.total)
        end
        # Distance from the data is one-sided (0 is best): sequential map. Bias and the volume
        # score have a sign, so they keep the diverging map centred on zero.
        sym(Z) = (m = maximum(abs, filter(isfinite, Z)); (-m, m))
        f2 = Figure(size=(1500, 950))
        Label(f2[0, 1:4], "$(basename(sweep)): distance from the Table2 samples, log10 units " *
              "(0 = inside the measured range; ±0.3 = a factor of 2)", fontsize=15, font=:bold)
        Zb = gridS(s -> s.lava_bias); Zv = gridS(s -> s.vol)
        sweep_heatmap!(f2, (1, 1), xs, ys, xk, yname, gridS(s -> s.lava_rms),
                       "Lava H₂O: RMS distance from the sample ranges",
                       v -> @sprintf("%.2f", v), cmap=MISFIT_CMAP)
        sweep_heatmap!(f2, (1, 3), xs, ys, xk, yname, Zb,
                       "Lava H₂O: mean signed distance (+ = too wet)",
                       v -> @sprintf("%+.2f", v), crange=sym(Zb))
        sweep_heatmap!(f2, (2, 1), xs, ys, xk, yname, gridS(s -> s.mantle_rms),
                       "Mantle H₂O: RMS distance from the source estimates",
                       v -> @sprintf("%.2f", v), cmap=MISFIT_CMAP)
        sweep_heatmap!(f2, (2, 3), xs, ys, xk, yname, Zv,
                       "Erupted volume vs $(@sprintf("%.0e", obs_volume)) km³ (+ = too much)",
                       v -> @sprintf("%+.2f", v), crange=sym(Zv))
        display(f2)
        (savein != "") && save(savein * "_stats.png", f2)
    end

# ------------------------------------------------------------------------
# Water mode: hydration state and storage capacity of the selected runs
# ------------------------------------------------------------------------
# Panels: (a) water and capacity against depth at five times, (b) saturation in depth and time,
#  (c) what each shell holds, (d) the mantle budget against what left it, (e) the free fraction.
elseif MODE == :water
    main = cs[1]                      # the run panel (b) is drawn for
    ages(t) = 4.5 .- t                # model time -> age, so the x axis reads like the samples
    CC_ = [length(cs) == 1 ? :black : Makie.wong_colors()[mod1(n, 7)] for n in eachindex(cs)]
    TIMES  = [0.5, 1.0, 2.0, 3.0, 4.4]      # Gyr; nearest rprof frame is used
    SHELLS = [(0.0, 40.0, "0-40 km (crust)"), (40.0, 300.0, "40-300 km (lid)"),
              (300.0, 450.0, "300-450 km (melting)"), (450.0, 900.0, "450-900 km"),
              (900.0, Inf, "900 km - CMB")]
    TC_ = cgrad(:viridis, length(TIMES), categorical=true)
    SC_ = [:firebrick, :darkorange, :seagreen, :steelblue, :purple]
    LS_ = [:solid, :dash, :dot, :dashdot]
    SUP_ = Dict(-3 => "⁻³", -2 => "⁻²", -1 => "⁻¹", 0 => "⁰",
                1 => "¹", 2 => "²", 3 => "³", 4 => "⁴", 5 => "⁵",
                6 => "⁶", 7 => "⁷", 8 => "⁸")
    # log axes reject zeros and the dry shells are exactly zero, so blank them instead
    pos(v) = [(x > 0 && isfinite(x)) ? x : NaN for x in v]

    # per-run profile handles: column index, depth, frame times and the base of the wet layer
    RP = Dict(c => let db = sims[c]
                  iTc, iRc, _ = data_encoding(db.timeheader, db.rprofheader, db.platesheader)
                  rr = db.rprofdata[:, 1, iRc["r"]]
                  wl = let v = get(parse_pf(joinpath(sweep, c, "pf")), "water_layer_depth", "")
                           something(tryparse(Float64, replace(v, r"[dD]" => "e")), outer_default)
                       end
                  (iT = iTc, iR = iRc, zc = (maximum(rr) .- rr) ./ 1e3, tp = db.rproftime, wl = wl / 1e3)
              end for c in cs)

    fig = Figure(size=(1100, 950))
    Label(fig[0, 1:3], "Mantle hydration and storage capacity: " * join(cs, ", "),
          fontsize=15, font=:bold)

    # -------- (a) water against the capacity that can hold it --------
    axa = Axis(fig[1, 1], xscale=log10, yreversed=true, titlesize=13,
               title="(a) Water and storage capacity vs depth",
               xlabel="H₂O [μg/g]", ylabel="depth [km]")
    for (n, c) in enumerate(cs), (k, tt) in enumerate(TIMES)
        p = RP[c]; j = argmin(abs.(p.tp .- tt))
        # one run: colour says the time. several runs: colour says the run, transparency the time.
        col = length(cs) == 1 ? TC_[k] : (CC_[n], 0.22 + 0.78 * (k - 1) / (length(TIMES) - 1))
        lines!(axa, pos(sims[c].rprofdata[:, j, p.iR["Water"]] .* WT2UGG), p.zc,
               color=col, linewidth=2)
        lines!(axa, pos(sims[c].rprofdata[:, j, p.iR["Wsol"]] .* WT2UGG), p.zc,
               color=col, linewidth=1.5, linestyle=:dash)
    end
    for (n, c) in enumerate(cs)
        hlines!(axa, [RP[c].wl], color=(length(cs) == 1 ? :grey25 : CC_[n], 0.8),
                linestyle=:dot, linewidth=1)
    end
    text!(axa, 1.3e-1, RP[main].wl, text=" water layer base", fontsize=8,
          color=:grey25, align=(:left, :bottom))
    axa.subtitle = "solid = water, dashed = capacity, dotted = base of the initial wet layer" *
                   (length(cs) == 1 ? "" : "; time shown by transparency (pale = early)")
    axa.subtitlesize = 9
    ela = length(cs) == 1 ? [LineElement(color=TC_[k], linewidth=2) for k in eachindex(TIMES)] :
                            [LineElement(color=CC_[n], linewidth=2) for n in eachindex(cs)]
    laba = length(cs) == 1 ? [@sprintf("%.1f Gyr", t) for t in TIMES] : cs
    xlims!(axa, 1e-1, 1e6)   # the dry crust tail runs to 1e-25 and would flatten everything else
    Legend(fig[1, 1], ela, laba, tellwidth=false, tellheight=false, halign=0.02, valign=0.02,
           framevisible=false, labelsize=9, patchsize=(18, 8))

    # -------- (b) how close to saturation, everywhere and always --------
    dbm = sims[main]; pm = RP[main]
    jm = findall(>(tstart), pm.tp)          # the same dry start the other panels drop
    agm = ages(pm.tp[jm]); ot = sortperm(agm); oz = sortperm(pm.zc)
    # saturation is 0.1 % over most of the mantle and 100 % in thin layers: log colour, or the
    #  whole panel reads as empty
    Ssat = clamp.(dbm.rprofdata[:, jm, pm.iR["H2Osat"]], 1e-2, 100)
    Fmlt = dbm.rprofdata[:, jm, pm.iR["fmeltmean"]]
    axb = Axis(fig[2, 1:2], xreversed=true, yreversed=true, titlesize=13,
               title="(b) $main: saturation of the storage capacity",
               xlabel="age [Ga]", ylabel="depth [km]")
    hmb = heatmap!(axb, agm[ot], pm.zc[oz], permutedims(Ssat[oz, ot]),
                   colormap=:viridis, colorrange=(1e-2, 100), colorscale=log10)
    contour!(axb, agm[ot], pm.zc[oz], permutedims(Fmlt[oz, ot]), levels=[1e-4],
             color=:white, linewidth=1.2)
    hlines!(axb, [pm.wl], color=(:white, 0.6), linestyle=:dot, linewidth=1)
    # the shells of (c): a line on each boundary and a colour tab at the right edge, so the two
    #  panels read against each other
    zmax = maximum(pm.zc); xr = maximum(agm) - minimum(agm)
    hlines!(axb, [sh[1] for sh in SHELLS[2:end]], color=(:white, 0.9), linewidth=1)
    for (k, (z1, z2, _)) in enumerate(SHELLS)
        z2 = min(z2, zmax)
        poly!(axb, Point2f[(minimum(agm) + 0.012xr, z1), (minimum(agm), z1),
                           (minimum(agm), z2), (minimum(agm) + 0.012xr, z2)], color=SC_[k], strokewidth=0)
        text!(axb, minimum(agm) + 0.016xr, (z1 + z2) / 2, text=split(SHELLS[k][3], " (")[1],
              fontsize=8, color=:white, font=:bold, align=(:right, :center))
    end
    axb.subtitle = "white contour = fmelt 1e-4, the melting region; tabs at the right = the shells of (c)"
    axb.subtitlesize = 9
    Colorbar(fig[2, 3], hmb, label="bound H₂O / capacity [%]", width=12, labelsize=10,
             ticklabelsize=9, ticks=([1e-2, 1e-1, 1e0, 1e1, 1e2], ["0.01", "0.1", "1", "10", "100"]))

    # -------- (c) what each shell holds, against the bulk --------
    axc = Axis(fig[1, 2], yscale=log10, xreversed=true, titlesize=13,
               title="(c) Shell inventories vs time", xlabel="age [Ga]", ylabel="H₂O [μg/g]")
    for (n, c) in enumerate(cs)
        p = RP[c]; db = sims[c]
        Wd = db.rprofdata[:, :, p.iR["Water"]]; dMd = db.rprofdata[:, :, p.iR["dM"]]
        ls = length(cs) == 1 ? :solid : LS_[mod1(n, length(LS_))]
        jj = findall(>(tstart), p.tp)          # the dry start would cost ten decades of the y axis
        for (s, (z1, z2, _)) in enumerate(SHELLS)
            m = (p.zc .>= z1) .& (p.zc .< z2)
            any(m) || continue
            v = [wmean(Wd[m, j], dMd[m, j]) * WT2UGG for j in jj]
            lines!(axc, ages(p.tp[jj]), pos(v), color=SC_[s], linewidth=2, linestyle=ls)
        end
        kt = findall(>(tstart), db.timedata[:, p.iT["time"]])
        lines!(axc, ages(db.timedata[kt, p.iT["time"]]),
               pos(db.timedata[kt, p.iT["cH2O_mean"]] .* WT2UGG), color=:black, linewidth=2,
               linestyle=ls)
    end
    elc = vcat([LineElement(color=SC_[s], linewidth=2) for s in eachindex(SHELLS)],
               [LineElement(color=:black, linewidth=2)])
    labc = vcat([sh[3] for sh in SHELLS], ["bulk mantle"])
    if length(cs) > 1
        elc = vcat(elc, [LineElement(color=:grey40, linewidth=2, linestyle=LS_[mod1(n, length(LS_))])
                         for n in eachindex(cs)])
        labc = vcat(labc, cs)
    end
    # decades only: over a range of one or two the default ticks read 10^0.5
    axc.yticks = ([10.0^k for k in -3:8], ["10" * SUP_[k] for k in -3:8])
    Legend(fig[1, 2], elc, labc, tellwidth=false, tellheight=false, halign=0.99, valign=0.02,
           framevisible=false, labelsize=9, nbanks=2, patchsize=(18, 8))

    rowsize!(fig.layout, 2, Relative(0.42))

    display(fig)
    (savein != "") && save(savein * ".png", fig)
# ------------------------------------------------------------------------
# Fit mode: one offset plane per sample stage, every run of every sweep
# ------------------------------------------------------------------------
# Per stage: x = how far the lava erupted in that window sits from the stage's primary-melt
#  range, y = the same for the melting-zone source against the stage's mantle range; 0,0 = both in.
elseif MODE == :fit
    obs  = observations()
    ogrp = String.(obs[1, :])
    omlo, omhi = Float64.(obs[3, :]), Float64.(obs[4, :])
    oslo, oshi = Float64.(obs[5, :]), Float64.(obs[6, :])
    GC_ = Dict("Apollo" => :steelblue, "CE5" => :firebrick, "CE6" => :darkorange)
    CC_ = [length(cs) == 1 ? :black : Makie.wong_colors()[mod1(n, 7)] for n in eachindex(cs)]
    root = dirname(sweep)
    swdirs = vcat([joinpath(root, d) for d in readdir(root) if isdir(joinpath(root, d)) && !(d in TRADE_SKIP)],
                  filter(isdir, TRADE_EXTRA))
    allc = [(sw, c) for sw in swdirs for c in readdir(sw)
            if isfile(joinpath(sw, c, "+op", "Moon_rprof.dat")) && isfile(joinpath(sw, c, "+op", "Moon_time.dat"))]
    winsrc(r, w) = (k = (r.t .> w[1]) .& (r.t .<= w[2]) .& isfinite.(r.src); any(k) ? mean(r.src[k]) : NaN)
    winlava(s, w) = (k = (s.t .> w[1]) .& (s.t .<= w[2]) .& (s.dE .> 0) .& (s.lava .> 1e-6);
                     any(k) ? wmean(s.lava[k], s.dE[k]) : NaN)
    runs = [(c = c, sw = basename(sw), sel = sw == sweep && c in cs,
             r = cached(() -> read_rprof_src(joinpath(sw, c, "+op", "Moon_rprof.dat")),
                        ("src", joinpath(sw, c, "+op", "Moon_rprof.dat")), joinpath(sw, c, "+op", "Moon_rprof.dat")),
             s = cached(() -> read_timedat(joinpath(sw, c, "+op", "Moon_time.dat")),
                        ("time", joinpath(sw, c, "+op", "Moon_time.dat")), joinpath(sw, c, "+op", "Moon_time.dat"))) for (sw, c) in allc]
    rng(v, j, f) = (u = filter(isfinite, v[j]); isempty(u) ? NaN : f(u))

    fig = Figure(size=(1500, 540))
    axs = Axis[]; allxy = Float64[]
    for (i, (name, w)) in enumerate(WINDOWS)
        j = findall(==(name), ogrp)
        flo, fhi = FIT_RANGE == :central ? (median, median) : (minimum, maximum)
        llo, lhi = rng(omlo, j, flo), rng(omhi, j, fhi)   # lava range of the stage across its estimates
        slo, shi = rng(oslo, j, flo), rng(oshi, j, fhi)   # source range of the stage
        # distance to the CENTRE of the range (geometric mean), so runs inside the box still rank
        mid(m, lo, hi) = (isfinite(m) && m > 0) ? log10(m / sqrt(lo * hi)) : NaN
        hx, hy = 0.5 * log10(lhi / llo), 0.5 * log10(shi / slo)   # half-widths of the ranges in dex
        pts = [(c = q.c, sw = q.sw, sel = q.sel,
                x = mid(winlava(q.s, w), llo, lhi), y = mid(winsrc(q.r, w), slo, shi)) for q in runs]
        filter!(p -> isfinite(p.x) && isfinite(p.y), pts)   # no eruption in the window = nothing to compare
        append!(allxy, vcat([p.x for p in pts], [p.y for p in pts]))
        ord = sortperm([hypot(p.x, p.y) for p in pts])
        @printf("(%s) stage %d %-6s lava %.0f-%.0f, source %.0f-%.0f ug/g; closest runs (lava, source offset):\n",
                'a' + i - 1, i, name, llo, lhi, slo, shi)
        for p in pts[ord[1:min(6, end)]]
            @printf("   %-10s %-18s %+.2f %+.2f\n", p.sw, p.c, p.x, p.y)
        end
        ax = Axis(fig[1, i], titlesize=13, xgridvisible=false, ygridvisible=false, aspect=1,
                  title=@sprintf("(%s) stage %d: %s, %.1f-%.1f Ga", 'a' + i - 1, i, name, 4.5 - w[2], 4.5 - w[1]),
                  xlabel="log₁₀ (window lava / centre of the $name melt range)",
                  ylabel=i == 1 ? "log₁₀ (window source / centre of the mantle range)" : "")
        ax.subtitle = @sprintf("melt %.0f-%.0f, mantle %.0f-%.0f μg/g; box = the ranges, 0 = their centres", llo, lhi, slo, shi)
        ax.subtitlesize = 9
        vspan!(ax, -hx, hx, color=(GC_[name], 0.10)); hspan!(ax, -hy, hy, color=(GC_[name], 0.10))
        poly!(ax, Point2f[(-hx, -hy), (hx, -hy), (hx, hy), (-hx, hy)],
              color=(GC_[name], 0.25), strokecolor=GC_[name], strokewidth=1.5, linestyle=:dash)
        hlines!(ax, [0.0], color=:grey60, linewidth=0.5); vlines!(ax, [0.0], color=:grey60, linewidth=0.5)
        oth = [p for p in pts if !p.sel]
        scatter!(ax, [p.x for p in oth], [p.y for p in oth], markersize=7, color=(:grey45, 0.8), strokewidth=0)
        for (n, c) in enumerate(cs)
            p = findfirst(q -> q.sel && q.c == c, pts)
            isnothing(p) && continue
            scatter!(ax, [pts[p].x], [pts[p].y], marker=:diamond, markersize=13, color=CC_[n],
                     strokewidth=1.2, strokecolor=:black)
            text!(ax, pts[p].x, pts[p].y, text=" " * c, fontsize=9, color=CC_[n], font=:bold,
                  offset=(6, 12 * (isodd(n) ? 1 : -1)), align=(:left, :center))
        end
        text!(ax, hx + 0.02, hy, text="$name ranges", fontsize=9, font=:bold, color=GC_[name], align=(:left, :top))
        push!(axs, ax)
    end
    # one scale for the three planes, so a run's drift from stage to stage reads across the row
    lim = max(maximum(abs, allxy), 0.5) + 0.1
    foreach(ax -> (xlims!(ax, -lim, lim); ylims!(ax, -lim, lim)), axs)
    display(fig)
    (savein != "") && save(savein * ".png", fig)
# ------------------------------------------------------------------------
# Sector mode: sample-site sectors straight from Moon_sector.dat, no StagPP load needed
# ------------------------------------------------------------------------
# Panels: (a) lava H2O by sample sector against the sample boxes, (b) the same for the melt
#  source, (c) cumulative erupted volume near vs far sectors, (d) window means and last eruptions.
elseif MODE == :sector
    obs  = observations()
    ogrp = String.(obs[1, :])
    omlo, omhi = Float64.(obs[3, :]), Float64.(obs[4, :])
    oslo, oshi = Float64.(obs[5, :]), Float64.(obs[6, :])
    rngf(v, j, f) = (u = filter(isfinite, v[j]); isempty(u) ? NaN : f(u))
    flo, fhi = FIT_RANGE == :central ? (median, median) : (minimum, maximum)
    GC_ = Dict("Apollo" => :steelblue, "CE6" => :darkorange, "CE5" => :firebrick)
    LS_ = [:solid, :dash, :dot, :dashdot]
    ages(t) = 4.5 .- t
    posv(v) = [(x > 0 && isfinite(x)) ? x : NaN for x in v]
    secdat = Dict(c => cached(() -> read_sectordat(joinpath(sweep, c, "+op", "Moon_sector.dat")),
                              ("sector", joinpath(sweep, c, "+op", "Moon_sector.dat")),
                              joinpath(sweep, c, "+op", "Moon_sector.dat")) for c in cs)

    fig = Figure(size=(1500, 900))

    # -------- (a) lava H2O by sample sector, against the sample boxes --------
    axa = Axis(fig[1, 1], yscale=log10, xreversed=true, titlesize=13,
               title="(a) lava H₂O by sample sector", xlabel="age [Ga]", ylabel="magma H₂O [μg/g]")
    for (name, (a, b)) in WINDOWS
        haskey(SECTOR_OF, name) || continue
        j = findall(==(name), ogrp)
        llo, lhi = rngf(omlo, j, flo), rngf(omhi, j, fhi)
        poly!(axa, Point2f[(4.5 - b, llo), (4.5 - a, llo), (4.5 - a, lhi), (4.5 - b, lhi)],
              color=(GC_[name], 0.20), strokecolor=GC_[name], strokewidth=1.2, linestyle=:dash)
    end
    hasa = false
    for (n, c) in enumerate(cs), (name, k) in SECTOR_OF
        s = secdat[c]; ls = length(cs) == 1 ? :solid : LS_[mod1(n, length(LS_))]
        e = tstart:0.1:4.5; bt, bv = Float64[], Float64[]
        keep = (s.t .> tstart) .& (s.dE[k] .> 0) .& (s.lava[k] .> 1e-6)
        for i in 1:length(e)-1
            m = keep .& (s.t .> e[i]) .& (s.t .<= e[i+1])
            any(m) || continue
            push!(bt, (e[i] + e[i+1]) / 2); push!(bv, wmean(s.lava[k][m], s.dE[k][m]))
        end
        isempty(bt) && continue
        lines!(axa, ages(bt), bv, color=GC_[name], linestyle=ls, linewidth=2,
               label=length(cs) == 1 ? name : "$c: $name")
        global hasa = true    # script-level for loop: plain assignment would shadow, not update
    end
    xlims!(axa, 0.0, 4.0)
    hasa && axislegend(axa, position=:rt, framevisible=false, labelsize=8, nbanks=2, merge=true)

    # -------- (b) source H2O of the same sectors, no boxes --------
    axb = Axis(fig[1, 2], yscale=log10, xreversed=true, titlesize=13,
               title="(b) source H₂O of the erupted magma, by sample sector",
               xlabel="age [Ga]", ylabel="source H₂O [μg/g]")
    hasb = false
    for (n, c) in enumerate(cs), (name, k) in SECTOR_OF
        s = secdat[c]; ls = length(cs) == 1 ? :solid : LS_[mod1(n, length(LS_))]
        e = tstart:0.1:4.5; bt, bv = Float64[], Float64[]
        keep = (s.t .> tstart) .& (s.dE[k] .> 0) .& (s.src[k] .> 1e-6)
        for i in 1:length(e)-1
            m = keep .& (s.t .> e[i]) .& (s.t .<= e[i+1])
            any(m) || continue
            push!(bt, (e[i] + e[i+1]) / 2); push!(bv, wmean(s.src[k][m], s.dE[k][m]))
        end
        isempty(bt) && continue
        lines!(axb, ages(bt), bv, color=GC_[name], linestyle=ls, linewidth=2,
               label=length(cs) == 1 ? name : "$c: $name")
        global hasb = true
    end
    xlims!(axb, 0.0, 4.0)
    hasb && axislegend(axb, position=:rt, framevisible=false, labelsize=8, nbanks=2, merge=true)

    # -------- (c) cumulative erupted volume, nearside vs SPA --------
    axc = Axis(fig[2, 1], yscale=log10, xreversed=true, titlesize=13,
               title="(c) erupted volume, nearside vs SPA", xlabel="age [Ga]", ylabel="cumulative volume [km³]")
    hspan!(axc, 0.5e7, 2e7, color=(:seagreen, 0.12))
    text!(axc, 3.9, sqrt(0.5e7 * 2e7), text="nearside target (Head & Wilson 1992, ±2x)",
          fontsize=9, color=:seagreen, align=(:left, :center))
    hlines!(axc, [1e5], color=:firebrick, linestyle=:dash, linewidth=1.5)
    text!(axc, 3.9, 1e5, text="SPA ceiling (~1 % of nearside)", fontsize=9, color=:firebrick,
          align=(:left, :bottom))
    for (n, c) in enumerate(cs)
        s = secdat[c]; ls = length(cs) == 1 ? :solid : LS_[mod1(n, length(LS_))]
        i0 = something(findlast(<=(tstart), s.t), 1)
        near = sum(s.E[k] for k in NEAR_SECTORS) .- sum(s.E[k][i0] for k in NEAR_SECTORS)
        far  = sum(s.E[k] for k in FAR_SECTORS)  .- sum(s.E[k][i0] for k in FAR_SECTORS)
        lines!(axc, ages(s.t), posv(near .* 4.784e9), color=:black, linestyle=ls, linewidth=2,
               label=length(cs) == 1 ? "near" : "$c: near")
        lines!(axc, ages(s.t), posv(far .* 4.784e9), color=:darkorange, linestyle=ls, linewidth=2,
               label=length(cs) == 1 ? "far" : "$c: far")
    end
    xlims!(axc, 0.0, 4.0)
    axislegend(axc, position=:lt, framevisible=false, labelsize=8, nbanks=2, merge=true)

    # -------- (d) window means, far/near ratio, last eruption age, per case --------
    axd = Axis(fig[2, 2], titlesize=13, title="(d) sample sector summary")
    hidedecorations!(axd); hidespines!(axd)
    lines_ = String[]
    for c in cs
        s = secdat[c]
        push!(lines_, c)
        for (name, (a, b)) in WINDOWS
            haskey(SECTOR_OF, name) || continue
            k = SECTOR_OF[name]
            keep = (s.t .> a) .& (s.t .<= b) .& (s.dE[k] .> 0) .& (s.lava[k] .> 1e-6)
            lv = any(keep) ? wmean(s.lava[k][keep], s.dE[k][keep]) : NaN
            sv = any(keep) ? wmean(s.src[k][keep], s.dE[k][keep]) : NaN
            push!(lines_, @sprintf("  %-6s  lava %6.1f  src %6.1f  ug/g", name, lv, sv))
        end
        i0 = something(findlast(<=(tstart), s.t), 1)
        vnear = sum(s.E[k][end] for k in NEAR_SECTORS) - sum(s.E[k][i0] for k in NEAR_SECTORS)
        vfar  = sum(s.E[k][end] for k in FAR_SECTORS)  - sum(s.E[k][i0] for k in FAR_SECTORS)
        push!(lines_, @sprintf("  far/near volume ratio = %.2e", vfar / max(vnear, eps())))
        anear = (dEg = sum(s.dE[k] for k in NEAR_SECTORS); ik = findlast(>(0), dEg); isnothing(ik) ? NaN : 4.5 - s.t[ik])
        afar  = (dEg = sum(s.dE[k] for k in FAR_SECTORS);  ik = findlast(>(0), dEg); isnothing(ik) ? NaN : 4.5 - s.t[ik])
        push!(lines_, @sprintf("  last eruption: near %.2f Ga, far %.2f Ga", anear, afar))
    end
    text!(axd, 0.02, 0.98, text=join(lines_, "\n"), fontsize=11, align=(:left, :top), space=:relative)
    xlims!(axd, 0, 1); ylims!(axd, 0, 1)

    display(fig)
    (savein != "") && save(savein * "_sector.png", fig)
# ------------------------------------------------------------------------
# Normal mode: one figure per selected run, built to be read against the samples
# ------------------------------------------------------------------------
# Panels: (a) lava H2O against every sample, (b) when the mare was emplaced, (c) magma against
#  its melting-zone source, (d) the match target by target, (e) the CE5 source fit against the late source.
else
    obs  = observations()
    ogrp = String.(obs[1, :]);  oage = 4.5 .- Float64.(obs[2, :])
    omlo, omhi = Float64.(obs[3, :]), Float64.(obs[4, :])
    oslo, oshi = Float64.(obs[5, :]), Float64.(obs[6, :])
    SHAPES_ = [:circle, :rect, :utriangle, :dtriangle, :pentagon, :hexagon, :star5, :cross]
    GC_ = Dict("Apollo" => :steelblue, "CE5" => :firebrick, "CE6" => :darkorange, "Rest" => :grey40)
    GM_ = Dict("Apollo" => :circle, "CE5" => :dtriangle, "CE6" => :utriangle, "Rest" => :xcross)
    main = cs[1]                      # the run the scorecard title and the melting lines refer to
    ages(t) = 4.5 .- t                # model time -> age, so the x axis reads like the samples
    # one colour per run, shared by (a), (b) and (d) so a line and its scorecard marker match
    CC_ = [length(cs) == 1 ? :black : Makie.wong_colors()[mod1(n, 7)] for n in eachindex(cs)]

    fig = Figure(size=(1700, 1120))

    # -------- (a) lava H2O: the panel everything else supports --------
    axa = Axis(fig[1, 1:2], yscale=log10, xreversed=true, title="(a) Erupted magma H₂O against sample estimates",
               xlabel="age [Ga]", ylabel="magma H₂O [μg/g]", titlesize=13)
    for (name, (a, b)) in WINDOWS          # the intervals the samples are averaged over
        vspan!(axa, 4.5 - b, 4.5 - a, color=(:grey, 0.10))
        text!(axa, 4.5 - (a + b) / 2, 3.3e3, text=name, align=(:center, :top), fontsize=10, color=:grey30)
    end
    for (n, c) in enumerate(cs)
        s = diag[c].series
        k = (s.t .> tstart) .& (s.dE .> 0) .& (s.lava .> 1e-6)
        # every eruption as a point sized by what it delivered, so a wet outlier that erupted
        #  nothing cannot masquerade as a trend
        w = s.dE[k] ./ maximum(s.dE[k])
        n == 1 && scatter!(axa, ages(s.t[k]), s.lava[k], markersize=2 .+ 9 .* sqrt.(w),
                           color=(:black, 0.16), strokewidth=0)
        e = tstart:0.1:4.5                 # volume-weighted mean in 100 Myr bins
        bt, bv = Float64[], Float64[]
        for i in 1:length(e)-1
            m = k .& (s.t .> e[i]) .& (s.t .<= e[i+1])
            any(m) || continue
            push!(bt, (e[i] + e[i+1]) / 2); push!(bv, wmean(s.lava[m], s.dE[m]))
        end
        lines!(axa, ages(bt), bv, linewidth=3, label=c, color=CC_[n])
    end
    for g in ("Apollo", "CE6", "CE5", "Rest")     # samples: the range each estimate allows
        j = findall(i -> ogrp[i] == g && isfinite(oage[i]) && isfinite(omlo[i]) && isfinite(omhi[i]), eachindex(ogrp))
        isempty(j) && continue
        rangebars!(axa, oage[j], omlo[j], omhi[j], color=GC_[g], linewidth=2.5, whiskerwidth=8)
        scatter!(axa, oage[j], sqrt.(omlo[j] .* omhi[j]), color=GC_[g], marker=GM_[g],
                 markersize=11, strokewidth=0.8, strokecolor=:white, label=g)
    end
    # stop the age axis a little past the youngest eruption; empty space hides the detail
    tyoung = minimum(4.5 - maximum(diag[c].series.t[diag[c].series.dE .> 0]) for c in cs)
    xend = max(0.0, tyoung - 0.3)
    ylims!(axa, 1, 4e3); xlims!(axa, 4.5 - tstart, xend)
    axislegend(axa, position=:lb, framevisible=false, labelsize=10, nbanks=2, merge=true)

    # row 2 gets its own grid, so the three panels share the width independently of row 1
    g2 = fig[2, 1:3] = GridLayout()

    # -------- (b) when the mare was emplaced --------
    axb = Axis(g2[1, 1], yscale=log10, xreversed=true, title="(b) Eruption history",
               xlabel="age [Ga]", ylabel="erupted volume per 100 Myr [km³]", titlesize=13)
    vspan!(axb, 3.2, 3.8, color=(:seagreen, 0.12))   # Hiesinger et al. (2011) main mare phase
    text!(axb, 3.5, 40, text="main mare\nphase", align=(:center, :bottom), fontsize=9, color=:seagreen)
    # crater-count ages of the youngest mare units end near 1.2 Ga (Hiesinger et al. 2011)
    vspan!(axb, MARE_END, 3.2, color=(:seagreen, 0.05))
    vlines!(axb, [MARE_END], color=(:seagreen, 0.7), linestyle=:dash, linewidth=1)
    text!(axb, (3.2 + max(MARE_END, xend)) / 2, 40, text="waning activity, to ~$(MARE_END) Ga",
          align=(:center, :bottom), fontsize=9, color=:seagreen)
    axb2 = Axis(g2[1, 1], yscale=log10, xreversed=true, yaxisposition=:right,
                ylabel="cumulative [km³]", ygridvisible=false)
    hidexdecorations!(axb2); linkxaxes!(axb, axb2)
    # decade ticks: over a narrow log range Makie otherwise labels 10^6.5 and similar
    axb2.yticks = ([1e4, 1e5, 1e6, 1e7, 1e8], ["10⁴", "10⁵", "10⁶", "10⁷", "10⁸"])
    for (n, c) in enumerate(cs)
        s = diag[c].series
        k = (s.t .> tstart) .& (s.dE .> 0)
        e = tstart:0.1:4.5
        bt = [(e[i] + e[i+1]) / 2 for i in 1:length(e)-1]
        bv = [sum(s.dE[k .& (s.t .> e[i]) .& (s.t .<= e[i+1])]) * KM3_PER_KG for i in 1:length(e)-1]
        col = CC_[n]
        m = bv .> 0
        stairs!(axb, ages(bt[m]), bv[m], color=col, linewidth=2, step=:center)
        lines!(axb2, ages(s.t[k]), cumsum(s.dE[k]) .* KM3_PER_KG, color=col, linewidth=2, linestyle=:dash)
    end
    hspan!(axb2, obs_volume * (1 - VOL_BAND), obs_volume * (1 + VOL_BAND), color=(:firebrick, 0.13))
    hlines!(axb2, [obs_volume], color=:firebrick, linestyle=:dot, linewidth=2)
    text!(axb2, xend + 0.1, obs_volume, text="mare ≈ $(@sprintf("%.0e", obs_volume)) km³ ± $(round(Int, 100VOL_BAND)) %",
          fontsize=9, color=:firebrick, align=(:right, :bottom))
    ylims!(axb, 30, nothing)     # floor under the band labels
    xlims!(axb, 4.5 - tstart, xend); xlims!(axb2, 4.5 - tstart, xend)

    # -------- (c) the magma against the source it melted from --------
    axc = Axis(g2[1, 2], xscale=log10, yscale=log10, titlesize=13,
               title="(c) Magma H₂O against its melting-zone source",
               xlabel="source H₂O [μg/g]", ylabel="magma H₂O [μg/g]",
               xgridvisible=false, ygridvisible=false)
    for g in ("Apollo", "CE6", "CE5")
        j = findall(i -> ogrp[i] == g && all(isfinite, (omlo[i], omhi[i], oslo[i], oshi[i])), eachindex(ogrp))
        isempty(j) && continue
        # envelope of the mission: dashed and faint, because most of it holds no sample
        x1, x2 = minimum(oslo[j]), maximum(oshi[j]); y1, y2 = minimum(omlo[j]), maximum(omhi[j])
        poly!(axc, Point2f[(x1, y1), (x2, y1), (x2, y2), (x1, y2)], color=(GC_[g], 0.06),
              strokecolor=(GC_[g], 0.75), strokewidth=1.2, linestyle=:dash)
        # union of the sample boxes on the grid of their own edges: one flat colour, no stacking
        xe = sort(unique(vcat(oslo[j], oshi[j]))); ye = sort(unique(vcat(omlo[j], omhi[j])))
        M = [any(i -> oslo[i] <= xe[a] && xe[a+1] <= oshi[i] && omlo[i] <= ye[q] && ye[q+1] <= omhi[i], j)
             for a in 1:length(xe)-1, q in 1:length(ye)-1]
        heatmap!(axc, xe, ye, [m ? 1.0 : NaN for m in M], colormap=[(GC_[g], 0.40), (GC_[g], 0.40)],
                 colorrange=(0, 1), nan_color=:transparent)
        # solid outline = the cell edges between covered and uncovered
        seg = Point2f[]
        for a in axes(M, 1), q in axes(M, 2)
            M[a, q] || continue
            (a == 1          || !M[a-1, q]) && push!(seg, (xe[a], ye[q]),   (xe[a], ye[q+1]))
            (a == size(M, 1) || !M[a+1, q]) && push!(seg, (xe[a+1], ye[q]), (xe[a+1], ye[q+1]))
            (q == 1          || !M[a, q-1]) && push!(seg, (xe[a], ye[q]),   (xe[a+1], ye[q]))
            (q == size(M, 2) || !M[a, q+1]) && push!(seg, (xe[a], ye[q+1]), (xe[a+1], ye[q+1]))
        end
        linesegments!(axc, seg, color=GC_[g], linewidth=1.5)
        # Apollo is named at its right end: the tracks start over its left corner
        text!(axc, g == "Apollo" ? x2 : x1, y2, text=g, fontsize=10, font=:bold, color=GC_[g],
              align=(g == "Apollo" ? :right : :left, :bottom), offset=(0, 1))
        g == "Apollo" && text!(axc, x2, y1, text="no sample envelope", fontsize=8, color=(GC_[g], 0.8),
                               align=(:right, :bottom), offset=(-3, 2))
    end
    # batch melting at fixed F, with the run's own partition coefficient
    Dw = let v = get(parse_pf(joinpath(sweep, main, "pf")), "dpartition_water", "")
             something(tryparse(Float64, replace(v, r"[dD]" => "e")), 0.01)
         end
    xsD = [0.3, 300.0]
    for F in (0.001, 0.01, 0.05, 0.2)
        lines!(axc, xsD, xsD ./ (Dw + F), color=(:grey40, 0.7), linestyle=:dash, linewidth=1)
        # label where the line leaves the box, so the four labels never stack
        xl = min(75.0, 2.2e3 * (Dw + F))
        text!(axc, xl, xl / (Dw + F), text=@sprintf("F=%g", F), fontsize=9, color=:grey30,
              align=(:right, :top), offset=(-2, -2))
    end
    # F -> 0: nothing left of this line can come from a source of that composition
    lines!(axc, xsD, xsD ./ Dw, color=(:grey25, 0.8), linewidth=1)
    text!(axc, 1.6, 1.6 / Dw, text="F → 0 limit ", fontsize=9, color=:grey25, align=(:right, :bottom))
    ecs = tstart:0.1:4.6
    for (n, c) in enumerate(cs)
        dbc = sims[c]
        _, iRc, _ = data_encoding(dbc.timeheader, dbc.rprofheader, dbc.platesheader)
        rrc = dbc.rprofdata[:, 1, iRc["r"]]; zc = (maximum(rrc) .- rrc) ./ 1e3
        Wc = dbc.rprofdata[:, :, iRc["Water"]]
        fmc = dbc.rprofdata[:, :, iRc["fmeltmean"]]; dMc = dbc.rprofdata[:, :, iRc["dM"]]
        tp = dbc.rproftime
        # source of a frame: the melting shells above the extraction cap, weighted by how much melts
        src = [let s = (zc .<= 420) .& (fmc[:, j] .> 0)
                   any(s) ? wmean(Wc[s, j] .* WT2UGG, fmc[s, j] .* dMc[s, j]) : NaN
               end for j in eachindex(tp)]
        s = diag[c].series
        k = (s.t .> tstart) .& (s.dE .> 0) .& (s.lava .> 1e-6)
        bt, bx, by = Float64[], Float64[], Float64[]
        for i in 1:length(ecs)-1
            m = k .& (s.t .> ecs[i]) .& (s.t .<= ecs[i+1])
            ms = (tp .> ecs[i]) .& (tp .<= ecs[i+1]) .& isfinite.(src)
            (any(m) && any(ms)) || continue
            sv = mean(src[ms])
            sv > 0 || continue
            push!(bt, (ecs[i] + ecs[i+1]) / 2); push!(bx, sv); push!(by, wmean(s.lava[m], s.dE[m]))
        end
        isempty(bt) && continue
        lines!(axc, bx, by, color=CC_[n], linewidth=2.5, label=c)
        # the boxes are not age-matched: fill the bins inside a mission's window with its colour
        inw = [findfirst(p -> p[2][1] < b <= p[2][2], collect(WINDOWS)) for b in bt]
        # line colour and shape say which run; the fill says which mission's box the bin belongs in
        shp = SHAPES_[mod1(n, length(SHAPES_))]
        q = isnothing.(inw)
        scatter!(axc, bx[q], by[q], color=CC_[n], marker=shp, markersize=6, strokewidth=0)
        scatter!(axc, bx[.!q], by[.!q], color=[GC_[WINDOWS[w][1]] for w in inw[.!q]], marker=shp,
                 markersize=13, strokewidth=1.2, strokecolor=:black)
        dy = n == 1 ? 10 : -10        # the tracks start on top of each other, so split the labels
        text!(axc, bx[1], by[1], text=@sprintf("%.1f Ga ", 4.5 - bt[1]), fontsize=9, color=CC_[n],
              align=(:right, :center), offset=(-3, dy))
        text!(axc, bx[end], by[end], text=@sprintf(" %.1f Ga", 4.5 - bt[end]), fontsize=9, color=CC_[n],
              align=(:left, :center), offset=(3, dy))
    end
    xlims!(axc, 0.7, 80); ylims!(axc, 8, 3e3)
    axc.xticks = ([1e0, 1e1, 1e2], ["10⁰", "10¹", "10²"])
    axc.yticks = ([1e1, 1e2, 1e3], ["10¹", "10²", "10³"])
    axc.subtitle = "line = run; filled marker = bin inside that mission's age window; dashed = fixed F (D = $(Dw))"
    axc.subtitlesize = 9
    # the missions are labelled on their boxes, so the legend only keys run to shape
    elc = [[LineElement(color=CC_[n], linewidth=2.5),
            MarkerElement(marker=SHAPES_[mod1(n, length(SHAPES_))], color=:white, markersize=10,
                          strokewidth=1, strokecolor=:black)] for n in eachindex(cs)]
    Legend(g2[1, 2], elc, cs, tellwidth=false, tellheight=false, halign=0.98, valign=0.04,
           framevisible=false, labelsize=8, margin=(0, 0, 0, 0))

    # -------- (e) the CE5 source fit against where the source ends up in the tail --------
    wce = Dict(WINDOWS)["CE5"]; wtail = (2.7, 4.5)
    j5 = findall(i -> ogrp[i] == "CE5" && isfinite(oslo[i]) && isfinite(oshi[i]), eachindex(ogrp))
    lo5, hi5 = minimum(oslo[j5]), maximum(oshi[j5])
    # every run of every sweep beside this one, read straight from rprof: no StagPP load
    root = dirname(sweep)
    swdirs = vcat([joinpath(root, d) for d in readdir(root) if isdir(joinpath(root, d)) && !(d in TRADE_SKIP)],
                  filter(isdir, TRADE_EXTRA))
    allc = [(sw, c) for sw in swdirs for c in readdir(sw) if isfile(joinpath(sw, c, "+op", "Moon_rprof.dat"))]
    winsrc(r, w) = (k = (r.t .> w[1]) .& (r.t .<= w[2]) .& isfinite.(r.src); any(k) ? mean(r.src[k]) : NaN)
    tro = [let r = cached(() -> read_rprof_src(joinpath(sw, c, "+op", "Moon_rprof.dat")),
                          ("src", joinpath(sw, c, "+op", "Moon_rprof.dat")), joinpath(sw, c, "+op", "Moon_rprof.dat")),
               s5 = winsrc(r, wce)
               (c = c, sw = basename(sw), sel = sw == sweep && c in cs,
                x = isfinite(s5) ? log_offset(s5, lo5, hi5) : NaN, y = winsrc(r, wtail))
           end for (sw, c) in allc]
    filter!(p -> isfinite(p.x) && isfinite(p.y), tro)   # no melting in a window = nothing to compare
    # the runs that sit on CE5 with the driest late source, printed so the cloud can be read by name
    ord = sortperm([(abs(p.x) > 0.05 ? 1 : 0, p.y) for p in tro])
    println("(e) on the CE5 source range with the driest tail source (offset, tail source ug/g):")
    for p in tro[ord[1:min(8, end)]]
        @printf("   %-10s %-18s %+.2f %6.1f\n", p.sw, p.c, p.x, p.y)
    end
    xall = [p.x for p in tro]; yall = [p.y for p in tro]
    xlo = minimum(xall) - 0.1; xhi = max(maximum(xall), 0.1) + 0.2
    ylo = min(0.5 * minimum(yall), 0.5); yhi = 2 * maximum(yall)
    axe = Axis(g2[1, 3], yscale=log10, titlesize=13, xgridvisible=false, ygridvisible=false,
               title="(e) CE5 source fit against the late source, all sweeps",
               xlabel="log₁₀ offset of the CE5-window source from the CE5 range (0 = inside)",
               ylabel="source H₂O after 2.7 Gyr [μg/g]")
    vspan!(axe, -0.02, 0.02, color=(:seagreen, 0.15))
    # single digits: the source at CE5 age is 4-15 ug/g and nothing younger says it got wetter
    hspan!(axe, 1.0, 10.0, color=(:seagreen, 0.10))
    poly!(axe, Point2f[(-0.05, 1.0), (0.05, 1.0), (0.05, 10.0), (-0.05, 10.0)],
          color=(:seagreen, 0.30), strokecolor=:seagreen, strokewidth=1.5, linestyle=:dash)
    text!(axe, xhi - 0.02, sqrt(10.0), text="late source in single digits", fontsize=9, font=:bold,
          color=:seagreen, align=(:right, :center))
    text!(axe, 0.0, yhi / 1.1, text="CE5 source target", fontsize=10, font=:bold, color=:seagreen,
          rotation=pi / 2, align=(:right, :center))
    oth = [p for p in tro if !p.sel]
    scatter!(axe, [p.x for p in oth], [p.y for p in oth], markersize=7, color=(:grey45, 0.8), strokewidth=0)
    for (n, c) in enumerate(cs)
        p = findfirst(q -> q.sel && q.c == c, tro)
        isnothing(p) && continue
        scatter!(axe, [tro[p].x], [tro[p].y], marker=:diamond, markersize=13, color=CC_[n],
                 strokewidth=1.2, strokecolor=:black)
        text!(axe, tro[p].x, tro[p].y, text=" " * c, fontsize=9, color=CC_[n], font=:bold,
              offset=(6, 12 * (isodd(n) ? 1 : -1)), align=(:left, :center))
    end
    xlims!(axe, xlo, xhi); ylims!(axe, ylo, yhi)
    axe.yticks = ([0.3, 1, 3, 10, 30, 100], ["0.3", "1", "3", "10", "30", "100"])
    axe.subtitle = "source = melt-weighted water of the melting shells above 420 km, as in (c)"
    axe.subtitlesize = 9

    colsize!(g2, 1, Relative(0.28)); colsize!(g2, 2, Relative(0.36)); colsize!(g2, 3, Relative(0.36))

    # -------- (d) the match, target by target --------
    # volume gets its own strip: km3 on the same log axis as ug/g would squash every other row
    gd = fig[1, 3] = GridLayout()
    axd = Axis(gd[1, 1], xscale=log10, titlesize=13, yticksvisible=false,
               xgridvisible=false, ygridvisible=false,
               title=length(cs) == 1 ? "(d) $main: model vs measured range" :
                                       "(d) model vs measured range",
               xlabel="H₂O [μg/g]")
    axv = Axis(gd[2, 1], xscale=log10, yticksvisible=false,
               xgridvisible=false, ygridvisible=false, xlabel="erupted volume [km³]")
    # one row per target: the measured range, then one marker per run in its (a) colour
    rng(v, j, f) = (u = filter(isfinite, v[j]); isempty(u) ? NaN : f(u))
    tg = Tuple{String,Float64,Float64,Symbol,Vector{Float64}}[]
    for (name, (a, b)) in WINDOWS
        j = findall(i -> ogrp[i] == name && isfinite(oage[i]), eachindex(ogrp))
        isempty(j) && continue
        push!(tg, ("$name magma", rng(omlo, j, minimum), rng(omhi, j, maximum), GC_[name],
                   [lava_at(diag[c].series, (a + b) / 2, (b - a) / 2) for c in cs]))
        push!(tg, ("$name source", rng(oslo, j, minimum), rng(oshi, j, maximum), GC_[name],
                   [bulk_at(diag[c].series, b) for c in cs]))
    end
    filter!(r -> isfinite(r[2]) && isfinite(r[3]), tg)
    # runs are stacked within their row, so two that agree do not hide one another
    offs = length(cs) == 1 ? [0.0] : collect(range(-0.24, 0.24, length(cs)))
    # drawn the same way on both strips, so the volume row reads like the H2O rows
    function markrow!(ax, y, lo, hi, gcol, vals, fmt)
        lines!(ax, [lo, hi], [y, y], color=(gcol, 0.35), linewidth=13)
        for (n, val) in enumerate(vals)
            isfinite(val) || continue
            ok = lo <= val <= hi
            # drop a line to the bar's level: it lands on the bar when the run is inside the range
            # dense dash: the segment is only a fraction of a row, so a default dash barely shows
            lines!(ax, [val, val], [y + offs[n], y], color=(CC_[n], 0.9),
                   linestyle=Makie.Linestyle([0.0, 1.8, 3.6]), linewidth=1.3)
            # fill says which run; inside or outside is read off the bar the marker sits against
            scatter!(ax, [val], [y + offs[n]], marker=:diamond, markersize=13, color=CC_[n],
                     strokewidth=1.2, strokecolor=:black)
            # label away from the row centre, or it lands on the measured-range bar
            dy = offs[n] == 0 ? 0.20 : (offs[n] > 0 ? 0.13 : -0.13)
            length(cs) <= 3 && text!(ax, val, y + offs[n] + dy, fontsize=8, color=CC_[n],
                align=(:center, dy > 0 ? :bottom : :top),
                text=ok ? fmt(val) : fmt(val) * @sprintf(" (%+.2f)", log_offset(val, lo, hi)))
        end
    end
    rows = String[]
    for (y, (lab, lo, hi, gcol, vals)) in enumerate(tg)
        push!(rows, lab)
        markrow!(axd, y, lo, hi, gcol, vals, v -> @sprintf("%.0f", v))
    end
    axd.yticks = (1:length(rows), rows)
    ylims!(axd, 0.35, length(rows) + 0.75)

    # the volume strip: the bar is the factor-2 band around the mare estimate the score uses
    vv = [diag[c].erupt_km3 for c in cs]
    markrow!(axv, 1, obs_volume / 2, obs_volume * 2, :seagreen, vv, v -> @sprintf("%.1e", v))
    axv.yticks = ([1], ["mare volume"])
    axv.xticks = ([1e4, 1e5, 1e6, 1e7, 1e8, 1e9], ["10⁴", "10⁵", "10⁶", "10⁷", "10⁸", "10⁹"])
    ylims!(axv, 0.45, 1.6)
    xlims!(axv, min(obs_volume / 8, minimum(vv) / 2), max(obs_volume * 8, maximum(vv) * 2))
    rowsize!(gd, 2, Relative(0.16))
    rowgap!(gd, 16)
    # one summary line per run, truncated rather than silently dropping the rest
    nsub = min(4, length(cs))
    axd.subtitle = join(vcat([(sc = scores(diag[c], obs);
                               @sprintf("%s: %d/%d samples in, RMS %.2f dex", c, sc.hit, sc.n,
                                        sc.total)) for c in cs[1:nsub]],
                             length(cs) > nsub ? ["+ $(length(cs) - nsub) more, see (a)"] : String[]), "\n")
    axd.subtitlesize = 9

    colsize!(fig.layout, 3, Relative(0.27))
    display(fig)
    (savein != "") && save(savein * ".png", fig)
end

@printf("cache: %d reused, %d loaded\n", LW_CACHE_HITS, LW_CACHE_MISS)
