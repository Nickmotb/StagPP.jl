using StagPP, MAGEMin_C, GLMakie, DelimitedFiles

function load_data()
    # Load data
    path = "/Users/nickmotb/Desktop/LunarWater/Table2.csv"
    csv = readdlm(path, ','); replace!(csv, "" => NaN)
    n = Int(0.5(size(csv, 1)-2))

    replace!(csv, "AP11"=>"Apollo")
    replace!(csv, "AP12"=>"Apollo")
    replace!(csv, "AP15"=>"Apollo")
    replace!(csv, "AP17"=>"Apollo")
    replace!(csv, "Ap7"=>"Apollo")

    idxCE6 = findall(csv[2:2:end-1, 1].=="CE6")
    idxCE5 = findall(csv[2:2:end-1, 1].=="CE5")
    idxAp = findall(csv[2:2:end-1, 1].=="Apollo")

    # Store data as [i, :] -> i = [Study Time(Gyr) minH2O_pm maxH2O_pm minH2O_mtl maxH2O_mtl]
    grps = fill("Rest", n)
    grps[idxCE6]    .= "CE6"
    grps[idxAp]     .= "Apollo Samples"
    grps[idxCE5]    .= "CE5"
    mat = zeros(5, n)
    mat[1, :] .= 4.5.-csv[2:2:end-1, 2]
    mat[2, :] .= csv[2:2:end-1, 4]
    mat[3, :] .= csv[3:2:end-1, 4]
    mat[4, :] .= csv[2:2:end-1, 5]
    mat[5, :] .= csv[3:2:end-1, 5]
    mat = vcat(permutedims(grps), mat)

    return mat
end

mat = load_data()

# parameters
xlabsz, ylabsz = 17, 17
xtklabsz, ytklabsz = 14, 14
lw = 2.5
tstart = 0.6
cmap = :rainbow

load_sims = false
if load_sims
    # sim0003 = load_sim("/Volumes/Sinergia2/MoonWater/s0003_T1625_tc025_nt2w20/+op/moon", "s0003_T1625_tc025_nt2w20", plates=false);
    # sim0001 = load_sim("/Volumes/Sinergia2/MoonWater/s0001_T1625_tc025_nt1w50/+op/moon", "s0001_T1625_tc025_nt1w50", plates=false);
    # sim0001mid = load_sim("/Volumes/Sinergia2/MoonWater/s0001_T1625_tc025_nt1w50_mid/+op/moon", "s0001_T1625_tc025_nt1w50_mid", plates=false);
    # sim0005 = load_sim("/Volumes/Sinergia2/MoonWater/s0005_T1625_tc025_nt1w50/+op/moon", "s0005_T1625_tc025_nt1w50", plates=false);
    # sim0005mid = load_sim("/Volumes/Sinergia2/MoonWater/s0005_T1625_tc025_nt1w50_mid/+op/moon", "s0005_T1625_tc025_nt1w50_mid", plates=false);
    # sim0007 = load_sim("/Volumes/Sinergia2/MoonWater/s0007_T1625_tc025_nt1w50/+op/moon", "s0007_T1625_tc025_nt1w50", plates=false);
    # sim0007mid = load_sim("/Volumes/Sinergia2/MoonWater/s0007_T1625_tc025_nt1w50_mid/+op/moon", "s0007_T1625_tc025_nt1w50_mid", plates=false);

    sim00001 = load_sim("/Volumes/Sinergia2/MoonWater/is00001_os0001_T1625_tc025_nt1w50_mid/+op/moon", "is00001_os0001_T1625_tc025_nt1w50_mid", plates=false);
    # sim00005 = load_sim("/Volumes/Sinergia2/MoonWater/is00005_os0001_T1625_tc025_nt1w50_mid/+op/moon", "is00005_os0001_T1625_tc025_nt1w50_mid", plates=false);
    sim0002 = load_sim("/Volumes/Sinergia2/MoonWater/is0002_os0001_T1625_tc025_nt1w50_mid/+op/moon", "is0002_os0001_T1625_tc025_nt1w50_mid", plates=false);
    sim0004 = load_sim("/Volumes/Sinergia2/MoonWater/is0004_os0001_T1625_tc025_nt1w50_mid/+op/moon", "is0004_os0001_T1625_tc025_nt1w50_mid", plates=false);
    sim0006 = load_sim("/Volumes/Sinergia2/MoonWater/is0006_os0001_T1625_tc025_nt1w50_mid/+op/moon", "is0006_os0001_T1625_tc025_nt1w50_mid", plates=false);

end

sim  = [sim00001, sim0002, sim0004, sim0006]
simR = sim00001
# figure selection
fig_select = 1

# plots
fig = Figure(size=(1400, 900))
tstart = 0.6
if fig_select==1
    time_vs_field(sim, "cH2O_mean", yrange=(0.1, 1e2), color=:black, logy=true, fig=fig, fpos=(1,1), tstart=tstart, xlabelsize=xlabsz, ylabelsize=ylabsz, xticklabelsize=xtklabsz, yticklabelsize=ytklabsz, yconv=(1e4, "μg/g"), linewidth=lw, error_bar_leg=true, error_bar_legpos=:rb,
                   error_bars=vcat(mat[1:2,:], mat[5:6,:]), leg=false, scatter=false, cmap=cmap)
    time_vs_field(sim, "Wsol", yrange=(-1, -1), color=:black, logy=false, fig=fig, fpos=(2,1), tstart=tstart, xlabelsize=xlabsz, ylabelsize=ylabsz, xticklabelsize=xtklabsz, yticklabelsize=ytklabsz, yconv=(1e4, "μg/g"), linewidth=lw, cmap=cmap,
                    leg=true, scatter=false)
    time_vs_field(sim, "erupta", yrange=(-1, -1), color=:black, zero_at_tstart=true, logy=true, line=true, fig=fig, fpos=(1,2), tstart=tstart, xlabelsize=xlabsz, ylabelsize=ylabsz, xticklabelsize=xtklabsz, yticklabelsize=ytklabsz, yconv=(3.3e-12, "km³"), linewidth=lw, cmap=cmap,
                    leg=false, scatter=false)
    time_vs_field(sim, "∂erupta", yrange=(-1, -1), mavg_excl_floor=true, mavg_intensity=0.05, mov_avg=true, color=:black, logy=true, line=false, fig=fig, fpos=(1,3), tstart=tstart, xlabelsize=xlabsz, ylabelsize=ylabsz, xticklabelsize=xtklabsz, yticklabelsize=ytklabsz, yconv=(3.3e-12, "km³"), linewidth=lw, cmap=cmap, leg=false, mavg_contrast=0.1)
    time_vs_field(sim, "Tmean", color=:black, logy=false, line=true, fig=fig, fpos=(2,2), tstart=tstart, xlabelsize=xlabsz, ylabelsize=ylabsz, xticklabelsize=xtklabsz, yticklabelsize=ytklabsz, linewidth=lw, cmap=cmap,
                    leg=false, scatter=false)
    time_vs_field(sim, "cH2Oerupt", mov_avg=true, yrange=(1., 1e4), color=:firebrick, logy=true, line=false, fig=fig, fpos=(2,3), tstart=tstart, xlabelsize=xlabsz, ylabelsize=ylabsz, xticklabelsize=xtklabsz, yticklabelsize=ytklabsz, linewidth=lw, cmap=cmap,
                    savein="", yconv=(1e4, "μg/g"), mavg_excl_floor=true, mavg_intensity=0.05,
                    error_bars=mat[1:4,:], leg=false, mavg_contrast=0.15)
elseif fig_select==2
    rprof_vs_field(simR, "Wsol", tstart=tstart, cmap=:Greys, logy=true, pres=false, colorrange=(-2.5, -1.), fig=fig, fpos=(1,1), interpolate=true, leg=false)
    rprof_vs_field(simR, "boundH2O", tstart=tstart,  cmap=:Blues, cmap_reverse=false, colorrange=(-4, -2), logy=true, pres=false, fig=fig, fpos=(1,2), interpolate=true, leg=false)
    rprof_vs_field(simR, "fmeltmean", tstart=tstart,  cmap=:thermal, cmap_reverse=false, colorrange=(-6, -1.), logy=true, pres=false, fig=fig, fpos=(3,1), interpolate=true, leg=false)
    rprof_vs_field(simR, "freeH2O", tstart=tstart,  cmap=:Blues, cmap_reverse=false, colorrange=(-6, -2), logy=true, pres=false, fig=fig, fpos=(1,3), interpolate=true, leg=true)
    rprof_vs_field(simR, "bsmean", tstart=tstart,  cmap=:RdGy, cmap_reverse=true, colorrange=(0, 1), logy=false, pres=false, fig=fig, fpos=(3,3), interpolate=true, leg=false)
    rprof_vs_field(simR, "H2Osat", tstart=tstart,  cmap=:vik100, cmap_reverse=false, colorrange=(-3, 0.), logy=true, pres=false, fig=fig, fpos=(3,2), interpolate=true, yconv=(1e-2, ""), leg=false, savein="rprof_s0001mid")
end
display(fig)

# time_vs_field(sim, "H2Odarcy", logy=true)

