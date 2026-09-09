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
    tmat = zeros(5, n)
    tmat[1, :] .= 4.5.-csv[2:2:end-1, 2]
    tmat[2, :] .= csv[2:2:end-1, 4]
    tmat[3, :] .= csv[3:2:end-1, 4]
    tmat[4, :] .= csv[2:2:end-1, 5]
    tmat[5, :] .= csv[3:2:end-1, 5]
    tmat = vcat(permutedims(grps), tmat)

    return tmat
end

tmat = load_data()

# parameters
xlabsz, ylabsz = 17, 17
xtklabsz, ytklabsz = 14, 14
lw = 2.5
tstart = 0.6
cmap = :rainbow

load_sims = false
if load_sims
    # low water
    # sim00001 = load_sim("/Volumes/Sinergia2/MoonWater/is00001_os0001_T1625_tc025_nt1w50_mid/+op/moon", "is00001_os0001_T1625_tc025_nt1w50_mid", plates=false);
    # sim00005 = load_sim("/Volumes/Sinergia2/MoonWater/is00005_os0001_T1625_tc025_nt1w50_mid/+op/moon", "is00005_os0001_T1625_tc025_nt1w50_mid", plates=false);
    # sim0001 = load_sim("/Volumes/Sinergia2/MoonWater/is0001_os0001_T1625_tc025_nt1w50_mid/+op/moon", "is0001_os0001_T1625_tc025_nt1w50_mid", plates=false);
    # sim00015 = load_sim("/Volumes/Sinergia2/MoonWater/is00015_os0001_T1625_tc025_nt1w50_mid/+op/moon", "is00015_os0001_T1625_tc025_nt1w50_mid", plates=false);
    # sim0002 = load_sim("/Volumes/Sinergia2/MoonWater/is0002_os0001_T1625_tc025_nt1w50_mid/+op/moon", "is0002_os0001_T1625_tc025_nt1w50_mid", plates=false);
    # sim00025 = load_sim("/Volumes/Sinergia2/MoonWater/is00025_os0001_T1625_tc025_nt1w50_mid/+op/moon", "is00025_os0001_T1625_tc025_nt1w50_mid", plates=false);
    # simos00015 = load_sim("/Volumes/Sinergia2/MoonWater/is0001_os00015_T1625_tc025_nt1w50_mid/+op/moon", "is0001_os00015_T1625_tc025_nt1w50_mid", plates=false);
    # simos0002 = load_sim("/Volumes/Sinergia2/MoonWater/is0001_os0002_T1625_tc025_nt1w50_mid/+op/moon", "is0001_os0002_T1625_tc025_nt1w50_mid", plates=false);

    # high water inner-outer
    # is001os0001 = load_sim("/Volumes/Sinergia2/MoonWater/is001_os0001_T1625_tc025_nt1w50_mid/+op/moon", "is001_os0001_T1625_tc025_nt1w50_mid", plates=false);
    # is001os01 = load_sim("/Volumes/Sinergia2/MoonWater/is001_os01_T1625_tc025_nt1w50_mid/+op/moon", "is001_os01_T1625_tc025_nt1w50_mid", plates=false);
    # is001os001 = load_sim("/Volumes/Sinergia2/MoonWater/old/is001_os001_T1625_tc025_nt2w20_mid/+op/moon", "is001_os001_T1625_tc025_nt2w20_mid", plates=false);
    # is0001os001 = load_sim("/Volumes/Sinergia2/MoonWater/is0001_os001_T1625_tc025_nt1w50_mid/+op/moon", "is0001_os001_T1625_tc025_nt1w50_mid", plates=false);
    # is01os001 = load_sim("/Volumes/Sinergia2/MoonWater/old/is01_os001_T1625_tc025_nt2w20_mid/+op/moon", "is01_os001_T1625_tc025_nt1w50_mid", plates=false);

    # inner-water less range nt1w50
    # is001os001 = load_sim("/Volumes/Sinergia2/MoonWater/is001_os001_T1625_tc025_nt1w50_mid/+op/moon", "is001_os001_T1625_tc025_nt1w50_mid", plates=false);
    # is001os0015 = load_sim("/Volumes/Sinergia2/MoonWater/is001_os0015_T1625_tc025_nt1w50_mid/+op/moon", "is001_os0015_T1625_tc025_nt1w50_mid", plates=false);
    # is001os002 = load_sim("/Volumes/Sinergia2/MoonWater/is001_os002_T1625_tc025_nt1w50_mid/+op/moon", "is001_os002_T1625_tc025_nt1w50_mid", plates=false);
    # is001os0025 = load_sim("/Volumes/Sinergia2/MoonWater/is001_os0025_T1625_tc025_nt1w50_mid/+op/moon", "is001_os0025_T1625_tc025_nt1w50_mid", plates=false);
    # is001os003 = load_sim("/Volumes/Sinergia2/MoonWater/is001_os003_T1625_tc025_nt1w50_mid/+op/moon", "is001_os003_T1625_tc025_nt1w50_mid", plates=false);
    # is0015os0015 = load_sim("/Volumes/Sinergia2/MoonWater/is0015_os0015_T1625_tc025_nt1w50_mid/+op/moon", "is0015_os0015_T1625_tc025_nt1w50_mid", plates=false);
    # is0015os002 = load_sim("/Volumes/Sinergia2/MoonWater/is0015_os002_T1625_tc025_nt1w50_mid/+op/moon", "is0015_os002_T1625_tc025_nt1w50_mid", plates=false);
    # is0015os0025 = load_sim("/Volumes/Sinergia2/MoonWater/is0015_os0025_T1625_tc025_nt1w50_mid/+op/moon", "is0015_os0025_T1625_tc025_nt1w50_mid", plates=false);
    # is0015os003 = load_sim("/Volumes/Sinergia2/MoonWater/is0015_os003_T1625_tc025_nt1w50_mid/+op/moon", "is0015_os003_T1625_tc025_nt1w50_mid", plates=false);

    # similar with erupt threshold
    # is0015os003_ef002_eth000 = load_sim("/Volumes/Sinergia2/MoonWater/is0015_os003_T1625_tc025_nt1w50_mid_ef002_eth000/+op/moon", "is0015_os003_T1625_tc025_nt1w50_mid_ef002_eth000", plates=false);
    # is0015os003_ef002_eth0005 = load_sim("/Volumes/Sinergia2/MoonWater/is0015_os003_T1625_tc025_nt1w50_mid_ef002_eth0005/+op/moon", "is0015_os003_T1625_tc025_nt1w50_mid_ef002_eth0005", plates=false);
    # is0015os003_ef005_eth000 = load_sim("/Volumes/Sinergia2/MoonWater/is0015_os003_T1625_tc025_nt1w50_mid_ef005_eth000/+op/moon", "is0015_os003_T1625_tc025_nt1w50_mid_ef005_eth000", plates=false);
    # is0015os003_ef005_eth0005 = load_sim("/Volumes/Sinergia2/MoonWater/is0015_os003_T1625_tc025_nt1w50_mid_ef005_eth0005/+op/moon", "is0015_os003_T1625_tc025_nt1w50_mid_ef005_eth0005", plates=false);

    # od - Default is0015_T1625_os0020_od300_tc025_nt1w50_mid_ef005
    od150 = load_sim("/Volumes/Sinergia2/MoonWater/old/T1625/is0015_os0020_od150_tc025_nt1w50_mid_ef005/+op/moon", "od150", plates=false);
    od200 = load_sim("/Volumes/Sinergia2/MoonWater/old/T1625/is0015_os0020_od200_tc025_nt1w50_mid_ef005/+op/moon", "od200", plates=false);
    # def = load_sim("/Volumes/Sinergia2/MoonWater/old/T1625/is0015_os0020_od300_tc025_nt1w50_mid_ef005/+op/moon", "od300", plates=false);
    od150_T1610 = load_sim("/Volumes/Sinergia2/MoonWater/old/T1610/is0015_os0020_od150_tc025_nt1w50_mid_ef005/+op/moon", "od150_T1610", plates=false);
    od200_T1610 = load_sim("/Volumes/Sinergia2/MoonWater/old/T1610/is0015_os0020_od200_tc025_nt1w50_mid_ef005/+op/moon", "od200_T1610", plates=false);
    od300_T1610 = load_sim("/Volumes/Sinergia2/MoonWater/old/T1610/is0015_os0020_od300_tc025_nt1w50_mid_ef005/+op/moon", "od300_T1610", plates=false);
    
    # tc
    tc010 = load_sim("/Volumes/Sinergia2/MoonWater/old/T1625/is0015_os0020_od300_tc010_nt1w50_mid_ef005/+op/moon", "tc010", plates=false);
    tc015 = load_sim("/Volumes/Sinergia2/MoonWater/old/T1625/is0015_os0020_od300_tc015_nt1w50_mid_ef005/+op/moon", "tc015", plates=false);
    tc010_T1610 = load_sim("/Volumes/Sinergia2/MoonWater/old/T1610/is0015_os0020_od300_tc010_nt1w50_mid_ef005/+op/moon", "tc010_T1610", plates=false);
    tc015_T1610 = load_sim("/Volumes/Sinergia2/MoonWater/old/T1610/is0015_os0020_od300_tc015_nt1w50_mid_ef005/+op/moon", "tc015_T1610", plates=false);
    # def
    # os
    os0015 = load_sim("/Volumes/Sinergia2/MoonWater/old/T1625/is0015_os0015_od300_tc025_nt1w50_mid_ef005/+op/moon", "os0015", plates=false);
    os0020 = load_sim("/Volumes/Sinergia2/MoonWater/old/T1625/is0015_os0020_od300_tc025_nt1w50_mid_ef005/+op/moon", "os0020", plates=false);
    os0020_T1610 = load_sim("/Volumes/Sinergia2/MoonWater/old/T1610/is0015_os0020_od300_tc025_nt1w50_mid_ef005/+op/moon", "os0020_T1610", plates=false);
    os0020_T1610 = load_sim("/Volumes/Sinergia2/MoonWater/old/T1610/is0015_os0020_od300_tc025_nt1w50_mid_ef005/+op/moon", "os0020_T1610", plates=false);
    # ef
    # def
    ef010 = load_sim("/Volumes/Sinergia2/MoonWater/old/T1625/is0015_os0020_od300_tc025_nt1w50_mid_ef010/+op/moon", "ef010", plates=false);
    ef020 = load_sim("/Volumes/Sinergia2/MoonWater/old/T1625/is0015_os0020_od300_tc025_nt1w50_mid_ef020/+op/moon", "ef020", plates=false);
    ef010_T1610 = load_sim("/Volumes/Sinergia2/MoonWater/old/T1610/is0015_os0020_od300_tc025_nt1w50_mid_ef010/+op/moon", "ef010_T1610", plates=false);
    ef020_T1610 = load_sim("/Volumes/Sinergia2/MoonWater/old/T1610/is0015_os0020_od300_tc025_nt1w50_mid_ef020/+op/moon", "ef020_T1610", plates=false);
    # Combo
    os0015_od200_tc015_ef010 = load_sim("/Volumes/Sinergia2/MoonWater/old/T1610/is0015_os0015_od200_tc015_nt1w50_mid_ef010/+op/moon", "os0015_od200_tc015_ef010", plates=false);

end

sim  = [od150, tc010, os0015, ef010] # 
simR = od150

# Thickness [m] of the initially water-enriched outer layer, read per simulation from the
# parameters StagYY wrote out (water_layer_depth). The layer occupies the top od metres of the
# 1340 km lunar mantle, so it starts at r = rkm - od measured from the CMB.
function outer_layer_depth(db::DataBlock)
    pf = joinpath(db.metadata.outdir, db.metadata.sroot * "_parameters.dat")
    m = match(r"WATER_LAYER_DEPTH\s*=\s*([0-9.eE+-]+)", read(pf, String))
    return parse(Float64, m.captures[1])
end
# figure selection
fig_select = 1

# plots
fig = Figure(size=(1400, 900))
tstart = 0.6
if fig_select==1
    time_vs_field(sim, "cH2O_mean", yrange=(0.1, 1e2), color=:black, logy=true, fig=fig, fpos=(1,1), tstart=tstart, xlabelsize=xlabsz, ylabelsize=ylabsz, xticklabelsize=xtklabsz, yticklabelsize=ytklabsz, yconv=(1e4, "μg/g"), linewidth=lw, error_bar_leg=true, error_bar_legpos=:rb,
                   error_bars=vcat(tmat[1:2,:], tmat[5:6,:]), leg=false, scatter=false, cmap=cmap)
    time_vs_field(sim, "Water", rprof_top=outer_layer_depth.(sim), ylabel=L"Outer\;layer\;c^{H_2O}", yrange=(0.1, 1e2), color=:black, logy=true, fig=fig, fpos=(2,1), tstart=tstart, xlabelsize=xlabsz, ylabelsize=ylabsz, xticklabelsize=xtklabsz, yticklabelsize=ytklabsz, yconv=(1e4, "μg/g"), linewidth=lw,
                    error_bars=vcat(tmat[1:2,:], tmat[5:6,:]), leg=false, scatter=false, cmap=cmap)
    time_vs_field(sim, "Wsol", yrange=(-1, -1), color=:black, logy=false, fig=fig, fpos=(2,2), tstart=tstart, xlabelsize=xlabsz, ylabelsize=ylabsz, xticklabelsize=xtklabsz, yticklabelsize=ytklabsz, yconv=(1e4, "μg/g"), linewidth=lw, cmap=cmap,
                    leg=true, scatter=false)
    time_vs_field(sim, "erupta", yrange=(-1, -1), color=:black, zero_at_tstart=true, logy=true, line=true, fig=fig, fpos=(1,2), tstart=tstart, xlabelsize=xlabsz, ylabelsize=ylabsz, xticklabelsize=xtklabsz, yticklabelsize=ytklabsz, yconv=(3.3e-12, "km³"), linewidth=lw, cmap=cmap,
                    leg=false, scatter=false)
    time_vs_field(sim, "∂erupta", yrange=(-1, -1), mavg_excl_floor=true, mavg_intensity=0.05, mov_avg=true, color=:black, logy=true, line=false, fig=fig, fpos=(1,3), tstart=tstart, xlabelsize=xlabsz, ylabelsize=ylabsz, xticklabelsize=xtklabsz, yticklabelsize=ytklabsz, yconv=(3.3e-12, "km³"), linewidth=lw, cmap=cmap, leg=false, mavg_contrast=0.1)
    time_vs_field(sim, "cH2Oerupt", mov_avg=true, yrange=(1., 1e4), color=:firebrick, logy=true, line=false, fig=fig, fpos=(2,3), tstart=tstart, xlabelsize=xlabsz, ylabelsize=ylabsz, xticklabelsize=xtklabsz, yticklabelsize=ytklabsz, linewidth=lw, cmap=cmap,
                    savein="test", yconv=(1e4, "μg/g"), mavg_excl_floor=true, mavg_intensity=0.05,
                    error_bars=tmat[1:4,:], leg=false, mavg_contrast=0.15)
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

