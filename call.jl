using StagPP, GLMakie, MAGEMin_C

load_sims = false
# Sims
    if load_sims
        # Full Runs
        y20 = load_sim("/Volumes/Sinergia/EarthWater/y20e19/+op/EW", "γ20");
        y25 = load_sim("/Volumes/Sinergia/EarthWater/y25e19/+op/EW", "γ25");
        y30 = load_sim("/Volumes/Sinergia/EarthWater/y30e19/+op/EW", "γ30");
        y35 = load_sim("/Volumes/Sinergia/EarthWater/y35e19/+op/EW", "γ35");
        y40 = load_sim("/Volumes/Sinergia/EarthWater/y40e19/+op/EW", "γ40");
        y45 = load_sim("/Volumes/Sinergia/EarthWater/y45e19/+op/EW", "γ45");
        y50 = load_sim("/Volumes/Sinergia/EarthWater/y50e19/+op/EW", "γ50");
        y55 = load_sim("/Volumes/Sinergia/EarthWater/y55e19/+op/EW", "γ55");

        # No Water
        y20nw = load_sim("/Volumes/Sinergia/EarthWater/y20e19nw/+op/EW", "γ20nw");
        y25nw = load_sim("/Volumes/Sinergia/EarthWater/y25e19nw/+op/EW", "γ25nw");
        y30nw = load_sim("/Volumes/Sinergia/EarthWater/y30e19nw/+op/EW", "γ30nw");
        y35nw = load_sim("/Volumes/Sinergia/EarthWater/y35e19nw/+op/EW", "γ35nw");
        y40nw = load_sim("/Volumes/Sinergia/EarthWater/y40e19nw/+op/EW", "γ40nw");
        y45nw = load_sim("/Volumes/Sinergia/EarthWater/y45e19nw/+op/EW", "γ45nw");
        y50nw = load_sim("/Volumes/Sinergia/EarthWater/y50e19nw/+op/EW", "γ50nw");
        y55nw = load_sim("/Volumes/Sinergia/EarthWater/y55e19nw/+op/EW", "γ55nw");

        # Ingas depth
        id4 = load_sim("/Volumes/Sinergia/EarthWater/id4/+op/EW", "id4");
        id6 = load_sim("/Volumes/Sinergia/EarthWater/id6/+op/EW", "id6");
        id8 = load_sim("/Volumes/Sinergia/EarthWater/id8/+op/EW", "id8");
        
        # No RK
        noR = load_sim("/Volumes/Sinergia/EarthWater/noR/+op/EW", "noR");
        noK = load_sim("/Volumes/Sinergia/EarthWater/noK/+op/EW", "noK");
        noKR = load_sim("/Volumes/Sinergia/EarthWater/noKR/+op/EW", "noKR");

        # Total Water Budget
        tot1 = load_sim("/Volumes/Sinergia/EarthWater/tot1/+op/EW", "1 OM");
        tot2 = load_sim("/Volumes/Sinergia/EarthWater/tot2/+op/EW", "2 OM");
        tot3 = load_sim("/Volumes/Sinergia/EarthWater/tot3/+op/EW", "3 OM");
        tot4 = load_sim("/Volumes/Sinergia/EarthWater/tot4/+op/EW", "4 OM");
        tot5 = load_sim("/Volumes/Sinergia/EarthWater/tot5/+op/EW", "5 OM");

        # # sets
        y = [y20, y25, y30, y35, y40, y45, y50, y55];
        ynw = [y20nw, y25nw, y30nw, y35nw, y40nw, y45nw, y50nw, y55nw];
        yt = vcat(y, ynw)
        id = [id4, id6, id8, y20];
        tot = [tot1, tot2, tot3, tot4, tot5, y20];
        no = [y20, noR, noK, noKR]
    end

# fig = Figure(size=(1600, 800))
# lbsz = 20
# time_vs_field(y, "Wsol_um", line=true, fig=fig, fpos=(1,1), framevisible=false, linewidth=1.7, markersize=7, ylabelsize=lbsz, xlabelsize=lbsz)
# time_vs_field(y, "Water_um", line=true, fig=fig, fpos=(1,2), leg=false, linewidth=1.7, markersize=7, ylabelsize=lbsz, xlabelsize=lbsz)
# time_vs_field(y, "Wsol_tz", line=true, fig=fig, fpos=(2,1), leg=false, linewidth=1.7, markersize=7, ylabelsize=lbsz, xlabelsize=lbsz)
# time_vs_field(y, "Water_tz", line=true, fig=fig, fpos=(2,2), leg=false, linewidth=1.7, markersize=7, ylabelsize=lbsz, xlabelsize=lbsz)
# time_vs_field(y, "Wsol_lm", line=true, fig=fig, fpos=(3,1), leg=false, linewidth=1.7, markersize=7, ylabelsize=lbsz, xlabelsize=lbsz)
# time_vs_field(y, "Water_lm", line=true, fig=fig, fpos=(3,2), leg=false, linewidth=1.7, markersize=7, ylabelsize=lbsz, xlabelsize=lbsz)

# mantle_water(y25)

# H2O_memory_time(y20, fig=fig, fpos=(1,1))
# H2O_memory_time(y25, fig=fig, fpos=(1,2), legend=false)
# H2O_memory_time(y30, fig=fig, fpos=(2,1), legend=false)
# H2O_memory_time(y35, fig=fig, fpos=(2,2), legend=false)
# time_vs_field(tot2, "satH2O_um")

# fig = Figure(size=(1000, 700))
# rprof_series(y50, "fmeltmean", [0.5, 1.5, 2.5, 3.5, 4.4], fig=fig, fpos=(1,1), logscale=false)
# rprof_series(y25, "Water", 1.5, fig=fig, fpos=(1,2))
# rprof_series(y25, "Water", 2.5, fig=fig, fpos=(2,1))
# rprof_series(y25, "Water", 3.5, fig=fig, fpos=(2,2))
# GLMakie.save("ywater.png", fig)

P = 0.5
T = 1300.
p = 0.0             
ϕ = 0.05
TOex = 0.05
TC = 0.1

niter=100

# # Dump oxygen in oxidizing C → CaCO₃?

# data    = Initialize_MAGEMin("sb24", verbose=false);
a, b = partition_Oₑₓ(P, T, p, ϕ, TOex, TC, debugging=true, verbose=true, plotevo=true, data=data, damp=damp, niter=niter);
# println("Done.")
# P_T_ϕ_TOₑₓ_TC_topoplogy(data, "T_ϕ", ns=5)
# Finalize_MAGEMin(data)