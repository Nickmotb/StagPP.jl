using StagPP, GLMakie

load_sims = true
### Sims
    if load_sims
        # # e20 old
        # y20e20o = load_sim("/Volumes/Sinergia/EarthWater/old/HighPerm_HighS/y20e20/+op/EW", "γ20η20");
        # y25e20o = load_sim("/Volumes/Sinergia/EarthWater/old/HighPerm_HighS/y25e20/+op/EW", "γ25η20");
        # y30e20o = load_sim("/Volumes/Sinergia/EarthWater/old/HighPerm_HighS/y30e20/+op/EW", "γ30η20");
        # y35e20o = load_sim("/Volumes/Sinergia/EarthWater/old/HighPerm_HighS/y35e20/+op/EW", "γ35η20");
        # y40e20o = load_sim("/Volumes/Sinergia/EarthWater/old/HighPerm_HighS/y40e20/+op/EW", "γ40η20");
        # y45e20o = load_sim("/Volumes/Sinergia/EarthWater/old/HighPerm_HighS/y45e20/+op/EW", "γ45η20");
        # y50e20o = load_sim("/Volumes/Sinergia/EarthWater/old/HighPerm_HighS/y50e20/+op/EW", "γ50η20");
        # # e21 old
        # y20e21o = load_sim("/Volumes/Sinergia/EarthWater/old/HighPerm_HighS/y20e21/+op/EW", "γ20η21");
        # y25e21o = load_sim("/Volumes/Sinergia/EarthWater/old/HighPerm_HighS/y25e21/+op/EW", "γ25η21");
        # y30e21o = load_sim("/Volumes/Sinergia/EarthWater/old/HighPerm_HighS/y30e21/+op/EW", "γ30η21");
        # y35e21o = load_sim("/Volumes/Sinergia/EarthWater/old/HighPerm_HighS/y35e21/+op/EW", "γ35η21");
        # y40e21o = load_sim("/Volumes/Sinergia/EarthWater/old/HighPerm_HighS/y40e21/+op/EW", "γ40η21");
        # y45e21o = load_sim("/Volumes/Sinergia/EarthWater/old/HighPerm_HighS/y45e21/+op/EW", "γ45η21");
        # y50e21o = load_sim("/Volumes/Sinergia/EarthWater/old/HighPerm_HighS/y50e21/+op/EW", "γ50η21");

        # e20 new
        y20e20 = load_sim("/Volumes/Sinergia/EarthWater/y20e20/+op/EW", "γ20η20");
        y25e20 = load_sim("/Volumes/Sinergia/EarthWater/y25e20/+op/EW", "γ25η20");
        y30e20 = load_sim("/Volumes/Sinergia/EarthWater/y30e20/+op/EW", "γ30η20");
        y35e20 = load_sim("/Volumes/Sinergia/EarthWater/y35e20/+op/EW", "γ35η20");
        y40e20 = load_sim("/Volumes/Sinergia/EarthWater/y40e20/+op/EW", "γ40η20");
        y45e20 = load_sim("/Volumes/Sinergia/EarthWater/y45e20/+op/EW", "γ45η20");
        y50e20 = load_sim("/Volumes/Sinergia/EarthWater/y50e20/+op/EW", "γ50η20");
        y55e20 = load_sim("/Volumes/Sinergia/EarthWater/y55e20/+op/EW", "γ55η20");
        # e21 new
        y20e21 = load_sim("/Volumes/Sinergia/EarthWater/y20e21/+op/EW", "γ20η21");
        y25e21 = load_sim("/Volumes/Sinergia/EarthWater/y25e21/+op/EW", "γ25η21");
        y30e21 = load_sim("/Volumes/Sinergia/EarthWater/y30e21/+op/EW", "γ30η21");
        y35e21 = load_sim("/Volumes/Sinergia/EarthWater/y35e21/+op/EW", "γ35η21");
        y40e21 = load_sim("/Volumes/Sinergia/EarthWater/y40e21/+op/EW", "γ40η21");
        y45e21 = load_sim("/Volumes/Sinergia/EarthWater/y45e21/+op/EW", "γ45η21");
        y50e21 = load_sim("/Volumes/Sinergia/EarthWater/y50e21/+op/EW", "γ50η21");
        y55e21 = load_sim("/Volumes/Sinergia/EarthWater/y55e21/+op/EW", "γ55η21");
        # id
        id4 = load_sim("/Volumes/Sinergia/EarthWater/id4/+op/EW", "id4");
        id6 = load_sim("/Volumes/Sinergia/EarthWater/id6/+op/EW", "id6");
        id8 = load_sim("/Volumes/Sinergia/EarthWater/id8/+op/EW", "id8");
        # tot
        tot1 = load_sim("/Volumes/Sinergia/EarthWater/tot1/+op/EW", "tot1");
        tot2 = load_sim("/Volumes/Sinergia/EarthWater/tot2/+op/EW", "tot2");
        tot3 = load_sim("/Volumes/Sinergia/EarthWater/tot3/+op/EW", "tot3");
        tot4 = load_sim("/Volumes/Sinergia/EarthWater/tot4/+op/EW", "tot4");
        tot5 = load_sim("/Volumes/Sinergia/EarthWater/tot5/+op/EW", "tot5");
        # no
        noK = load_sim("/Volumes/Sinergia/EarthWater/noK/+op/EW", "noK");
        noR = load_sim("/Volumes/Sinergia/EarthWater/noR/+op/EW", "noR");
        noKR = load_sim("/Volumes/Sinergia/EarthWater/noKR/+op/EW", "noKR");
        # sets
        e20 = [y20e20, y25e20, y30e20, y35e20, y40e20, y45e20, y50e20, y55e20];
        e21 = [y20e21, y25e21, y30e21, y35e21, y40e21, y45e21, y50e21, y55e21];
        id  = [id4, id6, id8, y20e20];
        tot = [tot1, tot2, tot3, tot4, tot5, y20e20];
        no  = [y20e20, noK, noR, noKR];
        # e20o = [y20e20o, y25e20o, y30e20o, y35e20o, y40e20o, y45e20o, y50e20o];
        # e21o = [y20e21o, y25e21o, y30e21o, y35e21o, y40e21o, y45e21o, y50e21o];
    end
###

fig = Figure(size = (1600, 700))
xlabsz, ylabsz = 15, 15
colorrange=(0.0, 100)
setlabs = [L"\eta_0 = 10^{20}\;set", L"\eta_0 = 10^{21}\;set", L"Total\;H_2O\;set", L"d_{depth}\;set", L"No\;rheo\;set"]
# rprof_vs_field(y50e20, "satH2O", fig=fig, fpos=(1,1), tstart=0.1, logscale=false, cmap_reverse=true, disp=false, xlabelsize=xlabsz, ylabelsize=ylabsz, colorrange=colorrange, title=y55e20.metadata.name)
# rprof_vs_field(y40e20, "satH2O", fig=fig, fpos=(1,3), tstart=0.1, logscale=false, cmap_reverse=true, disp=false, xlabelsize=xlabsz, ylabelsize=ylabsz, colorrange=colorrange, title=y40e20.metadata.name)
# rprof_vs_field(y35e20, "satH2O", fig=fig, fpos=(1,5), tstart=0.1, logscale=false, cmap_reverse=true, disp=false, xlabelsize=xlabsz, ylabelsize=ylabsz, colorrange=colorrange, title=y35e20.metadata.name)
# rprof_vs_field(y30e20, "satH2O", fig=fig, fpos=(1,7), tstart=0.1, logscale=false, cmap_reverse=true, disp=false, xlabelsize=xlabsz, ylabelsize=ylabsz, colorrange=colorrange, title=y30e20.metadata.name)
# rprof_vs_field(y25e20, "satH2O", fig=fig, fpos=(1,9), tstart=0.1, logscale=false, cmap_reverse=true, disp=false, xlabelsize=xlabsz, ylabelsize=ylabsz, colorrange=colorrange, title=y25e20.metadata.name)
# rprof_vs_field(y20e20, "satH2O", fig=fig, fpos=(1,11), tstart=0.1, logscale=false, cmap_reverse=true, disp=false, xlabelsize=xlabsz, ylabelsize=ylabsz, colorrange=colorrange, title=y25e20.metadata.name)

# rprof_vs_field(y50e20, "Wsol", fig=fig, fpos=(2,1), tstart=0.1, logscale=true, cmap=:greys, cmap_reverse=true, disp=false, ylabelsize=ylabsz, xlabelsize=xlabsz)
# rprof_vs_field(y40e20, "Wsol", fig=fig, fpos=(2,3), tstart=0.1, logscale=true, cmap=:greys, cmap_reverse=true, disp=false, ylabelsize=ylabsz, xlabelsize=xlabsz)
# rprof_vs_field(y35e20, "Wsol", fig=fig, fpos=(2,5), tstart=0.1, logscale=true, cmap=:greys, cmap_reverse=true, disp=false, ylabelsize=ylabsz, xlabelsize=xlabsz)
# rprof_vs_field(y30e20, "Wsol", fig=fig, fpos=(2,7), tstart=0.1, logscale=true, cmap=:greys, cmap_reverse=true, disp=false, ylabelsize=ylabsz, xlabelsize=xlabsz)
# rprof_vs_field(y25e20, "Wsol", fig=fig, fpos=(2,9), tstart=0.1, logscale=true, cmap=:greys, cmap_reverse=true, disp=false, ylabelsize=ylabsz, xlabelsize=xlabsz)
# rprof_vs_field(y20e20, "Wsol", fig=fig, fpos=(2,11), tstart=0.1, logscale=true, cmap=:greys, cmap_reverse=true, disp=true, ylabelsize=ylabsz, xlabelsize=xlabsz)

# time_vs_field(e20, "Vrms_lm", mov\_avg=true, logscale=false, line=false)
# ta_field_vs_field([e20, e21, tot, id, no], "Mob", "Water", fig=fig, fpos=(1,2), tstart=0.25, line=false, ylabelsize=ylabsz, strokewidth=1.5, markersize=15,marker=:circle, 
#                     setlabs=[L"\eta_0 = 10^{20}\;set", L"\eta_0 = 10^{21}\;set", L"Total\;H_2O\;set", L"d_{depth}\;set", L"No\;rheo\;set"])
# ta_field_vs_field([e20, e21], "Mob", "fmeltmean_um", fig=fig, fpos=(3,5:10), tstart=0.2, line=false, ylabelsize=ylabsz, strokewidth=1.5, markersize=15,marker=:rect)

# mantle_water(y55e20)

# rprof_vs_field(tot1, "Water", fig=fig, fpos=(1,1), tstart=0.1, logscale=true, cmap_reverse=true, disp=false, ylabelsize=ylabsz)
# rprof_vs_field(tot2, "Water", fig=fig, fpos=(1,3), tstart=0.1, logscale=true, cmap_reverse=true, disp=false, ylabelsize=ylabsz)
# rprof_vs_field(tot3, "Water", fig=fig, fpos=(1,5), tstart=0.1, logscale=true, cmap_reverse=true, disp=false, ylabelsize=ylabsz)
# rprof_vs_field(tot4, "Water", fig=fig, fpos=(1,7), tstart=0.1, logscale=true, cmap_reverse=true, disp=false, ylabelsize=ylabsz)
# rprof_vs_field(tot5, "Water", fig=fig, fpos=(1,9), tstart=0.1, logscale=true, cmap_reverse=true, disp=false, ylabelsize=ylabsz)

# rprof_vs_field(tot1, "Wsol", fig=fig, fpos=(2,1), tstart=0.1, logscale=true, cmap=:greys, cmap_reverse=true, disp=false, ylabelsize=ylabsz)
# rprof_vs_field(tot2, "Wsol", fig=fig, fpos=(2,3), tstart=0.1, logscale=true, cmap=:greys, cmap_reverse=true, disp=false, ylabelsize=ylabsz)
# rprof_vs_field(tot3, "Wsol", fig=fig, fpos=(2,5), tstart=0.1, logscale=true, cmap=:greys, cmap_reverse=true, disp=false, ylabelsize=ylabsz)
# rprof_vs_field(tot4, "Wsol", fig=fig, fpos=(2,7), tstart=0.1, logscale=true, cmap=:greys, cmap_reverse=true, disp=false, ylabelsize=ylabsz)
# rprof_vs_field(tot5, "Wsol", fig=fig, fpos=(2,9), tstart=0.1, logscale=true, cmap=:greys, cmap_reverse=true, disp=true, ylabelsize=ylabsz)

# time_vs_field(no, "Water_tz", fig=fig, fpos=(1,1), mov_avg=true)
# snapshot(y20e20, 100, "C:Water")



# ta_field_vs_field(e20, "Mob", "SurfOceanMass3D", line=false)
# mantle_water_at_t(tot3, 0.0)

# Rprofiles for H2O and Wsol.
# ta correlation between mob and cH2O_mean, and Ingas/Outgas fluxes
# ACF -> e-folding timescales memory vs reversals + integrated memory times

loadings = H2O_PCA(e20, "twindow", "ys", fig=fig, fpos=(1,1), disp=true)
# loadings = H2O_PCA(e20, "twindow", "t", fig=fig, fpos=(1,1), disp=true)
#  H2O_memory_time(e21, fig=fig, fpos=(1,1), disp=true)
# mantle_water_at_t(y20e20, 4.5)


