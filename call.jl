using StagPP, CairoMakie

load_sims = false
### Sims
    if load_sims
        # e20
        y20e20 = load_sim("/Volumes/Sinergia/EarthWater/y20e20/+op/EW", "y20e20");
        y25e20 = load_sim("/Volumes/Sinergia/EarthWater/y25e20/+op/EW", "y25e20");
        y30e20 = load_sim("/Volumes/Sinergia/EarthWater/y30e20/+op/EW", "y30e20");
        y35e20 = load_sim("/Volumes/Sinergia/EarthWater/y35e20/+op/EW", "y35e20");
        y40e20 = load_sim("/Volumes/Sinergia/EarthWater/y40e20/+op/EW", "y40e20");
        y45e20 = load_sim("/Volumes/Sinergia/EarthWater/y45e20/+op/EW", "y45e20");
        y50e20 = load_sim("/Volumes/Sinergia/EarthWater/y50e20/+op/EW", "y50e20");
        y55e20 = load_sim("/Volumes/Sinergia/EarthWater/y55e20/+op/EW", "y55e20");
        # e21
        y20e21 = load_sim("/Volumes/Sinergia/EarthWater/y20e21/+op/EW", "y20e21");
        y25e21 = load_sim("/Volumes/Sinergia/EarthWater/y25e21/+op/EW", "y25e21");
        y30e21 = load_sim("/Volumes/Sinergia/EarthWater/y30e21/+op/EW", "y30e21");
        y35e21 = load_sim("/Volumes/Sinergia/EarthWater/y35e21/+op/EW", "y35e21");
        y40e21 = load_sim("/Volumes/Sinergia/EarthWater/y40e21/+op/EW", "y40e21");
        y45e21 = load_sim("/Volumes/Sinergia/EarthWater/y45e21/+op/EW", "y45e21");
        y50e21 = load_sim("/Volumes/Sinergia/EarthWater/y50e21/+op/EW", "y50e21");
        y55e21 = load_sim("/Volumes/Sinergia/EarthWater/y55e21/+op/EW", "y55e21");
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
    end
###

fig = Figure(size = (2500, 1100))
ylabsz = 15
# rprof_vs_field(y55e20, "Water", fig=fig, fpos=(1,1), tstart=0.1, logscale=true, cmap_reverse=true, disp=false, ylabelsize=ylabsz)
# rprof_vs_field(y45e20, "Water", fig=fig, fpos=(1,3), tstart=0.1, logscale=true, cmap_reverse=true, disp=false, ylabelsize=ylabsz)
# rprof_vs_field(y35e20, "Water", fig=fig, fpos=(1,5), tstart=0.1, logscale=true, cmap_reverse=true, disp=false, ylabelsize=ylabsz)
# rprof_vs_field(y30e20, "Water", fig=fig, fpos=(1,7), tstart=0.1, logscale=true, cmap_reverse=true, disp=false, ylabelsize=ylabsz)

# rprof_vs_field(y55e20, "Wsol", fig=fig, fpos=(2,1), tstart=0.1, logscale=true, cmap=:greys, cmap_reverse=true, disp=false, ylabelsize=ylabsz)
# rprof_vs_field(y45e20, "Wsol", fig=fig, fpos=(2,3), tstart=0.1, logscale=true, cmap=:greys, cmap_reverse=true, disp=false, ylabelsize=ylabsz)
# rprof_vs_field(y35e20, "Wsol", fig=fig, fpos=(2,5), tstart=0.1, logscale=true, cmap=:greys, cmap_reverse=true, disp=false, ylabelsize=ylabsz)
# rprof_vs_field(y30e20, "Wsol", fig=fig, fpos=(2,7), tstart=0.1, logscale=true, cmap=:greys, cmap_reverse=true, disp=false, ylabelsize=ylabsz)

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

# time_vs_field(e20, "cH2O_mean", fig=fig, fpos=(3,3:4), tstart=0.1, tend=3.0, ylabelsize=ylabsz, mov_avg=false)
# ta_field_vs_field([e20, e21], "Mob", "cH2O_mean", fig=fig, fpos=(3,1:2), tstart=0.1, line=false, ylabelsize=ylabsz, marker=:rect, strokewidth=1.5, markersize=25,setlabs=[L"\eta_0 = 10^{20}", L"\eta_0 = 10^{21}"])

#field_vs_field(y20e20, "Wsol_tz", "Tmean")
# mantle_water_at_t(tot3, 0.0)

# Rprofiles for H2O and Wsol.
# ta correlation between mob and cH2O_mean, and Ingas/Outgas fluxes
# ACF -> e-folding timescales memory vs reversals + integrated memory times

H2O_memory_time(y30e20)