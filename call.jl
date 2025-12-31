using StagPP, CairoMakie

y40e21 = load_sim("/Volumes/Sinergia/EarthWater/y40e21/+op/EW", "test")
# y30e21 = load_sim("/Volumes/Sinergia/EarthWater/y30e21/+op/EW", "test2")
# y40e21 = load_sim("/Volumes/Sinergia/EarthWater/y40e21/+op/EW", "test3")
# y60e21 = load_sim("/Volumes/Sinergia/EarthWater/y60e21/+op/EW", "test5")

# fig = Figure(size=(2200,1000))
# rprof_vs_field(y20e21, "H2Ofree", cmap=:vik100, cmap_reverse=true, logscale=true, fig=fig, fpos=(1,1), interpolate=true, disp=false, colorrange=(-5, 1))
# rprof_vs_field(y30e21, "H2Ofree", cmap=:vik100, cmap_reverse=true, logscale=true, fig=fig, fpos=(1,3), interpolate=true, disp=false, colorrange=(-5, 1))
# rprof_vs_field(y40e21, "H2Ofree", cmap=:vik100, cmap_reverse=true, logscale=true, fig=fig, fpos=(1,5), interpolate=true, disp=false, colorrange=(-5, 1))
# rprof_vs_field(y60e21, "H2Ofree", cmap=:vik100, cmap_reverse=true, logscale=true, fig=fig, fpos=(1,7), interpolate=true, disp=false)

# rprof_vs_field(y20e21, "bsmean", cmap=:RdGy, cmap_reverse=true, logscale=false, fig=fig, fpos=(2,1), interpolate=true, disp=false)
# rprof_vs_field(y30e21, "bsmean", cmap=:RdGy, cmap_reverse=true, logscale=false, fig=fig, fpos=(2,3), interpolate=true, disp=false)
# rprof_vs_field(y40e21, "bsmean", cmap=:RdGy, cmap_reverse=true, logscale=false, fig=fig, fpos=(2,5), interpolate=true, disp=false)
# rprof_vs_field(y60e21, "bsmean", cmap=:RdGy, cmap_reverse=true, logscale=false, fig=fig, fpos=(2,7), interpolate=true, disp=false)

# mantle_water(y20e21, y20e21.metadata.tend, fig=fig, fpos=(3,1:2), disp=false)
# mantle_water(y30e21, y30e21.metadata.tend, fig=fig, fpos=(3,3:4), disp=false)
# mantle_water(y40e21, y40e21.metadata.tend, fig=fig, fpos=(3,5:6), disp=false)
# mantle_water(y60e21, y60e21.metadata.tend, fig=fig, fpos=(3,7:8), disp=true, savein="e21_Hfree_bsmean")

# snapshot(Dblock, 100, "C:Water", logscale=true, fig=fig, fpos=(2,1))

IOplot(y40e21)
