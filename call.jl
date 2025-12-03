using StagPP

# Dblock = load_sim("/Volumes/Sinergia/EarthWater/y20e20/+op/EW", "test")
Dblock = load_sim("/home/nickmb/Desktop/StagYY/+op/EW", "test")
idxT, idxR, idxP = data_encoding(Dblock)

# fig = Figure(size=(1900,1200))
# rprof_vs_field(Dblock, "Water", cmap=:vik100, cmap_reverse=true, logscale=true, fig=fig, fpos=(1,1), interpolate=true, disp=false)
# rprof_vs_field(Dblock, "Wsol", cmap=:grays, cmap_reverse=true, logscale=true, fig=fig, fpos=(1,3), interpolate=true, disp=false)
# mantle_water(Dblock, Dblock.timedata[end,2], fig=fig, fpos=(2,3), disp=false)
snapshot(Dblock, 100, "C:Water", logscale=true)
# IOplot(Dblock)

