using Revise, StagPP, BenchmarkTools, CairoMakie

# Dblock = load_sim(joinpath("/home/nickmb/Desktop/StagYY/+op/EW"), "test")
# idxT, idxR, idxP = data_encoding(Dblock)
# a = solve_sH2O_fO2(200, 200)
# plot_sᴴ²ᴼ(a, cmap=:davos, cmap_reverse=true)

# fig = Figure(size=(1900,900))
# rprof_vs_field(Dblock, "Water", cmap=:Blues, logscale=true, fig=fig, fpos=(1,1), interpolate=true, disp=false)
# rprof_vs_field(Dblock, "Wsol", cmap=:grays, cmap_reverse=true, logscale=true, fig=fig, fpos=(1,3), interpolate=true, disp=false)
# rprof_vs_field(Dblock, "Tmean", cmap=:vik100, logscale=true, fig=fig, fpos=(2,1), interpolate=true, disp=false)
# mantle_water(Dblock, Dblock.timedata[end,2], fig=fig, fpos=(2,3))


@btime readVTK("/home/nickmb/Desktop/StagYY/+op/EW_ALL00000.vts")