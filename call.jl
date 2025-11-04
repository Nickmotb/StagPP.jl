using Revise, StagPP, BenchmarkTools, CairoMakie

Dblock = load_sim(joinpath("/home/nickmb/Desktop/StagYY/+op/EW"), "test")
idxT, idxR, idxP = data_encoding(Dblock)
# a = solve_sH2O_fO2(200, 200)
# plot_sᴴ²ᴼ(a, cmap=:davos, cmap_reverse=true)

fig = Figure(size=(1900,1200))
rprof_vs_field(Dblock, "Water", cmap=:vik100, cmap_reverse=true, logscale=true, fig=fig, fpos=(1,1), interpolate=true, disp=false, tstart=0.1)
rprof_vs_field(Dblock, "Wsol", cmap=:grays, cmap_reverse=true, logscale=true, fig=fig, fpos=(1,3), interpolate=true, disp=false, tstart=0.1)
mantle_water(Dblock, Dblock.timedata[end,2], fig=fig, fpos=(2,3), disp=false)
snapshot(Dblock, 100, "Water", fig=fig, fpos=(2,1), logscale=true)

