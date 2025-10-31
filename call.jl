using Revise, StagPP, BenchmarkTools, CairoMakie

# Dblock = load_sim(joinpath("/home/nickmb/Desktop/StagYY/+op/EW"), "test")
# idxT, idxR, idxP = data_encoding(Dblock)
a = solve_sH2O_fO2(200, 200)
plot_sᴴ²ᴼ(a, cmap=:davos, cmap_reverse=true)

# rprof_vs_field(Dblock, "Tmean", cmap=:Blues, logscale=false)