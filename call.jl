using Revise, StagPP, BenchmarkTools, CairoMakie

# Dblock = load_sim(joinpath(@__DIR__, "test_data/EW"))
# idxT, idxR, idxP = data_encoding(y20e20)

a = solve_sH2O_fO2(50, 50; plt=true, interp=true, cmap=:davos, cmap_reverse=true);
# plot_sᴴ²ᴼ(a; cmap=:davos, interp=true, cmap_reverse=true, logscale=true)