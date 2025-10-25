using Revise, StagPP, BenchmarkTools, CairoMakie

# Dblock = load_sim(joinpath(@__DIR__, "test_data/EW"))
# idxT, idxR, idxP = data_encoding(y20e20)

a = solve_sH2O_fO2(150, 150);
plot_sᴴ²ᴼ(a; cmap=:davos, interp=true, cmap_reverse=true, logscale=true, savein="1")
# minmap("tz", "XB", savein="mm", ncols=3)
# a = solve_point(13., 3000., "XH")

# min_s = min_sᴴ²ᴼ_assembler(50);