using Revise, StagPP, CairoMakie, BenchmarkTools

# Dblock = load_sim(joinpath(@__DIR__, "test_data/EW"))
# idxT, idxR, idxP = data_encoding(y20e20)

a = solve_sH2O_fO2(50, 50, DHMS=true);
plot_sᴴ²ᴼ(a, cmap=:davos, interp=false, cmap_reverse=true, logscale=true)
# minmap("um", "XH", ncols=3, nP=100, nT=100, savein="100x100mm_um")
# a = solve_point(17., 1500., "XH")

# min_s = min_sᴴ²ᴼ_assembler(50);