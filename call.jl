using Revise, StagPP, CairoMakie

# Dblock = load_sim(joinpath(@__DIR__, "test_data/EW"))
# idxT, idxR, idxP = data_encoding(y20e20)

a = solve_sH2O_fO2(300, 300, DHMS=true);
plot_sᴴ²ᴼ(a, cmap=:davos, interp=false, cmap_reverse=true, logscale=false)
# minmap("um", "XH", ncols=3, nP=30, nT=30, savein="300x300mm_um")
# a = solve_point(17., 1500., "XH")

# min_s = min_sᴴ²ᴼ_assembler(50);