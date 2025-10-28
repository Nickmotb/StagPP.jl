using Revise, StagPP, CairoMakie

# Dblock = load_sim(joinpath(@__DIR__, "test_data/EW"))
# idxT, idxR, idxP = data_encoding(y20e20)

a = solve_sH2O_fO2(150, 150; DHMS=true);
plot_sᴴ²ᴼ(a; cmap=:davos, interp=false, cmap_reverse=true, logscale=true)
# minmap("tz", "XH", ncols=3)
# a = solve_point(17., 1500., "XH")

# min_s = min_sᴴ²ᴼ_assembler(50);