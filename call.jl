using Revise, StagPP, BenchmarkTools, CairoMakie

# Dblock = load_sim(joinpath(@__DIR__, "test_data/EW"))
# idxT, idxR, idxP = data_encoding(y20e20)


a = solve_sH2O_fO2(40, 40);