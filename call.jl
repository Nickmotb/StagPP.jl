using Revise, StagPP, BenchmarkTools, CairoMakie

Dblock = load_sim(joinpath(@__DIR__, "test_data/EW"))
idxT, idxR, idxP = data_encoding(y20e20)
