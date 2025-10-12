using StagPP, BenchmarkTools

path = joinpath(@__DIR__, "test_data", "EW")
Stag, Dblock = aggregate_StagData(path, "y20e20");
