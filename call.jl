using StagPP, BenchmarkTools

Stag, Dblock = aggregate_StagData(joinpath(@__DIR__, "test_data", "EW"), "y20e20");

time_vs_field(Dblock, "cH2O_mean"; scatter=true)