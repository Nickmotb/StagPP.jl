using StagPP, BenchmarkTools

# Stag, Dblock = aggregate_StagData(joinpath(@__DIR__, "test_data", "EW"), "y20e20");

rprof_vs_field(Dblock, "fO2"; fsize=(1200, 700), cmap=:blues)