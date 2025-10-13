using Revise, StagPP, BenchmarkTools

Stag, Dblock = load_sim(joinpath(@__DIR__, "test_data", "EW"), "y20e20");

rprof_vs_field(Dblock, "bsmean"; fsize=(1200, 700), cmap=:vik, log=false)

# Read also plates. time vector is the same as rprof and contains mobility.