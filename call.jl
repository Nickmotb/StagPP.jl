using Revise, StagPP, BenchmarkTools

# Dblock = load_sim(joinpath(@__DIR__, "test_data", "EW"), "y20e20";)

time_vs_field(Dblock, "cH2O_mean"; fsize=(1000, 500), color=:blue, mov_avg=false, savein="ch2o", tstart=1.0)
# rprof_vs_field(Dblock, "Water"; log=true, cmap=:vik100, cmap_reverse=true)
