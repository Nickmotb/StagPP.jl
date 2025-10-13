using Revise, StagPP, BenchmarkTools

Dblock = load_sim(joinpath(@__DIR__, "test_data", "EW"), "y20e20";)

# time_vs_field(Dblock, "Mob"; fsize=(1200, 700))
# rprof_vs_field(Dblock, "Water"; Paxis=true, log=true, cmap=:vik100, cmap_reverse=true)

# Try to modify reading routines such that parsers acts direclty on memory map, avoiding split-string