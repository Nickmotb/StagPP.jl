using Revise, StagPP, BenchmarkTools, CairoMakie

Dblock = load_sim(joinpath(@__DIR__, "test_data", "EW"), "y20e20";)
idxT, idxR, idxP = data_encoding(Dblock)

# fig = Figure(size = (1900, 800))
# time_vs_field(Dblock, "SurfOceanMass3D"; fig=fig, fpos=(1,1), color=:firebrick, mov_avg=true, savein="oceanmass", disp=false, tstart=0.1)
mantle_water(Dblock, 3.0)
# rprof_vs_field(Dblock, "Water"; fig=fig, fpos=(1,2), log=true, cmap=:vik100, cmap_reverse=true, Paxis=false)


# field_vs_field(Dblock, "OutgassedH2O"; fsize=(1000, 500), color=:firebrick, mov_avg=false, savein="ch2o", tstart=1.0, disp=true)
