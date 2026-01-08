using StagPP, CairoMakie

# y20e20 = load_sim("/Volumes/Sinergia/EarthWater/old/HighPerm_HighS/y20e20/+op/EW", "y20e20");
# y25e20 = load_sim("/Volumes/Sinergia/EarthWater/old/HighPerm_HighS/y25e20/+op/EW", "y25e20");
# y30e20 = load_sim("/Volumes/Sinergia/EarthWater/old/HighPerm_HighS/y30e20/+op/EW", "y30e20");
# y35e20 = load_sim("/Volumes/Sinergia/EarthWater/old/HighPerm_HighS/y35e20/+op/EW", "y35e20");
# y40e20 = load_sim("/Volumes/Sinergia/EarthWater/old/HighPerm_HighS/y40e20/+op/EW", "y40e20");
# y45e20 = load_sim("/Volumes/Sinergia/EarthWater/old/HighPerm_HighS/y45e21/+op/EW", "y45e20");
# y50e20 = load_sim("/Volumes/Sinergia/EarthWater/old/HighPerm_HighS/y50e21/+op/EW", "y50e20");

# y20e21 = load_sim("/Volumes/Sinergia/EarthWater/old/HighPerm_HighS/y20e21/+op/EW", "y20e21");
# y25e21 = load_sim("/Volumes/Sinergia/EarthWater/old/HighPerm_HighS/y25e21/+op/EW", "y25e21");
# y30e21 = load_sim("/Volumes/Sinergia/EarthWater/old/HighPerm_HighS/y30e21/+op/EW", "y30e21");
# y35e21 = load_sim("/Volumes/Sinergia/EarthWater/old/HighPerm_HighS/y35e21/+op/EW", "y35e21");
# y40e21 = load_sim("/Volumes/Sinergia/EarthWater/old/HighPerm_HighS/y40e21/+op/EW", "y40e21");
# y45e21 = load_sim("/Volumes/Sinergia/EarthWater/old/HighPerm_HighS/y45e21/+op/EW", "y45e21");
# y50e21 = load_sim("/Volumes/Sinergia/EarthWater/old/HighPerm_HighS/y50e21/+op/EW", "y50e21");

# id4 = load_sim("/Volumes/Sinergia/EarthWater/id4/+op/EW", "id4");

# dblocks20 = [y20e20, y25e20, y30e20, y35e20, y40e20, y45e20, y50e20];
# dblocks21 = [y20e21, y25e21, y30e21, y35e21, y40e21, y45e21, y50e21];
# blocks = [dblocks20, dblocks21];

fig = Figure(size = (2200, 1100))
ylabsz = 18
rprof_vs_field(y50e20, "Water", fig=fig, fpos=(1,1), tstart=0.5, logscale=true, cmap_reverse=true, disp=false, ylabelsize=ylabsz)
rprof_vs_field(y40e20, "Water", fig=fig, fpos=(1,3), tstart=0.5, logscale=true, cmap_reverse=true, disp=false, ylabelsize=ylabsz)
rprof_vs_field(y30e20, "Water", fig=fig, fpos=(1,5), tstart=0.5, logscale=true, cmap_reverse=true, disp=false, ylabelsize=ylabsz)
rprof_vs_field(y20e20, "Water", fig=fig, fpos=(1,7), tstart=0.5, logscale=true, cmap_reverse=true, disp=false, ylabelsize=ylabsz)

mantle_water_at_t(y50e20, 4.5, fig=fig, fpos=(2,1:2), disp=false, ylabelsize=ylabsz)
mantle_water_at_t(y40e20, 4.5, fig=fig, fpos=(2,3:4), disp=false, ylabelsize=ylabsz)
mantle_water_at_t(y30e20, 4.5, fig=fig, fpos=(2,5:6), disp=false, ylabelsize=ylabsz)
mantle_water_at_t(y20e20, 4.5, fig=fig, fpos=(2,7:8), disp=false, ylabelsize=ylabsz)

ta_field_vs_field(blocks, "Mob", "cH2O_mean", fig=fig, fpos=(3,1:2), tstart=0.5, line=false, ylabelsize=ylabsz)
ta_field_vs_field(blocks, "Mob", "SurfOceanMass3D", fig=fig, fpos=(3,3:4), tstart=0.5, line=false, ylabelsize=ylabsz)
ta_field_vs_field(blocks, "Mob", "e/eH2O_norm", fig=fig, fpos=(3,5:6), tstart=0.5, line=false, ylabelsize=ylabsz)
ta_field_vs_field(blocks, "Mob", "erupta", fig=fig, fpos=(3,7:8), tstart=0.5, line=false, ylabelsize=ylabsz)

# mantle_water(y50e20, fig=fig, fpos=(2,1:2), tstart=0.5, logscale=true, disp=false)
# mantle_water(y40e20, fig=fig, fpos=(2,3:4), tstart=0.5, logscale=true, disp=false)
# mantle_water(y30e20, fig=fig, fpos=(2,5:6), tstart=0.5, logscale=true, disp=false)
# mantle_water(y20e20, fig=fig, fpos=(2,7:8), tstart=0.5, logscale=true)

