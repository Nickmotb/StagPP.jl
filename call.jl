using StagPP, GLMakie

load_sims = false
# Sims
    if load_sims
        # Full Runs
        y20 = load_sim("/Volumes/Sinergia/EarthWater/y20e19/+op/EW", "γ20");
        y25 = load_sim("/Volumes/Sinergia/EarthWater/y25e19/+op/EW", "γ25");
        y30 = load_sim("/Volumes/Sinergia/EarthWater/y30e19/+op/EW", "γ30");
        y35 = load_sim("/Volumes/Sinergia/EarthWater/y35e19/+op/EW", "γ35");
        y40 = load_sim("/Volumes/Sinergia/EarthWater/y40e19/+op/EW", "γ40");
        y45 = load_sim("/Volumes/Sinergia/EarthWater/y45e19/+op/EW", "γ45");
        y50 = load_sim("/Volumes/Sinergia/EarthWater/y50e19/+op/EW", "γ50");
        y55 = load_sim("/Volumes/Sinergia/EarthWater/y55e19/+op/EW", "γ55");

        # # No Water
        # y20nw = load_sim("/Volumes/Sinergia/EarthWater/y20e19nw/+op/EW", "γ20nw");
        # y25nw = load_sim("/Volumes/Sinergia/EarthWater/y25e19nw/+op/EW", "γ25nw");
        # y30nw = load_sim("/Volumes/Sinergia/EarthWater/y30e19nw/+op/EW", "γ30nw");
        # y35nw = load_sim("/Volumes/Sinergia/EarthWater/y35e19nw/+op/EW", "γ35nw");
        # y40nw = load_sim("/Volumes/Sinergia/EarthWater/y40e19nw/+op/EW", "γ40nw");
        # y45nw = load_sim("/Volumes/Sinergia/EarthWater/y45e19nw/+op/EW", "γ45nw");
        # y50nw = load_sim("/Volumes/Sinergia/EarthWater/y50e19nw/+op/EW", "γ50nw");
        # y55nw = load_sim("/Volumes/Sinergia/EarthWater/y55e19nw/+op/EW", "γ55nw");

        logm05 = load_sim("/Volumes/Sinergia/EarthWater/darcy-0.5logm/+op/EW", "κ₀-0.5log");
        logp05 = load_sim("/Volumes/Sinergia/EarthWater/darcy-0.5logp/+op/EW", "κ₀+0.5log");
        logm = load_sim("/Volumes/Sinergia/EarthWater/darcy-logm/+op/EW", "κ₀-log");
        logp = load_sim("/Volumes/Sinergia/EarthWater/darcy-logp/+op/EW", "κ₀+log");
        d1_16 = load_sim("/Volumes/Sinergia/EarthWater/darcy1-16/+op/EW", "um_κ₀=1×10⁻¹⁶");
        d3_16 = load_sim("/Volumes/Sinergia/EarthWater/darcy3-16/+op/EW", "um_κ₀=3×10⁻¹⁶");
        d7_16 = load_sim("/Volumes/Sinergia/EarthWater/darcy7-16/+op/EW", "um_κ₀=7×10⁻¹⁶");
        d9_16 = load_sim("/Volumes/Sinergia/EarthWater/darcy9-16/+op/EW", "um_κ₀=9×10⁻¹⁶");
        d5_17 = load_sim("/Volumes/Sinergia/EarthWater/darcy5-17/+op/EW", "um_κ₀=5×10⁻¹⁷");
        d7_17 = load_sim("/Volumes/Sinergia/EarthWater/darcy7-17/+op/EW", "um_κ₀=7×10⁻¹⁷");
        d9_17 = load_sim("/Volumes/Sinergia/EarthWater/darcy9-17/+op/EW", "um_κ₀=9×10⁻¹⁷");
        d2_15 = load_sim("/Volumes/Sinergia/EarthWater/darcy2-15/+op/EW", "um_κ₀=2×10⁻¹⁵");
        d4_15 = load_sim("/Volumes/Sinergia/EarthWater/darcy4-15/+op/EW", "um_κ₀=4×10⁻¹⁵");
        d6_15 = load_sim("/Volumes/Sinergia/EarthWater/darcy6-15/+op/EW", "um_κ₀=6×10⁻¹⁵");

        # Ingas depth
        id4 = load_sim("/Volumes/Sinergia/EarthWater/id4/+op/EW", "id4");
        id6 = load_sim("/Volumes/Sinergia/EarthWater/id6/+op/EW", "id6");
        id8 = load_sim("/Volumes/Sinergia/EarthWater/id8/+op/EW", "id8");
        
        # No RK
        noR = load_sim("/Volumes/Sinergia/EarthWater/noR/+op/EW", "noR");
        noK = load_sim("/Volumes/Sinergia/EarthWater/noK/+op/EW", "No Katz");
        # noKR = load_sim("/Volumes/Sinergia/EarthWater/noKR/+op/EW", "noKR");

        # Total Water Budget
        tot1 = load_sim("/Volumes/Sinergia/EarthWater/tot1/+op/EW", "1 OM");
        tot2 = load_sim("/Volumes/Sinergia/EarthWater/tot2/+op/EW", "2 OM");
        tot3 = load_sim("/Volumes/Sinergia/EarthWater/tot3/+op/EW", "3 OM");
        tot4 = load_sim("/Volumes/Sinergia/EarthWater/tot4/+op/EW", "4 OM");
        tot5 = load_sim("/Volumes/Sinergia/EarthWater/tot5/+op/EW", "5 OM");

        # # sets
        y = [y20, y25, y30, y35, y40, y45, y50, y55];
        fulldarcy = [logm, logm05, logp05, logp]
        umdarcy = [d5_17, d7_17, d9_17, d1_16, d3_16, d7_16, d9_16, d2_15, d4_15, d6_15]
        # ynw = [y20nw, y25nw, y30nw, y35nw, y40nw, y45nw, y50nw, y55nw];
        # yt = vcat(y, ynw)
        id = [id4, id6, id8, y20];
        tot = [tot1, tot2, tot3, tot4, tot5, y20];
        # no = [y20, noR, noK, noKR]
    
    end

fig = Figure(size=(1000,600))
# mantle_water(y20, tstart=0.2, tend=1.8, fluxes=true, Tg_Myr=true, leg=false, fig=fig, fpos=(1,1), disp=false, mavg=true)
# mantle_water(y25, tstart=0.2, tend=1.8, fluxes=true, Tg_Myr=true, leg=false, fig=fig, fpos=(2,1), disp=false, mavg=true)
# mantle_water(y30, tstart=0.2, tend=1.8, fluxes=true, Tg_Myr=true, leg=false, fig=fig, fpos=(3,1), mavg=false)

# Fig 1
# ts=0.2; yr=0.75; y1=(-5.0e19, 5.0e19)
# IOplot(y20, fig=fig, fpos=(1,1:2), D2D3=true, yrange1=y1, yrange2=yr, tstart=ts,xticklabelsize=15,yticklabelsize=15,xlabelsize=18,ylabelsize=18, disp=false)
# IOplot(y25, fig=fig, fpos=(2,1:2), D2D3=true, yrange1=y1, yrange2=yr, tstart=ts,xticklabelsize=15,yticklabelsize=15,xlabelsize=18,ylabelsize=18, disp=false)
# IOplot(y30, fig=fig, fpos=(3,1:2), D2D3=true, yrange1=y1, yrange2=yr, tstart=ts,xticklabelsize=15,yticklabelsize=15,xlabelsize=18,ylabelsize=18, disp=false)
# IOplot(y35, fig=fig, fpos=(4,1:2), D2D3=true, yrange1=y1, yrange2=yr, tstart=ts, savein="IOplots", disp=false,xticklabelsize=15,yticklabelsize=15,xlabelsize=18,ylabelsize=18)
# time_vs_field(y20, "influxH2O", fig=fig, fpos=(1,3), tstart=ts, disp=false, mov_avg=true, color=:black, logy=true, yrange=(5e7, 5e10), mavg_contrast=0.05)
# time_vs_field(y25, "influxH2O", fig=fig, fpos=(2,3), tstart=ts, disp=false, mov_avg=true, color=:black, logy=true, yrange=(5e7, 5e10), mavg_contrast=0.05)
# time_vs_field(y30, "influxH2O", fig=fig, fpos=(3,3), tstart=ts, disp=false, mov_avg=true, color=:black, logy=true, mavg_contrast=0.05, yrange=(5e7, 5e11))
# time_vs_field(y35, "influxH2O", fig=fig, fpos=(4,3), tstart=ts, mov_avg=true, color=:black, logy=true, mavg_contrast=0.05, yrange=(5e7, 5e11), savein="IOplots")


# Fig 2
# crange=(-3.5, 0.25); itp=true; ts=0.5
# rprof_vs_field(y20, "Water", fig=fig, fpos=(1,1), cmap_reverse=true, logy=true, colorrange=crange, interpolate=itp, tstart=ts, xlabelvisible=false,xticklabelsize=15,yticklabelsize=15,xlabelsize=18,ylabelsize=18)
# rprof_vs_field(y25, "Water", fig=fig, fpos=(3,1), cmap_reverse=true, logy=true, colorrange=crange, cbar=false, interpolate=itp, tstart=ts, xlabelvisible=false,xticklabelsize=15,yticklabelsize=15,xlabelsize=18,ylabelsize=18)
# rprof_vs_field(noK, "Water", fig=fig, fpos=(4,1), cmap_reverse=true, logy=true, colorrange=crange, cbar=false, interpolate=itp, tstart=ts, xlabelvisible=false,xticklabelsize=15,yticklabelsize=15,xlabelsize=18,ylabelsize=18)
# rprof_vs_field(y35, "Water", fig=fig, fpos=(5,1), cmap_reverse=true, logy=true, colorrange=crange, cbar=false, interpolate=itp, tstart=ts,xticklabelsize=15,yticklabelsize=15,xlabelsize=18,ylabelsize=18)

# Fig 3
# ts = 3.0
# mantle_water(y25, fluxes=true, Tg_Myr=true, tstart=ts, fig=fig, fpos=(1,1), yrange=(-7e9, 7e9), disp=false)
# mantle_water(y30, fluxes=true, Tg_Myr=true, tstart=ts, fig=fig, fpos=(2,1), yrange=(-7e9, 7e9), disp=false, leg=false)
# mantle_water(y35, fluxes=true, Tg_Myr=true, tstart=ts, fig=fig, fpos=(3,1), yrange=(-7e9, 7e9), leg=false, savein="fluxes")

# # Fig 4
# time_vs_ratio([y20, y25, y30, y35], "EruptedH2O", "erupta", tstart=0.5, mov_avg=true, line=false, linewidth=4.0, mavg_contrast=0.2, fig=fig, fpos=(1,1), logy=true, ylabelsize=13, xlabelvisible=false)
# time_vs_ratio([y20, y25, y30, y35], "∂EruptedH2O", "∂erupta", tstart=0.5, mov_avg=true, line=false, linewidth=4.0, mavg_contrast=0.2, fig=fig, fpos=(2,1), logy=true, yrange=(1e-6, 1e-2), leg=false, ylabelsize=13, xlabelvisible=false)
# time_vs_field(y20, "SurfOceanMass3D", logy=false, tstart=2.5, fig=fig, fpos=(3,1), ylabelsize=15)

# Fig 5
# ts = 0.1
# sim = noK
# sim2 = y25
# time_vs_ratio([noK, y25], "EruptedH2O", "erupta", tstart=0.5, mov_avg=true, line=false, linewidth=4.0, mavg_contrast=0.2, fig=fig, fpos=(1,1), logy=true)
# time_vs_ratio([sim, sim2], "OutgassedH2O", "IngassedH2O", tstart=ts, mov_avg=true, line=false, linewidth=4.0, mavg_contrast=0.2, fig=fig, fpos=(1,2), logy=false, legpos=:rb, framevisible=false, legflat=true, disp=false, xlabelvisible=false)
# time_vs_field([sim, sim2], "Water_tz", tstart=ts, mov_avg=true, line=false, linewidth=4.0, mavg_contrast=0.15, fig=fig, fpos=(2,3), logy=true, leg=false, disp=false, xlabelvisible=false)
# time_vs_field([sim, sim2], "Wsol_tz", tstart=ts, mov_avg=true, line=false, linewidth=4.0, mavg_contrast=0.15, fig=fig, fpos=(2,2), logy=true, leg=false, disp=false, xlabelvisible=false)
# time_vs_ratio([sim, sim2], "freeH2O", "Water", tstart=ts, mov_avg=true, line=false, linewidth=4.0, mavg_contrast=0.2, fig=fig, fpos=(1,3), logy=false, leg=false, disp=false, xlabelvisible=false)
# rprof_vs_field(y25, "Water", fig=fig, fpos=(3,2:3), cmap_reverse=true, logy=true, colorrange=crange, cbar=false, interpolate=itp, tstart=ts, xlabelvisible=false,xticklabelsize=15,yticklabelsize=15,xlabelsize=18,ylabelsize=18)
# rprof_vs_field(noK, "Water", fig=fig, fpos=(4,2:3), cmap_reverse=true, logy=true, colorrange=crange, cbar=false, interpolate=itp, tstart=ts, xlabelvisible=false,xticklabelsize=15,yticklabelsize=15,xlabelsize=18,ylabelsize=18)
# mantle_water(sim, tstart=ts, fig=fig, fpos=(3,2:3), legflat=true, fluxes=false, Tg_Myr=true, framevisible=false, yrange=(0.0, 1.25), disp=false, mavg=true, legfontsize=10, xlabelvisible=false)
# mantle_water(sim2, tstart=ts, fig=fig, fpos=(4,2:3), leg=false, fluxes=false, Tg_Myr=true, yrange=(0.0, 1.25), disp=false, mavg=true)
# time_vs_field([sim, sim2], "fmeltmean_crust", tstart=ts, mov_avg=true, line=false, linewidth=4.0, mavg_contrast=0.2, fig=fig, fpos=(1,1), logy=true, leg=false, disp=false, xlabelvisible=false)
# time_vs_field([sim, sim2], "fmeltmean_um", tstart=ts, mov_avg=true, line=false, linewidth=4.0, mavg_contrast=0.2, fig=fig, fpos=(2,1), logy=true, leg=false, disp=false, xlabelvisible=false)
# time_vs_field([sim, sim2], "fmeltmean_tz", tstart=ts, mov_avg=true, line=false, linewidth=4.0, mavg_contrast=0.2, fig=fig, fpos=(3,1), logy=true, leg=false, disp=false, xlabelvisible=false)
# time_vs_field([sim, sim2], "fmeltmean_lm", tstart=ts, mov_avg=true, line=false, linewidth=4.0, mavg_contrast=0.2, fig=fig, fpos=(4,1), logy=true, leg=false, savein="noK")


# Fig 6
# time_vs_field([logp, logp05, logm05, logm], "Water_um", tstart=ts, mov_avg=true, line=false, linewidth=4.0, mavg_contrast=0.1, fig=fig, fpos=(1,1), logy=true, leg=false, disp=false, xlabelvisible=false)
# time_vs_field([logp, logp05, logm05, logm], "Water_tz", tstart=ts, mov_avg=true, line=false, linewidth=4.0, mavg_contrast=0.1, fig=fig, fpos=(2,1), logy=true, leg=false, disp=false, xlabelvisible=false)
# time_vs_field([logp, logp05, logm05, logm], "Water_lm", tstart=ts, mov_avg=true, line=false, linewidth=4.0, mavg_contrast=0.1, fig=fig, fpos=(3,1), logy=true, leg=false)
# time_vs_field([logp, logp05, logm05, logm], "Tpotl", tstart=ts, mov_avg=true, line=false, linewidth=4.0, mavg_contrast=0.2, fig=fig, fpos=(1,2), legfontsize=15, logy=false, legpos=:rt, framevisible=false, legflat=true, disp=false, xlabelvisible=false)
# time_vs_field([logp, logp05, logm05, logm], "SurfOceanMass3D", tstart=ts, mov_avg=true, line=false, linewidth=4.0, mavg_contrast=0.10, fig=fig, fpos=(2,2), logy=false, leg=false, disp=false, xlabelvisible=false)
# time_vs_ratio([logp, logp05, logm05, logm], "freeH2O", "Water", tstart=ts, mov_avg=true, line=false, linewidth=4.0, mavg_contrast=0.1, fig=fig, fpos=(3,2), logy=false, leg=false, disp=false, xlabelvisible=false, savein="darcy")
# mantle_water(logp, tstart=ts, fig=fig, fpos=(3,2:3), legflat=true, fluxes=false, Tg_Myr=true, framevisible=false, yrange=(0.0, 1.5), disp=false, mavg=true, legfontsize=10, xlabelvisible=false)
# mantle_water(logm, tstart=ts, fig=fig, fpos=(4,2:3), leg=false, fluxes=false, Tg_Myr=true, yrange=(0.0, 1.5), disp=false, mavg=true)

# Fig 7
# mix = [y, umdarcy, fulldarcy]
# ta_field_vs_field(mix, "Mob", "satH2O_um", tstart=0.2, tend=4.5, labels=false, line=false, markersize=20, savein="", clrby="Water_um", cmap=:Blues, # setticks=["γₛ", "UM κ₀", "full κ₀"], 
#                     logx=false, strokewidth=0.5, marker=:utriangle, deg="root", dfit=true, color=:gray, fig=fig, fpos=(1,1), yrange=(35, 65))
# ta_field_vs_field(mix, "Mob", "satH2O_tz", tstart=0.2, tend=4.5, labels=false, line=false, markersize=20, savein="", clrby="Water_tz", cmap=:Blues, # setticks=["γₛ", "UM κ₀", "full κ₀"], 
#                     logx=false, strokewidth=0.5, marker=:utriangle, deg="log", dfit=true, color=:gray, fig=fig, fpos=(2,1), yrange=(12, 28.))
# ta_field_vs_field(mix, "Mob", "satH2O_lm", tstart=0.2, tend=4.5, labels=false, line=false, markersize=20, savein="sat", clrby="Water_lm", cmap=:Blues, # setticks=["γₛ", "UM κ₀", "full κ₀"], 
#                     logx=false, strokewidth=0.5, marker=:utriangle, deg="root", dfit=true, color=:gray, fig=fig, fpos=(3,1), yrange=(20., 60.))

# time_vs_field([y20, y25], "freeH2O_um", ylabelsize=15, yticklabelsize=15, xticklabelsize=15, xlabelvisible=false, logy=true, fig=fig, fpos=(1,1), tstart=0.5, mov_avg=true, mavg_contrast=0.12, linewidth=3.0, legpos=:rb, framevisible=false, legflat=true, legfontsize=25)
# time_vs_field([y20, y25], "boundH2O_um", ylabelsize=20, yticklabelsize=15, xticklabelsize=15, leg=false, logy=true, fig=fig, fpos=(2,1), tstart=0.5, mov_avg=true, mavg_contrast=0.12, linewidth=3.0, framevisible=false, legflat=true, legfontsize=25)
# time_vs_field([y20, y25], "Water_um", ylabelsize=20, yticklabelsize=15, xticklabelsize=15, leg=false, logy=true, fig=fig, fpos=(3,1), tstart=0.5, mov_avg=true, mavg_contrast=0.12, linewidth=3.0, framevisible=false, legflat=true, legfontsize=25)
# time_vs_field([y20, y25], "vzmax_um", ylabelsize=20, yticklabelsize=15, xticklabelsize=15, yrange=(1e-8, 1e-6), leg=false, logy=true, fig=fig, fpos=(4,1), tstart=0.5, mov_avg=true, mavg_contrast=0.12, linewidth=3.0, framevisible=false, legflat=true, legfontsize=25, savein="y20y25")

# ta_field_vs_field([y25, y30, y35, y40, y45, y50, y55], "Mob", "satH2O_um", clrby="Wsol_um", line=false, tstart=2.0, tend=4.5, labels=true, markersize=20, cmap_reverse=false, fig=fig, fpos=(1,1),
#                     ylabelsize=20, yticklabelsize=15, xticklabelsize=15, xlabelsize=20, clabelsize=20, color=:black, logx=true)

# crange=(-3.5, 0.25); itp=true; ts=0.5
# f = "bsmean"
# rprof_vs_field(y35, f, fig=fig, fpos=(1,1), cmap_reverse=false, logy=false, interpolate=itp, tstart=ts, xlabelvisible=false,xticklabelsize=15,yticklabelsize=15,xlabelsize=18,ylabelsize=18)
# rprof_vs_field(y40, f, fig=fig, fpos=(3,1), cmap_reverse=false, logy=false, cbar=false, interpolate=itp, tstart=ts, xlabelvisible=false,xticklabelsize=15,yticklabelsize=15,xlabelsize=18,ylabelsize=18)
# rprof_vs_field(y45, f, fig=fig, fpos=(4,1), cmap_reverse=false, logy=false, cbar=false, interpolate=itp, tstart=ts, xlabelvisible=false,xticklabelsize=15,yticklabelsize=15,xlabelsize=18,ylabelsize=18)
# rprof_vs_field(y50, f, fig=fig, fpos=(5,1), cmap_reverse=false, logy=false, cbar=false, interpolate=itp, tstart=ts,xticklabelsize=15,yticklabelsize=15,xlabelsize=18,ylabelsize=18)

# Fig 8
# mix = [umdarcy, fulldarcy]
# vmix = vcat(umdarcy, fulldarcy)
time_vs_ratio(umdarcy, "EruptedH2O", "erupta", cmap_reverse=true, hmap=true, nmap=3000, tstart=0.3, fig=fig, fpos=(1,1), disp=false)
time_vs_ratio(fulldarcy, "EruptedH2O", "erupta", cmap_reverse=true, hmap=true, nmap=3000, tstart=0.3, fig=fig, fpos=(1,3))
# ta_field_vs_field(mix, "fmeltmean", "influxH2O", deg="fr", line=false, labels=false, marker=:utriangle, strokewidth=1.0, markersize=20, color=:gray)