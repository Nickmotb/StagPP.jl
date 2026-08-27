# using StagPP, GLMakie, MAGEMin_C

# P = 4.0
# T = 1500.
# p = 0.0
# ϕ = 0.05
# TOex = 0.005
# TC = 0.05
# H = 1.0

# data = Initialize_MAGEMin("sb24", verbose=false);
# partition_Oₑₓ(P, T, p, ϕ, TOex, TC, plotevo=true, debugging=true, data=data)

using StagPP, GLMakie

load_sims = false
# Sims
    if load_sims
        # ys25
        y25 = load_sim("/Volumes/Sinergia/EarthWater/ys25/+op/EW", "γ25");
        y25_b02_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys25_blh02_pf1/+op/EW", "γ25_b02_p1");
        y25_b04_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys25_blh04_pf1/+op/EW", "γ25_b04_p1");
        y25_b06_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys25_blh06_pf1/+op/EW", "γ25_b06_p1");
        y25_b08_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys25_blh08_pf1/+op/EW", "γ25_b08_p1");
        y25_b1_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys25_blh1_pf1/+op/EW", "γ25_b1_p1");
        y25_b1_p1_Int = load_sim("/Volumes/Sinergia/EarthWater/ys25_blh1_pf1_Int/+op/EW", "γ25_b1_p1_Int");
        y25_p1 = [y25, y25_b02_p1, y25_b04_p1, y25_b06_p1, y25_b08_p1, y25_b1_p1]
        # ys30
        y30 = load_sim("/Volumes/Sinergia/EarthWater/ys30/+op/EW", "γ30");
        y30_b02_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys30_blh02_pf1/+op/EW", "γ30_b02_p1");
        y30_b04_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys30_blh04_pf1/+op/EW", "γ30_b04_p1");
        y30_b06_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys30_blh06_pf1/+op/EW", "γ30_b06_p1");
        y30_b08_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys30_blh08_pf1/+op/EW", "γ30_b08_p1");
        y30_b1_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys30_blh1_pf1/+op/EW", "γ30_b1_p1");
        y30_b1_p1_Int = load_sim("/Volumes/Sinergia/EarthWater/ys30_blh1_pf1_Int/+op/EW", "γ30_b1_p1_Int");
        y30_p1 = [y30, y30_b02_p1, y30_b04_p1, y30_b06_p1, y30_b08_p1, y30_b1_p1]
        # ys35
        y35 = load_sim("/Volumes/Sinergia/EarthWater/ys35/+op/EW", "γ35");
        y35_b02_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys35_blh02_pf1/+op/EW", "γ35_b02_p1");
        y35_b04_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys35_blh04_pf1/+op/EW", "γ35_b04_p1");
        y35_b06_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys35_blh06_pf1/+op/EW", "γ35_b06_p1");
        y35_b08_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys35_blh08_pf1/+op/EW", "γ35_b08_p1");
        y35_b1_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys35_blh1_pf1/+op/EW", "γ35_b1_p1");
        y35_b1_p1_Int = load_sim("/Volumes/Sinergia/EarthWater/ys35_blh1_pf1_Int/+op/EW", "γ35_b1_p1_Int");
        y35_p1 = [y35, y35_b02_p1, y35_b04_p1, y35_b06_p1, y35_b08_p1, y35_b1_p1]
        # ys40
        y40 = load_sim("/Volumes/Sinergia/EarthWater/ys40/+op/EW", "γ40");
        y40_b02_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys40_blh02_pf1/+op/EW", "γ40_b02_p1");
        y40_b04_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys40_blh04_pf1/+op/EW", "γ40_b04_p1");
        y40_b06_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys40_blh06_pf1/+op/EW", "γ40_b06_p1");
        y40_b08_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys40_blh08_pf1/+op/EW", "γ40_b08_p1");
        y40_b1_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys40_blh1_pf1/+op/EW", "γ40_b1_p1");
        y40_b1_p1_Int = load_sim("/Volumes/Sinergia/EarthWater/ys40_blh1_pf1_Int/+op/EW", "γ40_b1_p1_Int");
        y40_p1 = [y40, y40_b02_p1, y40_b04_p1, y40_b06_p1, y40_b08_p1, y40_b1_p1]
        # ys45
        y45 = load_sim("/Volumes/Sinergia/EarthWater/ys45/+op/EW", "γ45");
        y45_b02_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys45_blh02_pf1/+op/EW", "γ45_b02_p1");
        y45_b04_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys45_blh04_pf1/+op/EW", "γ45_b04_p1");
        y45_b06_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys45_blh06_pf1/+op/EW", "γ45_b06_p1");
        y45_b08_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys45_blh08_pf1/+op/EW", "γ45_b08_p1");
        y45_b1_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys45_blh1_pf1/+op/EW", "γ45_b1_p1");
        y45_b1_p1_Int = load_sim("/Volumes/Sinergia/EarthWater/ys45_blh1_pf1_Int/+op/EW", "γ45_b1_p1_Int");
        y45_p1 = [y45, y45_b02_p1, y45_b04_p1, y45_b06_p1, y45_b08_p1, y45_b1_p1]
        # ys50
        y50 = load_sim("/Volumes/Sinergia/EarthWater/ys50/+op/EW", "γ50");
        y50_b02_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys50_blh02_pf1/+op/EW", "γ50_b02_p1");
        y50_b04_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys50_blh04_pf1/+op/EW", "γ50_b04_p1");
        y50_b06_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys50_blh06_pf1/+op/EW", "γ50_b06_p1");
        y50_b08_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys50_blh08_pf1/+op/EW", "γ50_b08_p1");
        y50_b1_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys50_blh1_pf1/+op/EW", "γ50_b1_p1");
        y50_b1_p1_Int = load_sim("/Volumes/Sinergia/EarthWater/ys50_blh1_pf1_Int/+op/EW", "γ50_b1_p1_Int");
        y50_p1 = [y50, y50_b02_p1, y50_b04_p1, y50_b06_p1, y50_b08_p1, y50_b1_p1]
        # ys55
        y55 = load_sim("/Volumes/Sinergia/EarthWater/ys55/+op/EW", "γ55");
        y55_b02_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys55_blh02_pf1/+op/EW", "γ55_b02_p1");
        y55_b04_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys55_blh04_pf1/+op/EW", "γ55_b04_p1");
        y55_b06_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys55_blh06_pf1/+op/EW", "γ55_b06_p1");
        y55_b08_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys55_blh08_pf1/+op/EW", "γ55_b08_p1");
        y55_b1_p1 = load_sim("/Volumes/Sinergia/EarthWater/ys55_blh1_pf1/+op/EW", "γ55_b1_p1");
        y55_b1_p1_Int = load_sim("/Volumes/Sinergia/EarthWater/ys55_blh1_pf1_Int/+op/EW", "γ55_b1_p1_Int");
        y55_p1 = [y55, y55_b02_p1, y55_b04_p1, y55_b06_p1, y55_b08_p1, y55_b1_p1]

        block_p1 = [y25_b02_p1  y30_b02_p1  y35_b02_p1  y40_b02_p1  y45_b02_p1  y50_b02_p1  y55_b02_p1;
                    y25_b04_p1  y30_b04_p1  y35_b04_p1  y40_b04_p1  y45_b04_p1  y50_b04_p1  y55_b04_p1;
                    y25_b06_p1  y30_b06_p1  y35_b06_p1  y40_b06_p1  y45_b06_p1  y50_b06_p1  y55_b06_p1;
                    y25_b08_p1  y30_b08_p1  y35_b08_p1  y40_b08_p1  y45_b08_p1  y50_b08_p1  y55_b08_p1;
                    y25_b1_p1   y30_b1_p1   y35_b1_p1   y40_b1_p1   y45_b1_p1   y50_b1_p1   y55_b1_p1]
        
        block_int1 = [y25 y30 y35 y40 y45 y50 y55;
                      y25_b1_p1 y30_b1_p1 y35_b1_p1 y40_b1_p1 y45_b1_p1 y50_b1_p1 y55_b1_p1;
                      y25_b1_p1_Int y30_b1_p1_Int y35_b1_p1_Int y40_b1_p1_Int y45_b1_p1_Int y50_b1_p1_Int y55_b1_p1_Int]
    
        block_nw = [y25 y30 y35 y40 y45 y50 y55;
                    y25_b1_p1 y30_b1_p1 y35_b1_p1 y40_b1_p1 y45_b1_p1 y50_b1_p1 y55_b1_p1]
    end

#     a = solve_sH2O_fO2(300, 300, fO2=true, exprt=false, densities=false)
#     plot_sf(a, savein="maps")

    choose_figure = 2

    xticks = [25., 30., 35., 40., 45., 50., 55.]
    yticks = [0.2, 0.4, 0.6, 0.8, 1.0]

    # time_vs_field(y30_b1_p1, "Mob")
    
    # General inputs
    tstart=0.0; tend=1.0; save=false

    if choose_figure==1
        # ===== H2O Regimes =====
        fig = Figure(size=(900, 1600))
        sect="_lm"; hmalpha=0.57; txtsz=11; interp=true
        clrrng = [(-1, -1), (-1, -1), (-1, -1)]
        # clrrng = sect=="_um" ? [(-2.33, -0.985), (-1.41, -0.548), (-3.94, -3.44)] : 
        #             sect=="_tz" ? [(-1.5, -0.98), (-0.679, -0.26), (-3.84, -3.02)] :
        #               sect=="_lm" ? [(-2.1, -1.5),  (-1.32, -0.8), (-3.86, -1.62)] : 
        #                 [(-1.925, -1.284), (-1.16, -0.8), (-3.628, -1.83)]
        regime(block_p1, "boundH2O"*sect, xticks, yticks, interpolate=interp, cmap_reverse=false, logy=true, fig=fig, fpos=(1,1), tstart=tstart, tend=tend, hmalpha=hmalpha, strokewidth=1.0,
                xlabelsize=23, ylabelsize=23, xticklabelsize=15, yticklabelsize=15, colorrange=clrrng[1], cmap=:Blues, txtsz=txtsz, disp=false)

        regime(block_p1, "Wsol"*sect, xticks, yticks, interpolate=interp, cmap_reverse=false, logy=true, fig=fig, fpos=(3,1), tstart=tstart, tend=tend, hmalpha=hmalpha, strokewidth=1.0,
                xlabelsize=23, ylabelsize=23, xticklabelsize=15, yticklabelsize=15, colorrange=clrrng[2], cmap=:Greys, txtsz=txtsz, disp=false)

        regime(block_p1, "freeH2O"*sect, xticks, yticks, interpolate=interp, cmap_reverse=true, logy=true, fig=fig, fpos=(5,1), tstart=tstart, tend=tend, hmalpha=hmalpha, strokewidth=1.0,
                xlabelsize=23, ylabelsize=23, xticklabelsize=15, yticklabelsize=15, colorrange=clrrng[3], cmap=:vik, txtsz=txtsz, disp=false)

        save && GLMakie.save("./scripts/Paper1/imgs/H2Oregimes.png", fig)

    elseif choose_figure == 2
        # ===== Mobility regimes =====
        fig = Figure(size=(1300, 900))
        cmap = :tofino; cmap_reverse=false; colorrange=(-1, -1); interp=true; hmalpha=0.57
        regime(block_p1, "Mob", xticks, yticks, interpolate=interp, cmap_reverse=cmap_reverse, logy=false, fig=fig, fpos=(1,1), tstart=tstart, tend=tend, hmalpha=hmalpha, strokewidth=1.0,
                xlabelsize=23, ylabelsize=23, xticklabelsize=15, yticklabelsize=15, colorrange=colorrange, cmap=cmap, tectonic=true, disp=false,
                ylab=L"Hydrated\;boundary\;layer\;saturation", xlab=L"Surface\;yield\;stress\;\mathrm{[MPa]}")
        regime(block_int1, "Mob", xticks, 1:3, syticks=["No Water", "Water", "Water + Intrusion"], interpolate=interp, cmap_reverse=cmap_reverse, logy=false, fig=fig, fpos=(3,1), tstart=tstart, tend=tend, hmalpha=hmalpha, strokewidth=1.0,
                xlabelsize=23, ylabelsize=23, xticklabelsize=15, yticklabelsize=15, colorrange=colorrange, cmap=cmap, tectonic=true, disp=false,
                ylab="", xlab=L"Surface\;yield\;stress\;\mathrm{[MPa]}")

        GLMakie.save("./scripts/Paper1/imgs/Mobregimes.png", fig)

    elseif choose_figure == 3
        # ===== Time-averaged correlation =====
        fig = Figure(size=(900, 900))
        tstart = 0.0
        ta_field_vs_field([block_p1[1,:], block_p1[2,:], block_p1[3,:], block_p1[4,:], block_p1[5,:]], "Mob", "Water_um", logy=true, logx=true, labels=false, color=:grey60, cmap=:vik, deg="exp", dfit=true, 
                            markersize=15., marker=:rect, strokewidth=1.1, tstart=tstart, tend=4.0, setfit=true, cmap_reverse=true, fitalpha=0.7, xticklabelsize=15, yticklabelsize=15,
                            line=false, fig=fig, fpos=(1,1), clrby="Set", setticks=["20% sat", "40% sat", "60% sat", "80% sat", "100% sat"], xlabelsize=0., clab=false, ylabelsize=20)
        ta_field_vs_field([block_p1[1,:], block_p1[2,:], block_p1[3,:], block_p1[4,:], block_p1[5,:]], "Mob", "Water_tz", logy=true, logx=true, labels=false, color=:grey60, cmap=:vik, deg="exp", dfit=true, 
                            markersize=15., marker=:rect, strokewidth=1.1, tstart=tstart, tend=4.0, setfit=true, cmap_reverse=true, fitalpha=0.7, xticklabelsize=15, yticklabelsize=15,
                            line=false, fig=fig, fpos=(2,1), clrby="Set", xlabelsize=0., cbar=false, ylabelsize=20)
        ta_field_vs_field([block_p1[1,:], block_p1[2,:], block_p1[3,:], block_p1[4,:], block_p1[5,:]], "Mob", "Water_lm", logy=true, logx=true, labels=false, color=:grey60, cmap=:vik, deg="exp", dfit=true, 
                            markersize=15., marker=:rect, strokewidth=1.1, tstart=tstart, tend=4.0, setfit=true, cmap_reverse=true, fitalpha=0.7, xticklabelsize=15, yticklabelsize=15,
                            line=false, fig=fig, fpos=(3,1), clrby="Set", cbar=false, ylabelsize=20)
        save && GLMakie.save("./scripts/Paper1/imgs/tacorr.png", fig)

    elseif choose_figure == 4
        # ===== Radial Profiles ======
        fig = Figure(size=(1600, 950))
        interp = true
        # bH2O
        rprof_vs_field(block_p1[end-1, 1], "boundH2O",  cmap=:Blues,    cmap_reverse=false, fig=fig, fpos=(1, 1), logy=true, colorrange=(-2, 0.5),  xlabelsize=15,  interpolate=interp, disp=false)
        rprof_vs_field(block_p1[end-1, 2], "boundH2O",  cmap=:Blues,    cmap_reverse=false, fig=fig, fpos=(3, 1), logy=true, colorrange=(-2, 0.5),  ylabelsize=0.0, interpolate=interp, xlabelvisible=false, plab=false, cbar=false, disp=false)
        rprof_vs_field(block_p1[end-1, 3], "boundH2O",  cmap=:Blues,    cmap_reverse=false, fig=fig, fpos=(4, 1), logy=true, colorrange=(-2, 0.5),  ylabelsize=0.0, interpolate=interp, xlabelvisible=false, plab=false, cbar=false, disp=false)
        rprof_vs_field(block_p1[end-1, 6], "boundH2O",  cmap=:Blues,    cmap_reverse=false, fig=fig, fpos=(5, 1), logy=true, colorrange=(-2, 0.5),  ylabelsize=0.0, interpolate=interp, xlabelvisible=false, plab=false, cbar=false, disp=false)
        # sH2O
        rprof_vs_field(block_p1[end-1, 1], "Wsol",      cmap=:Greys,    cmap_reverse=false, fig=fig, fpos=(1, 2), logy=true, colorrange=(-2, 0.5),  xlabelsize=0.0, interpolate=interp, plab=false, ylabelsize=0.0, leg=false, disp=false)
        rprof_vs_field(block_p1[end-1, 2], "Wsol",      cmap=:Greys,    cmap_reverse=false, fig=fig, fpos=(3, 2), logy=true, colorrange=(-2, 0.5),  ylabelsize=0.0, interpolate=interp, xlabelvisible=false, plab=false, cbar=false, leg=false, disp=false)
        rprof_vs_field(block_p1[end-1, 3], "Wsol",      cmap=:Greys,    cmap_reverse=false, fig=fig, fpos=(4, 2), logy=true, colorrange=(-2, 0.5),  ylabelsize=0.0, interpolate=interp, xlabelvisible=false, plab=false, cbar=false, leg=false, disp=false)
        rprof_vs_field(block_p1[end-1, 6], "Wsol",      cmap=:Greys,    cmap_reverse=false, fig=fig, fpos=(5, 2), logy=true, colorrange=(-2, 0.5),  ylabelsize=0.0, interpolate=interp, xlabelvisible=false, plab=false, cbar=false, leg=false, disp=false)
        # Melt
        rprof_vs_field(block_p1[end-1, 1], "fmeltmean", cmap=:thermal,  cmap_reverse=false, fig=fig, fpos=(1, 3), logy=true, colorrange=(-4, -1.5), ylabelsize=0.0, interpolate=interp, xlabelvisible=false, plab=false, leg=false, disp=false)
        rprof_vs_field(block_p1[end-1, 2], "fmeltmean", cmap=:thermal,  cmap_reverse=false, fig=fig, fpos=(3, 3), logy=true, colorrange=(-4, -1.5), ylabelsize=0.0, interpolate=interp, xlabelvisible=false, cbar=false, leg=false, disp=false)
        rprof_vs_field(block_p1[end-1, 3], "fmeltmean", cmap=:thermal,  cmap_reverse=false, fig=fig, fpos=(4, 3), logy=true, colorrange=(-4, -1.5), ylabelsize=0.0, interpolate=interp, xlabelvisible=false, cbar=false, leg=false, disp=false)
        rprof_vs_field(block_p1[end-1, 6], "fmeltmean", cmap=:thermal,  cmap_reverse=false, fig=fig, fpos=(5, 3), logy=true, colorrange=(-4, -1.5), ylabelsize=0.0, interpolate=interp, xlabelvisible=false, cbar=false, leg=false, disp=false)
        # Density
        # rprof_vs_field(block_p1[end-1, 1], "etalog", cmap=:speed,  cmap_reverse=false, fig=fig, fpos=(1, 3), logy=true, colorrange=(18, 25), ylabelsize=0.0, interpolate=interp, xlabelvisible=false, plab=false, leg=false, disp=false)
        # rprof_vs_field(block_p1[end-1, 2], "etalog", cmap=:speed,  cmap_reverse=false, fig=fig, fpos=(3, 3), logy=true, colorrange=(18, 25), ylabelsize=0.0, interpolate=interp, xlabelvisible=false, cbar=false, leg=false, disp=false)
        # rprof_vs_field(block_p1[end-1, 3], "etalog", cmap=:speed,  cmap_reverse=false, fig=fig, fpos=(4, 3), logy=true, colorrange=(18, 25), ylabelsize=0.0, interpolate=interp, xlabelvisible=false, cbar=false, leg=false, disp=false)
        # rprof_vs_field(block_p1[end-1, 6], "etalog", cmap=:speed,  cmap_reverse=false, fig=fig, fpos=(5, 3), logy=true, colorrange=(18, 25), ylabelsize=0.0, interpolate=interp, xlabelvisible=false, cbar=false, leg=false, disp=false)
        # fH2O
        rprof_vs_field(block_p1[end-1, 1], "freeH2O",   cmap=:vik100,   cmap_reverse=true,  fig=fig, fpos=(1, 4), logy=true, colorrange=(-4, -1),   xlabelsize=0.0, interpolate=interp, plab=false, ylabelsize=0.0, leg=false, disp=false)
        rprof_vs_field(block_p1[end-1, 2], "freeH2O",   cmap=:vik100,   cmap_reverse=true,  fig=fig, fpos=(3, 4), logy=true, colorrange=(-4, -1),   ylabelsize=0.0, interpolate=interp, xlabelvisible=false, plab=false, cbar=false, leg=false, disp=false)
        rprof_vs_field(block_p1[end-1, 3], "freeH2O",   cmap=:vik100,   cmap_reverse=true,  fig=fig, fpos=(4, 4), logy=true, colorrange=(-4, -1),   ylabelsize=0.0, interpolate=interp, xlabelvisible=false, plab=false, cbar=false, leg=false)
        rprof_vs_field(block_p1[end-1, 6], "freeH2O",   cmap=:vik100,   cmap_reverse=true,  fig=fig, fpos=(5, 4), logy=true, colorrange=(-4, -1),   ylabelsize=0.0, interpolate=interp, xlabelvisible=false, plab=false, cbar=false, leg=false)
        save && GLMakie.save("./scripts/Paper1/imgs/radial.png", fig)

    elseif choose_figure == 5
        # ===== H2O concentration of magmas ======
        fig = Figure(size=(1600, 900))
        time_vs_ratio(block_p1[1,:], "∂SaturationOutgasH2O", "∂erupta", line=false, mov_avg=true, logy=true, fig=fig, fpos=(1,1), ylabelsize=0, xlabelsize=0, yrange=[1e-12, 1e-2], yticklabelsize=15, xticklabelsize=15)
        time_vs_ratio(block_p1[2,:], "∂SaturationOutgasH2O", "∂erupta", line=false, mov_avg=true, logy=true, fig=fig, fpos=(1,2), ylabelsize=0, xlabelsize=0, yrange=[1e-12, 1e-2], yticklabelsize=15, xticklabelsize=15)
        time_vs_ratio(block_p1[3,:], "∂SaturationOutgasH2O", "∂erupta", line=false, mov_avg=true, logy=true, fig=fig, fpos=(2,1), ylabelsize=0, xlabelsize=0, yrange=[1e-12, 1e-2], yticklabelsize=15, xticklabelsize=15)
        time_vs_ratio(block_p1[4,:], "∂SaturationOutgasH2O", "∂erupta", line=false, mov_avg=true, logy=true, fig=fig, fpos=(2,2), ylabelsize=0, xlabelsize=0, yrange=[1e-12, 1e-2], yticklabelsize=15, xticklabelsize=15)
        time_vs_ratio(block_p1[5,:], "∂SaturationOutgasH2O", "∂erupta", line=false, mov_avg=true, logy=true, fig=fig, fpos=(3,1), yrange=[1e-12, 1e-2], yticklabelsize=15, xticklabelsize=15)
        save && GLMakie.save("./scripts/Paper1/imgs/H2Oconc.png", fig)

    elseif choose_figure == 6
        # ===== Surface ocean ======
        fig = Figure(size=(1600, 1100))
        mantle_water_at_t(block_p1[1,1], 2.5, fig=fig, fpos=(1,1))
        mantle_water_at_t(block_p1[3,1], 2.5, fig=fig, fpos=(2,1))
        mantle_water_at_t(block_p1[5,1], 2.5, fig=fig, fpos=(3,1))

        mantle_water_at_t(block_p1[1,4], 2.5, fig=fig, fpos=(1,2))
        mantle_water_at_t(block_p1[3,4], 2.5, fig=fig, fpos=(2,2))
        mantle_water_at_t(block_p1[5,4], 2.5, fig=fig, fpos=(3,2))

        mantle_water_at_t(block_p1[1,7], 2.5, fig=fig, fpos=(1,3))
        mantle_water_at_t(block_p1[3,7], 2.5, fig=fig, fpos=(2,3))
        mantle_water_at_t(block_p1[5,7], 2.5, fig=fig, fpos=(3,3))
        # time_vs_field(block_p1[1,1:3], "SurfOceanMass3D", tstart=0.5, tend=4., fig=fig, fpos=(1,1), legflat=false, legpos=:lb, framevisible=false, legfontsize=10, xlabelsize=0., ylabelsize=0., xticklabelsize=15, yticklabelsize=15)
        # time_vs_field(block_p1[2,1:3], "SurfOceanMass3D", tstart=0.5, tend=4., fig=fig, fpos=(1,3), legflat=false, legpos=:lb, framevisible=false, legfontsize=10, xlabelsize=0., ylabelsize=0., xticklabelsize=15, yticklabelsize=15)
        # time_vs_field(block_p1[3,1:3], "SurfOceanMass3D", tstart=0.5, tend=4., fig=fig, fpos=(2,1), legflat=false, legpos=:lb, framevisible=false, legfontsize=10, xticklabelsize=15, yticklabelsize=15)
        # time_vs_field(block_p1[4,1:3], "SurfOceanMass3D", tstart=0.5, tend=4., fig=fig, fpos=(2,3), legflat=false, legpos=:lb, framevisible=false, legfontsize=10, xlabelsize=0., ylabelsize=0., xticklabelsize=15, yticklabelsize=15)

        # time_vs_field(block_p2[1,1:3], "SurfOceanMass3D", tstart=0.5, tend=4., fig=fig, fpos=(1,2), legflat=false, legpos=:lb, framevisible=false, legfontsize=10, xlabelsize=0., ylabelsize=0., xticklabelsize=15, yticklabelsize=15)
        # time_vs_field(block_p2[2,1:3], "SurfOceanMass3D", tstart=0.5, tend=4., fig=fig, fpos=(1,4), legflat=false, legpos=:lb, framevisible=false, legfontsize=10, xlabelsize=0., ylabelsize=0., xticklabelsize=15, yticklabelsize=15)
        # time_vs_field(block_p2[3,1:3], "SurfOceanMass3D", tstart=0.5, tend=4., fig=fig, fpos=(2,2), legflat=false, legpos=:lb, framevisible=false, legfontsize=10, xticklabelsize=15, yticklabelsize=15)
        # time_vs_field(block_p2[4,1:3], "SurfOceanMass3D", tstart=0.5, tend=4., fig=fig, fpos=(2,4), legflat=false, legpos=:lb, framevisible=false, legfontsize=10, xlabelsize=0., ylabelsize=0., xticklabelsize=15, yticklabelsize=15)
        
    
    end
    display(fig)