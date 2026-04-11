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

        logm05 = load_sim("/Volumes/Sinergia/EarthWater/darcy-0.5logm/+op/EW", "-0.5log");
        logp05 = load_sim("/Volumes/Sinergia/EarthWater/darcy-0.5logp/+op/EW", "+0.5log");
        logm = load_sim("/Volumes/Sinergia/EarthWater/darcy-logm/+op/EW", "-log");
        logp = load_sim("/Volumes/Sinergia/EarthWater/darcy-logp/+op/EW", "+log");
        d1_16 = load_sim("/Volumes/Sinergia/EarthWater/darcy1-16/+op/EW", "d1e-16");
        d3_16 = load_sim("/Volumes/Sinergia/EarthWater/darcy3-16/+op/EW", "d3e-16");
        d7_16 = load_sim("/Volumes/Sinergia/EarthWater/darcy7-16/+op/EW", "d7e-16");
        d9_16 = load_sim("/Volumes/Sinergia/EarthWater/darcy9-16/+op/EW", "d9e-16");
        d5_17 = load_sim("/Volumes/Sinergia/EarthWater/darcy5-17/+op/EW", "d5e-17");
        d7_17 = load_sim("/Volumes/Sinergia/EarthWater/darcy7-17/+op/EW", "d7e-17");
        d9_17 = load_sim("/Volumes/Sinergia/EarthWater/darcy9-17/+op/EW", "d9e-17");
        d3_16 = load_sim("/Volumes/Sinergia/EarthWater/darcy3-16/+op/EW", "d3e-16");
        d2_15 = load_sim("/Volumes/Sinergia/EarthWater/darcy2-15/+op/EW", "d2e-15");
        d4_15 = load_sim("/Volumes/Sinergia/EarthWater/darcy4-15/+op/EW", "d4e-15");
        d6_15 = load_sim("/Volumes/Sinergia/EarthWater/darcy6-15/+op/EW", "d6e-15");

        # Ingas depth
        id4 = load_sim("/Volumes/Sinergia/EarthWater/id4/+op/EW", "id4");
        id6 = load_sim("/Volumes/Sinergia/EarthWater/id6/+op/EW", "id6");
        id8 = load_sim("/Volumes/Sinergia/EarthWater/id8/+op/EW", "id8");
        
        # No RK
        # noR = load_sim("/Volumes/Sinergia/EarthWater/noR/+op/EW", "noR");
        noK = load_sim("/Volumes/Sinergia/EarthWater/noK/+op/EW", "noK");
        # noKR = load_sim("/Volumes/Sinergia/EarthWater/noKR/+op/EW", "noKR");

        # Total Water Budget
        tot1 = load_sim("/Volumes/Sinergia/EarthWater/tot1/+op/EW", "1 OM");
        tot2 = load_sim("/Volumes/Sinergia/EarthWater/tot2/+op/EW", "2 OM");
        tot3 = load_sim("/Volumes/Sinergia/EarthWater/tot3/+op/EW", "3 OM");
        tot4 = load_sim("/Volumes/Sinergia/EarthWater/tot4/+op/EW", "4 OM");
        tot5 = load_sim("/Volumes/Sinergia/EarthWater/tot5/+op/EW", "5 OM");

        # # sets
        y = [y20, y25, y30, y35, y40, y45, y50, y55];
        fulldarcy = [logp, logp05, logm05, logm]
        umdarcy = [d5_17, d7_17, d9_17, d1_16, d3_16, d7_16, d9_16, d2_15, d4_15, d6_15]
        # ynw = [y20nw, y25nw, y30nw, y35nw, y40nw, y45nw, y50nw, y55nw];
        # yt = vcat(y, ynw)
        id = [id4, id6, id8, y20];
        tot = [tot1, tot2, tot3, tot4, tot5, y20];
        # no = [y20, noR, noK, noKR]
    end

# test = load_sim("C:\\Users\\Nickm\\.julia\\dev\\StagPP\\test_data1\\EW", "test");
# test2 = load_sim("C:\\Users\\Nickm\\.julia\\dev\\StagPP\\test_data2\\EW", "test");
# t = [test, test2]

a = solve_sH2O_fO2(30, 30, nR=3)
plot_sf(a)