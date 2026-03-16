using StagPP, GLMakie

load_sims = false
# Sims
    if load_sims
        y20 = load_sim("/Volumes/Sinergia/EarthWater/y20e19/+op/EW", "γ20");
        y25 = load_sim("/Volumes/Sinergia/EarthWater/y25e19/+op/EW", "γ25");
        y30 = load_sim("/Volumes/Sinergia/EarthWater/y30e19/+op/EW", "γ30");
        y35 = load_sim("/Volumes/Sinergia/EarthWater/y35e19/+op/EW", "γ35");
        y40 = load_sim("/Volumes/Sinergia/EarthWater/y40e19/+op/EW", "γ40");
        y45 = load_sim("/Volumes/Sinergia/EarthWater/y45e19/+op/EW", "γ45");
        y50 = load_sim("/Volumes/Sinergia/EarthWater/y50e19/+op/EW", "γ50");
        y55 = load_sim("/Volumes/Sinergia/EarthWater/y55e19/+op/EW", "γ55");

        # # sets
        y = [y20, y25, y30, y35, y40, y45, y50, y55];
    end
#

# sim = load_local("office")
# rprof_vs_field(y45, "Water", logscale=true, cmap=:Blues, cmap_reverse=false, colorrange=(-3, 1))
# mantle_water_at_t(y35, 4.5)

Δ = Hirschmann_fO2_to_R(RF=0.01, FMQ=true, P=3.5, T=1200., p=0.0);
# a = solve_sH2O_fO2(30, 30, s=false, Rv=0.05); plot_sf(a, FMQ=true)
# 