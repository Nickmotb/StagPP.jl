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

# b = Hirschmann_fO2_to_R(RF=0.001, FMQ=true, P=13.0, T=1800, p=0.0, fsize=(1000, 600));
# a = solve_sH2O_fO2(70, 70, nR=4);
# a2 = solve_sH2O_fO2(70, 70, nR=1, Rv=-0.07, s=false);
a = solve_sH2O_fO2(30, 30, nR=10, Rrange=(0.0005, 0.01), s=false);
plot_sf(a, fsize=(1200,800), FMQ=false, fblck=true)

# out = solve_point(50., 1500., "XH")