using StagPP, GLMakie

load_sims = false
# Sims
    if load_sims
        y20 = load_sim("/Volumes/Sinergia/EarthWater/y20e20/+op/EW", "γ20");
        y25 = load_sim("/Volumes/Sinergia/EarthWater/y25e20/+op/EW", "γ25");
        y30 = load_sim("/Volumes/Sinergia/EarthWater/y30e20/+op/EW", "γ30");
        y35 = load_sim("/Volumes/Sinergia/EarthWater/y35e20/+op/EW", "γ35");
        y40 = load_sim("/Volumes/Sinergia/EarthWater/y40e20/+op/EW", "γ40");
        y45 = load_sim("/Volumes/Sinergia/EarthWater/y45e20/+op/EW", "γ45");
        y50 = load_sim("/Volumes/Sinergia/EarthWater/y50e20/+op/EW", "γ50");
        y55 = load_sim("/Volumes/Sinergia/EarthWater/y55e20/+op/EW", "γ55");

        # sets
        y = [y20, y25, y30, y35, y40, y45, y50, y55];
    end
#

# time_vs_field(y, "Water_lm")

a = solve_sH2O_fO2(70, 70, plt=true, cmap=:Blues)
