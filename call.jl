using StagPP, GLMakie, MAGEMin_C

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

        # No Water
        y20nw = load_sim("/Volumes/Sinergia/EarthWater/y20e19nw/+op/EW", "γ20nw");
        y25nw = load_sim("/Volumes/Sinergia/EarthWater/y25e19nw/+op/EW", "γ25nw");
        y30nw = load_sim("/Volumes/Sinergia/EarthWater/y30e19nw/+op/EW", "γ30nw");
        y35nw = load_sim("/Volumes/Sinergia/EarthWater/y35e19nw/+op/EW", "γ35nw");
        y40nw = load_sim("/Volumes/Sinergia/EarthWater/y40e19nw/+op/EW", "γ40nw");
        y45nw = load_sim("/Volumes/Sinergia/EarthWater/y45e19nw/+op/EW", "γ45nw");
        y50nw = load_sim("/Volumes/Sinergia/EarthWater/y50e19nw/+op/EW", "γ50nw");
        y55nw = load_sim("/Volumes/Sinergia/EarthWater/y55e19nw/+op/EW", "γ55nw");

        # Ingas depth
        id4 = load_sim("/Volumes/Sinergia/EarthWater/id4/+op/EW", "id4");
        id6 = load_sim("/Volumes/Sinergia/EarthWater/id6/+op/EW", "id6");
        id8 = load_sim("/Volumes/Sinergia/EarthWater/id8/+op/EW", "id8");
        
        # No RK
        noR = load_sim("/Volumes/Sinergia/EarthWater/noR/+op/EW", "noR");
        noK = load_sim("/Volumes/Sinergia/EarthWater/noK/+op/EW", "noK");
        noKR = load_sim("/Volumes/Sinergia/EarthWater/noKR/+op/EW", "noKR");

        # Total Water Budget
        tot1 = load_sim("/Volumes/Sinergia/EarthWater/tot1/+op/EW", "1 OM");
        tot2 = load_sim("/Volumes/Sinergia/EarthWater/tot2/+op/EW", "2 OM");
        tot3 = load_sim("/Volumes/Sinergia/EarthWater/tot3/+op/EW", "3 OM");
        tot4 = load_sim("/Volumes/Sinergia/EarthWater/tot4/+op/EW", "4 OM");
        tot5 = load_sim("/Volumes/Sinergia/EarthWater/tot5/+op/EW", "5 OM");

        # # sets
        y = [y20, y25, y30, y35, y40, y45, y50, y55];
        ynw = [y20nw, y25nw, y30nw, y35nw, y40nw, y45nw, y50nw, y55nw];
        yt = vcat(y, ynw)
        id = [id4, id6, id8, y20];
        tot = [tot1, tot2, tot3, tot4, tot5, y20];
        no = [y20, noR, noK, noKR]
    end

P = 2.0
T = 1300.
p = 0.2             
ϕ = 0.01
TOex = 0.01
TC = 0.01
damp = 0.25
niter=100

# Dump oxygen in oxidizing C → CaCO₃?

# data    = Initialize_MAGEMin("sb24", verbose=false);
partition_Oₑₓ(P, T, p, ϕ, TOex, TC, verbose=true, plotevo=true, data=data, damp=damp, niter=niter)
# Finalize_MAGEMin(data)

# nP  = 15; Pr  = LinRange(1.0, 10.0, nP)
# nT  = 15; Tr  = LinRange(900., 2200., nT)
# nϕ  = 15; ϕr  = LinRange(1e-2, 0.2, nϕ)
# nTO = 15; TOr = LinRange(7e-5, 3e-4, nTO)

# mapϕ, mapTO = P_T_ϕ_TOₑₓ_Rspace(Pr, Tr, ϕr, TOr, data, nres=80)