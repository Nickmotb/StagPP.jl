using StagPP, MAGEMin_C, Printf
const mmFeO, mmFe, mmO = 71.85, 55.85, 16.0
# McDonough & Sun (1995) pyrolite / BSE, the composition Fux et al. model
Xox  = ["SiO2","MgO","FeO","CaO","Al2O3","Na2O","Cr2O3","O"]
oxFe = ["SiO2","MgO","Fe" ,"CaO","Al2O3","Na2O","Cr2O3","O"]
BSE  = [45.00, 37.80, 8.05, 3.55, 4.45, 0.36, 0.38, 0.0]
isFe(p) = lowercase(p) in ("feg","fea","fee","iron","fe")
iw(T,Pb)  = StagPP.logfO2_IW(T, Pb)          # sb24's own IW, the reference TOex.jl now uses
iwf(T,Pb) = StagPP.logfO2_IW_frost(T, Pb)    # the retired Frost (1991) curve, shown for contrast

function bulk_at(X, r)
    nFe = X[3]/mmFeO
    x = copy(X); x[3] = nFe*mmFe; x[8] = nFe*(1 + 0.5r)*mmO
    return 1e2 .* x ./ sum(x)
end
metalfrac(o) = sum(o.ph_frac_wt[i] for i in eachindex(o.ph) if isFe(o.ph[i]); init=0.0)

d = Initialize_MAGEMin("sb24", verbose=false)

# Fux et al. (2024): near-adiabatic T, modelled bulk Fe3+/SumFe (S4.4, Fig. S15),
# and modelled mantle dIW (S4.6, Fig. 8 / Conclusion).
cases = [ (5.0, 1420.0, 0.036, +2.6),
          (8.0, 1440.0, 0.052,  0.0),
          (10.0,1460.0, 0.062, -0.5),
          (12.0,1490.0, 0.072, -1.0),
          (14.0,1540.0, 0.084, -1.0),
          (16.0,1590.0, 0.110, -1.0) ]

println("="^108)
println("TEST B — dIW at Fux et al.'s OWN modelled bulk Fe3+/SumFe")
println("="^108)
@printf("%6s %7s %10s %11s %11s %9s | %11s %9s | %11s\n",
        "P GPa","T C","bulk Fe3+","dIW sb24","dIW Fux","err","dIW Frost","err","metal wt%")
esb, efr = Float64[], Float64[]
for (P,TC,r,dfux) in cases
    T = TC + 273.15
    o = single_point_minimization(10P, TC, d, X=bulk_at(BSE,r), Xoxides=oxFe, sys_in="wt")
    dm, df = o.fO2 - iw(T, 1e4P), o.fO2 - iwf(T, 1e4P)
    push!(esb, abs(dm-dfux)); push!(efr, abs(df-dfux))
    @printf("%6.1f %7.0f %10.3f %+11.3f %+11.3f %+9.3f | %+11.3f %+9.3f | %11.5f\n",
            P, TC, r, dm, dfux, dm-dfux, df, df-dfux, 1e2metalfrac(o))
end
@printf("\nmean |error| against Fux:  sb24 IW %.3f   Frost IW %.3f  (log units)\n",
        sum(esb)/length(esb), sum(efr)/length(efr))

println()
println("="^108)
println("TEST C — ALLOY ABUNDANCE. Fux's Fe-Ni-S-C(-O) alloy vs sb24's pure-Fe metal, the latter at")
println("         bulk R = 0, i.e. the MOST metal sb24 can make without going below the FeO reference.")
println("="^108)
@printf("%6s %7s %14s %14s %10s | %s\n","P GPa","T C","sb24 ppm @R=0","Fux alloy ppm","sb24/Fux","saturation dIW")
for ((P,TC,_,_), fppm) in zip(cases, (675.0, 1200.0, 1600.0, 2200.0, 3000.0, 3800.0))
    o = single_point_minimization(10P, TC, d, X=bulk_at(BSE,0.0), Xoxides=oxFe, sys_in="wt")
    m = 1e6metalfrac(o)
    @printf("%6.1f %7.0f %14.1f %14.1f %10.3f | %+.3f\n", P, TC, m, fppm, m/fppm, o.fO2-iw(TC+273.15,1e4P))
end
println("Fux alloy abundances read off their S4.4/Fig. 8 discussion; treat as order-of-magnitude.")

println()
println("="^108)
println("TEST A — CEILING on bulk Fe3+/SumFe above which sb24 loses metal, vs Fux's modelled bulk")
println("="^108)
@printf("%6s %7s %14s %14s   %s\n","P GPa","T C","metal-out r","Fux bulk r","sb24/Fux")
for (P,TC,rfux,_) in cases
    lo, hi = 0.0, 0.30                       # r=0 is the fully reduced bulk (metal present); bisect UP to where metal disappears
    hasmetal(r) = metalfrac(single_point_minimization(10P, TC, d, X=bulk_at(BSE,r), Xoxides=oxFe, sys_in="wt")) > 1e-6
    if !hasmetal(lo)
        @printf("%6.1f %7.0f %14s %14.3f   no metal even at r=0 (fully reduced bulk)\n", P, TC, "none", rfux); continue
    end
    for _ in 1:12
        mid = 0.5(lo+hi); hasmetal(mid) ? (lo = mid) : (hi = mid)
    end
    @printf("%6.1f %7.0f %14.3f %14.3f   %.2f\n", P, TC, lo, rfux, lo/rfux)
end
Finalize_MAGEMin(d)
