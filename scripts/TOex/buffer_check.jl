# Verification of TOex.jl's reference buffers (Section 4c of src/TOex.jl).
#
# Three claims are checked, in order of how much they matter:
#   A. the hard-coded SLB (2024) FMQ reproduces MAGEMin's own dQFM. The formulation is exact — with
#      the full-precision HeFESTo parameters it matches to 1e-5 — so the residual here is purely the
#      rounding in the paper's printed Table C1, and must stay under ~0.01 log units;
#   B. a metal-BEARING assemblage sits at ΔIW ≤ 0 (thermodynamically mandatory — with Fe⁰ present
#      the rock is ON its iron-wüstite equilibrium), which the old empirical buffer violated;
#   C. how far the retired Frost (1991) curves sit from the sb24 ones, and in what shape.
using StagPP, MAGEMin_C, Printf

const IW  = StagPP.logfO2_IW
const FMQ = StagPP.logfO2_FMQ
const mmFeO, mmFe, mmO = 71.85, 55.85, 16.0
oxFe = ["SiO2","MgO","Fe","CaO","Al2O3","Na2O","Cr2O3","O"]
BSE  = [45.00, 37.80, 8.05, 3.55, 4.45, 0.36, 0.38, 0.0]
function bulk_at(X, r)
    nFe = X[3]/mmFeO; x = copy(X); x[3] = nFe*mmFe; x[8] = nFe*(1 + 0.5r)*mmO
    return 1e2 .* x ./ sum(x)
end
metalfrac(o) = sum(o.ph_frac_wt[i] for i in eachindex(o.ph)
                   if lowercase(o.ph[i]) in ("feg","fea","fee"); init=0.0)

d = Initialize_MAGEMin("sb24", verbose=false)
fail = 0

println("="^104)
println("A. sb24 Gibbs FMQ vs MAGEMin's own (fO2 - dQFM)")
println("="^104)
@printf("%7s %7s %14s %14s %12s\n", "P GPa", "T K", "QFM MAGEMin", "QFM sb24 Gibbs", "diff")
for (P, TC, r) in ((1.0,1200.0,0.03), (5.0,1420.0,0.036), (8.0,1440.0,0.052),
                   (12.0,1490.0,0.072), (16.0,1590.0,0.110), (20.0,1650.0,0.12))
    T = TC + 273.15
    o = single_point_minimization(10P, TC, d, X=bulk_at(BSE,r), Xoxides=oxFe, sys_in="wt")
    qm, qg = o.fO2 - o.dQFM, FMQ(T, 1e4P)
    abs(qm - qg) > 1e-2 && (global fail += 1)
    @printf("%7.1f %7.0f %14.4f %14.4f %+12.6f\n", P, T, qm, qg, qm - qg)
end

println()
println("="^104)
println("B. ΔIW where native iron is in the assemblage — must be ≤ 0")
println("="^104)
@printf("%7s %7s %8s %11s %11s %13s %13s\n",
        "P GPa", "T K", "metal%", "fO2", "ΔIW sb24", "ΔIW Frost", "verdict")
println("('Frost ✗' marks a row where the retired empirical IW puts a metal-bearing rock ABOVE the buffer)")
for (P, TC) in ((3.0,1350.0), (5.0,1420.0), (8.0,1440.0), (12.0,1490.0), (16.0,1590.0))
    T = TC + 273.15
    o = single_point_minimization(10P, TC, d, X=bulk_at(BSE, 0.0), Xoxides=oxFe, sys_in="wt")
    mf = 1e2metalfrac(o)
    dn = o.fO2 - IW(T, 1e4P)
    df = o.fO2 - StagPP.logfO2_IW_frost(T, 1e4P)
    hasmetal = mf > 1e-4
    bad = hasmetal && dn > 0
    bad && (global fail += 1)
    @printf("%7.1f %7.0f %8.4f %11.3f %+11.3f %+13.3f %13s\n",
            P, T, mf, o.fO2, dn, df,
            bad                     ? "sb24 VIOLATES" :
            hasmetal && df > 0      ? "ok (Frost ✗)"  :
            hasmetal                ? "ok"            : "no metal")
end

println()
println("="^104)
println("C. the retired empirical curves against sb24's")
println("="^104)
@printf("%7s %7s | %11s %11s %10s | %11s %11s %10s\n",
        "P GPa", "T K", "FMQ sb24", "FMQ Frost", "Frost-sb24", "IW sb24", "IW Frost", "Frost-sb24")
for (P, T) in ((1e-3,1000.0), (1e-3,1400.0), (1.0,1500.0), (3.0,1600.0), (5.0,1700.0),
               (8.0,1800.0), (12.0,1900.0), (16.0,2000.0), (20.0,2100.0))
    Pb = 1e4P
    fs, ff = FMQ(T, Pb), StagPP.logfO2_FMQ_frost(T, Pb)
    is, if_ = IW(T, Pb),  StagPP.logfO2_IW_frost(T, Pb)
    @printf("%7.3f %7.0f | %11.4f %11.4f %+10.3f | %11.4f %11.4f %+10.3f\n",
            P, T, fs, ff, ff - fs, is, if_, if_ - is)
end
println("\nFMQ: the Frost error GROWS with pressure (extrapolated 1-bar assemblage).")
println("IW:  the Frost offset is near-CONSTANT — a calibration difference, not an extrapolation one.")

Finalize_MAGEMin(d)
println("\n", fail == 0 ? "buffer_check: PASS" : "buffer_check: $fail FAILURE(S)")
fail == 0 || error("buffer_check failed")
