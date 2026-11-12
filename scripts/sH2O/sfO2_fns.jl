# Kernels for sfO2.jl. Kept as functions only where a global-scope loop would be type-unstable or allocate.
using ForwardDiff, LsqFit

const Rgas = 8.31446261815324   # J/mol/K

# ---- Phase classification. Any returned phase outside these sets stops the run.
const FLUIDS = Set(["fl", "flc", "H2O", "aH2O", "DEW", "liq"])
const NAMLAW = Dict("ol" => :ol, "opx" => :opx, "pig" => :opx, "aug" => :cpx, "dio" => :cpx, "omph" => :cpx, "cpx" => :cpx,
                    "Na-cpx" => :cpx, "g" => :grt, "gtmj" => :grt, "coe" => :coe, "q" => :coe, "qtz" => :coe, "trd" => :coe,
                    "pl4tr" => :pl, "pl" => :pl, "plg" => :pl,
                    "wa" => :wa, "ri" => :ri, "hpcpx" => :hpx, "ak" => :ak, "st" => :st, "stv" => :st, "apbo" => :apbo,
                    "pv" => :pv, "ppv" => :ppv, "mw" => :mw, "fper" => :mw, "cf" => :cf, "capv" => :capv, "crst" => :crst)
const DRY    = Set(["spi", "mgt", "hem", "abc", "ky", "cm", "spl", "ru", "ilm", "sph", "cor", "sill", "and", "mt",
                    "nal", "sp", "feg", "fea", "fee", "neph", "ne", "crn"])
# Fugacity exponent n per NAM law, c(a) = c(1) a^n: 1 for (2H)Mg or molecular water, 1/2 for single OH (Keppler & Bolfan-
#  Casanova 2006); constant-capacity laws take 1. opx: En term 1, Al term 1/2
const NAM_N = Dict(:ol => 1.0, :opx => 1.0, :opxAl => 0.5, :cpx => 0.5, :grt => 0.5, :coe => 1.0, :pl => 1.0, :wa => 1.0,
                   :ri => 1.0, :hpx => 0.5, :ak => 1.0, :st => 0.5, :apbo => 0.5, :pv => 0.5, :ppv => 0.5, :mw => 1.0,
                   :cf => 1.0, :capv => 1.0, :crst => 1.0)

# ---- Pitzer & Sterner (1994) H2O, fugacity in GPa; P capped at 100 GPa like sf_H2O.jl
const PS_C = [
    0.0               0.0             0.24657688e+6        0.51359951e+2        0.0                 0.0
    0.0               0.0             0.58638965e+0       -0.28646939e-2        0.31375577e-4       0.0
    0.0               0.0            -0.62783840e+1        0.14791599e-1        0.35779579e-3       0.15432925e-7
    0.0               0.0             0.0                 -0.42719875e+0       -0.16325155e-4       0.0
    0.0               0.0             0.56654978e+4       -0.16580167e+2        0.76560762e-1       0.0
    0.0               0.0             0.0                  0.10917883e+0        0.0                 0.0
    0.38878656e+13   -0.13494878e+9   0.30916564e+6        0.75591105e+1        0.0                 0.0
    0.0               0.0            -0.65537898e+5        0.18810675e+3        0.0                 0.0
    -0.14182435e+14   0.18165390e+9  -0.19769068e+6       -0.23530318e+2        0.0                 0.0
    0.0               0.0             0.92093375e+5        0.12246777e+3        0.0                 0.0]

function fH2O(PGPa::Float64, T::Float64)
    c = PS_C * [T^-4, T^-2, T^-1, 1.0, T, T^2]
    RT = 1e6Rgas*T                                            # cm3 Pa/mol
    eos(ρ) = (ρ + c[1]*ρ^2 - ρ^2*((c[3] + 2c[4]*ρ + 3c[5]*ρ^2 + 4c[6]*ρ^3) / (c[2] + c[3]*ρ + c[4]*ρ^2 + c[5]*ρ^3 + c[6]*ρ^4)^2)
              + c[7]*ρ^2*exp(-c[8]*ρ) + c[9]*ρ^2*exp(-c[10]*ρ))*RT
    P = 1e9min(PGPa, 100.0); ρ, lb, ub = 0.1, 1e-3, 1e3       # mol/cm3, Newton with bisection fallback
    for _ in 1:300
        Pg = eos(ρ); abs(Pg - P) < 1e-6 && break
        ρn = ρ - (Pg - P)/ForwardDiff.derivative(eos, ρ)
        Pg > P ? (ub = ρ) : (lb = ρ)
        ρ = (!isfinite(ρn) || ρn < lb || ρn > ub) ? 0.5(lb + ub) : ρn
    end
    return exp(log(ρ) + c[1]*ρ + (1/(c[2] + c[3]*ρ + c[4]*ρ^2 + c[5]*ρ^3 + c[6]*ρ^4) - 1/c[2])
               - c[7]/c[8]*(exp(-c[8]*ρ) - 1) - c[9]/c[10]*(exp(-c[10]*ρ) - 1) + P/(ρ*RT) + log(RT) - 1)/1e9
end

# ---- In-house fits of the laws without a published parameter set; run once, data as in sf_H2O.jl
function nam_fits()
    kGPa(n) = (x, p) -> [p[1]*fH2O(x[i,1], x[i,2])^n*exp(-(p[2] + p[3]*x[i,1])/(Rgas*x[i,2])) for i in axes(x, 1)]
    kPa(n)  = (x, p) -> [p[1]*fH2O(x[i,1], x[i,2])^n*exp(-(p[2] + 1e9p[3]*x[i,1])/(Rgas*x[i,2])) for i in axes(x, 1)]
    box(P, T) = (extrema(P), extrema(T))
    # wadsleyite, Dong + Litasov data, wt%
    Pw = [14.5, 15.2, 15.2, 15.5, 16, 16, 17, 18.5, 19, 19, 20, 20, 20, 14.0, 14.0, 14.0, 15.0, 15.5, 16.0, 18.5, 20.0, 20.0]
    Tw = [1873., 1673, 1873, 1673, 1773, 1773, 1673, 1673, 1723, 1923, 1773, 1873, 2173, 1473, 1573, 1673, 1573, 1673, 1773, 1573, 1773, 1873]
    Ww = [0.6, 2.49, 0.56, 1.06, 0.47, 1., 1.15, 1.31, 0.85, 0.35, 0.95, 0.8, 0.4, 2.07, 1.83, 0.81, 1.02, 0.58, 0.55, 1.13, 0.52, 0.44]
    wa = (p=curve_fit(kGPa(1.0), [Pw Tw], Ww, [1e-8, -1e4, 1e4]).param, box=box(Pw, Tw))
    # ringwoodite, Dong 2020 data, wt%
    Pr = [19., 21., 21., 23., 23., 23.]; Tr = [1823., 1823., 1973., 1723., 1800., 2000.]; Wr = [1.16, 1.04, 0.89, 1.60, 1.74, 0.87]
    ri = (p=curve_fit(kGPa(1.0), [Pr Tr], Wr, [1e-2, -1e4, 1e4]).param, box=box(Pr, Tr))
    # high-P clinoenstatite, ppm
    Ph = [8., 10, 10.2, 11.6, 13.4, 13.4]; Th = [1300., 1100, 1400, 1400, 1300, 1400] .+ 273.15; Wh = [0.12, 0.10, 0.19, 0.16, 0.49, 0.35] .* 1e4
    hpx = (p=curve_fit(kPa(0.5), [Ph Th], Wh, [1e-2, -1e2, 1e-9]).param, box=box(Ph, Th))
    # CaCl2-type stishovite, wt%
    Tc = [2900., 2950, 3000, 3050, 3250, 3300, 3350, 3400, 3400, 3450, 3500]; Pc = [25., 42., 41., 29, 58., 31., 51., 42., 59., 50., 100]
    Wc = [0.51, 0.79, 1.27, 0.82, 3.38, 0.76, 1.18, 1.51, 1.85, 1.13, 2.58]
    cacl = (p=curve_fit(kGPa(0.5), [Pc Tc], Wc, [1e-2, 1e4, 3e3]).param, box=box(Pc, Tc))    # ones(3) start lands in a negative-A minimum
    # alpha-PbO2 silica (seifertite), wt%
    Ta = [3650., 3800., 3850., 3850., 4100.]; Pa = [143., 144., 141., 128., 142.]; Wa = [1.99, 1.97, 1.84, 1.94, 1.74]
    apbo = (p=curve_fit(kGPa(0.5), [Pa Ta], Wa, [1e-6, -1e4, 1e3]).param, box=box(Pa, Ta))
    # post-perovskite / bridgmanite water ratio, with and without Al
    dmod(x, p) = @. exp((p[1] + p[2]*x[:,1] - x[:,2]*(p[3] + p[4]*x[:,1] + p[5]*x[:,2] + p[6]))/x[:,2])
    Tal = [1800., 2750, 3550, 1850, 2620, 3300, 3900, 1890, 2500, 3100, 3650, 1950, 2420, 2800, 3350, 3750, 1910, 2380, 2750, 3200, 3510, 4000]
    Pal = [90., 90, 90, 100, 100, 100, 100, 110, 110, 110, 110, 120, 120, 120, 120, 120, 130, 130, 130, 130, 130, 130]
    Dal = [0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.4, 0.1, 0.2, 0.3, 0.4, 0.1, 0.2, 0.3, 0.4, 0.5, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    Tno = [1080., 1200, 1310, 1670, 2820, 1075, 1130, 1230, 1380, 1700, 2900, 1030, 1100, 1170, 1240, 1420, 1770, 3000,
           1050, 1120, 1190, 1250, 1450, 1800, 3030, 1060, 1125, 1200, 1260, 1470, 1850, 3200]
    Pno = [90., 90, 90, 90, 90, 100, 100, 100, 100, 100, 100, 110, 110, 110, 110, 110, 110, 110, 120, 120, 120, 120, 120, 120, 120, 130, 130, 130, 130, 130, 130, 130]
    Dno = [12.5, 10.0, 7.5, 5.0, 2.5, 15.0, 12.5, 10.0, 7.5, 5.0, 2.5, 17.5, 15.0, 12.5, 10.0, 7.5, 5.0, 2.5, 17.5, 15.0, 12.5, 10.0, 7.5, 5.0, 2.5, 17.5, 15.0, 12.5, 10.0, 7.5, 5.0, 2.5]
    p0 = [1e-6, 1e-6, 1e-6, 1e-2, 1e-5, 1e-2]
    ppv = (al=curve_fit(dmod, [Pal Tal], Dal, p0).param, noal=curve_fit(dmod, [Pno Tno], Dno, p0).param)
    return (wa=wa, ri=ri, hpx=hpx, cacl=cacl, apbo=apbo, ppv=ppv)
end

# Bridgmanite H2O (wt%) against its own Al2O3 (wt%): Litasov et al. 2003 pure MgSiO3, peridotite Al-pv, MORB pv; linear, flat beyond
const BDG_AL = (0.0, 5.5, 15.0); const BDG_W = (1e-2, 0.14, 7.5e-3)
bdg_h2o(a) = (a = clamp(a, BDG_AL[1], BDG_AL[3]); q = a <= BDG_AL[2] ? 1 : 2; BDG_W[q] + (BDG_W[q+1] - BDG_W[q])*(a - BDG_AL[q])/(BDG_AL[q+1] - BDG_AL[q]))

# CaCl2-type transition of stishovite (Umemoto 2006 points), P in GPa against T in K
const CACL_T = [563.2, 706.0, 833.9, 945.3, 1055., 1141., 1212., 1282., 1347., 1400., 1476., 1544., 1652., 1894.]
const CACL_P = [19.00, 22.80, 26.60, 30.40, 34.20, 38.00, 41.80, 45.60, 49.40, 53.20, 57.00, 60.80, 64.60, 68.40]

# ---- NAM saturation laws at a_H2O = 1, wt%; (P, T) clamped to each law's calibration box.
#  alpv: bridgmanite from its Al2O3 (BDG law, in assemble!) instead of the 30 ppm constant
function nam_laws(P::Float64, T::Float64, F; alpv=true)
    keppler(f, Pc, Tc, A, n, ΔH, ΔV) = 1e-4A*(1e4f)^n*exp(-(ΔH + ΔV*1e9Pc)/(Rgas*Tc))          # f GPa; A ppm/bar^n, ΔH J/mol, ΔV m3/mol
    fitGPa(f, n) = (Pc = clamp(P, f.box[1]...); Tc = clamp(T, f.box[2]...); f.p[1]*fH2O(Pc, Tc)^n*exp(-(f.p[2] + f.p[3]*Pc)/(Rgas*Tc)))
    Po, To = clamp(P, 0.01, 14.2), clamp(T, 1200.0, 2400.0)
    ol = 1e-4exp(0.5*0.6447*log(fH2O(Po, To)) + 4905.5403/To)                                   # Dong 2020 fit, fug in GPa
    Px, Tx = clamp(P, 0.01, 14.0), clamp(T, 500.0, 1500.0); fx = fH2O(Px, Tx)
    opx   = keppler(fx, Px, Tx, 0.01354, 1.0, -4563.0, 12.1e-6)                                 # Mierdel et al. 2007, pure En
    Ta = max(Tx, 1073.0)                                                                         # Al term: its data start at 800 °C
    opxAl = keppler(Ta == Tx ? fx : fH2O(Px, Ta), Px, Ta, 0.042, 0.5, -79685.0, 11.3e-6)        # Mierdel et al. 2007, Al-saturated term
    Pc, Tc = clamp(P, 1e-4, 10.0), clamp(T, 873.0, 973.0); fc = fH2O(Pc, Tc)
    di = keppler(fc, Pc, Tc, 2.15, 0.5, 0.0, 7.43e-6)                                           # Bromiley et al. 2004
    jd = keppler(fc, Pc, Tc, 7.144, 0.5, 0.0, 8.019e-6)                                         # Bromiley & Keppler 2004
    Pg, Tg = clamp(P, 0.01, 10.0), clamp(T, 700.0, 1500.0)
    grt = keppler(fH2O(Pg, Tg), Pg, Tg, 0.866, 0.5, 2580.0, 5.71e-6)                            # Lu & Keppler 1997 form, bar and Pa
    coe = 1e-4max(-49.0 + 6.0clamp(P, 1e-4, 10.0) + 0.06clamp(T, 800.0, 1073.0), 0.0)
    wa, ri = fitGPa(F.wa, 1.0), fitGPa(F.ri, 1.0)
    h = F.hpx; Ph, Th = clamp(P, h.box[1]...), clamp(T, h.box[2]...)
    hpx = 1e-4h.p[1]*fH2O(Ph, Th)^0.5*exp(-(h.p[2] + 1e9h.p[3]*Ph)/(Rgas*Th))
    # stishovite: Chen et al. 2021 below the CaCl2 transition (Umemoto 2006 points, linear beyond), CaCl2-type fit above
    Ts, Ps = CACL_T, CACL_P
    k = clamp(searchsortedlast(Ts, T), 1, length(Ts) - 1)
    Pb = Ps[k] + (Ps[k+1] - Ps[k])*(T - Ts[k])/(Ts[k+1] - Ts[k])
    Pst, Tst = clamp(P, 0.01, 60.01), clamp(T, 800.0, 1500.0)
    st = P <= Pb ? 1e-4*24.0*fH2O(Pst, Tst)^0.5*exp(-(-3065.0 + 4.29e-6*1e9Pst)/(Rgas*Tst)) : fitGPa(F.cacl, 0.5)
    apbo = fitGPa(F.apbo, 0.5)
    dm(p, Pd, Td) = exp((p[1] + p[2]*Pd - Td*(p[3] + p[4]*Pd + p[5]*Td + p[6]))/Td)
    Pd, Td = clamp(P, 90.0, 135.0), clamp(T, 1000.0, 4000.0)
    pv = 3e-3                                                                                    # cf and Al_pv = false reference
    # ppv/pv water ratios: the fitted Townsend et al. (2016) D is C_pv/C_ppv, so ppv holds pv/D
    return (ol=ol, opx=opx, opxAl=opxAl, di=di, jd=jd, grt=grt, coe=coe, wa=wa, ri=ri, hpx=hpx, ak=ri/21.0, st=st, apbo=apbo,
            pv=alpv ? BDG_W[1] : pv, rAl=1/dm(F.ppv.al, Pd, Td), rNo=1/dm(F.ppv.noal, Pd, Td), mw=5e-3, cf=pv, capv=0.0, crst=1e-4,
            pl=0.0)                                                                              # plagioclase: set per use (sH2O_moon.jl)
end
function nam_laws(P::AbstractVector, T::AbstractVector, F; alpv=true)                           # threaded over a node list
    L = Vector{typeof(nam_laws(P[1], T[1], F; alpv=alpv))}(undef, length(P))
    Threads.@threads for k in eachindex(P); L[k] = nam_laws(P[k], T[k], F; alpv=alpv); end
    return L
end

# ---- One composition over the (P,T) grid: water in hydrous solids + NAM laws, per unit solid mass (wt%).
#  With sub, also stores the MSH-type phases' moles of MgO+FeO and SiO2, their mass and their NAM water (dhms! input).
#  a: water activity of the capacity; NAM laws scale as a^n (NAM_N), hydrous phases come from the run itself
function assemble!(cap, hyd, host, unknown, out, L; sub=nothing, alpv=true, a=1.0)
    lk = ReentrantLock()
    Threads.@threads for k in eachindex(out)
        o = out[k]; iH = findfirst(==("H2O"), o.oxides); iAl = findfirst(==("Al2O3"), o.oxides)
        iSi = findfirst(==("SiO2"), o.oxides); iMg = findfirst(==("MgO"), o.oxides); iFe = findfirst(o.oxides) do x; x == "Fe" || x == "FeO"; end
        mFe = (iFe !== nothing && o.oxides[iFe] == "Fe") ? 55.85 : 71.85; nM = 0.0; nS = 0.0; mMS = 0.0; wMS = 0.0
        nss = length(o.SS_vec); ms = 0.0
        for j in eachindex(o.ph); o.ph[j] in FLUIDS || (ms += o.ph_frac_wt[j]); end
        wh = 0.0; wn = 0.0; best = 0.0; bname = ""
        for j in eachindex(o.ph)
            ph = o.ph[j]; ph in FLUIDS && continue
            x = o.ph_frac_wt[j]/ms
            cw = iH === nothing ? 0.0 : (j <= nss ? o.SS_vec[j].Comp_wt[iH] : o.PP_vec[j - nss].Comp_wt[iH])
            if cw > 1e-6
                w = 100cw*x; wh += w
            elseif haskey(NAMLAW, ph)
                law = NAMLAW[ph]; l = L[k]
                if law === :cpx
                    em = o.SS_vec[j].emNames; ef = o.SS_vec[j].emFrac
                    ij = findfirst(e -> e == "jd" || e == "jdm", em); xjd = ij === nothing ? 0.0 : ef[ij]
                    w = x*(xjd*l.jd + (1.0 - xjd)*l.di)
                elseif law === :opx
                    w = x*(l.opx*a^NAM_N[:opx] + min(o.SS_vec[j].Comp[iAl]/0.1005, 1.0)*l.opxAl*a^NAM_N[:opxAl])
                elseif law === :pv && alpv
                    w = x*bdg_h2o(100o.SS_vec[j].Comp_wt[iAl])
                elseif law === :ppv                                         # coexisting-pv reference at the ppv's own Al2O3
                    s = o.SS_vec[j]; xal = s.Comp[iAl]; ref = alpv ? bdg_h2o(100s.Comp_wt[iAl]) : l.pv
                    w = x*ref*(xal*l.rAl + (1.0 - xal)*l.rNo)
                else
                    w = x*getfield(l, law)
                end
                law === :opx || (w *= a^NAM_N[law])
                wn += w
            else
                w = 0.0; ph in DRY || lock(() -> push!(unknown, ph), lk)
            end
            if sub !== nothing && ph in MSH_SET
                cc = j <= nss ? o.SS_vec[j].Comp_wt : o.PP_vec[j - nss].Comp_wt
                nM += x*(cc[iMg]/40.30 + cc[iFe]/mFe); nS += x*cc[iSi]/60.08; mMS += x; wMS += w
            end
            w > best && (best = w; bname = ph)
        end
        cap[k] = wh + wn; hyd[k] = wh; host[k] = bname
        sub === nothing || (sub[1, k] = nM; sub[2, k] = nS; sub[3, k] = mMS; sub[4, k] = wMS)
    end
end

# ---- MgO-SiO2-H2O with dense hydrous magnesium silicates: Komabayashi & Omori (2006, PEPI 156) Tables 3-4.
#  Row: MgO, SiO2, H2O per formula; ΔfH (kJ), S; Cp a..g; V0 (J/bar), K'0, K'T; β0-β3, α0-α3 as printed, scaled here
ko_row(m, s, h, H, S, a, b, c, d, e, f, g, V0, Kp, KpT, β0, β1, β2, β3, α0, α1, α2, α3) =
    (m=Float64(m), s=Float64(s), h=Float64(h), H=1e3H, S=Float64(S), a=Float64(a), b=1e-5b, c=Float64(c), d=1e-5d, e=Float64(e),
     f=Float64(f), g=Float64(g), V0=Float64(V0), Kp=Float64(Kp), KpT=1e-4KpT, β0=1e-7β0, β1=1e-10β1, β2=1e-15β2, β3=1e-17β3,
     α0=1e-5α0, α1=1e-9α1, α2=Float64(α2), α3=Float64(α3))
const KO_NAME = ("phA", "cHu", "phE", "phD", "shB", "hyWad", "hyRin", "br", "fo", "hen", "wd", "rg", "aki", "pv", "per", "stv", "phH")
const KO = (
    ko_row(7, 2, 3, -7139.989, 341, 648.79, 6869.1, -27008600, 0, 2347599591, -327, 9301.7, 15.44, 5.7, 0, 9.4816, 2.5538, -3.1507, 4.0619, 6.8482, -1.7420, -0.013093, 0.9623375),
    ko_row(9, 4, 1, -9608.049, 461.6, 782.53, 5471.7, -40286200, 0, 4695199864, -654, 19173.9, 19.785, 4, 0, 7.5706, 1.3905, 12.4692, 0.9887, 5.2104, -1.1913, -0.011055, 1.0429425),
    ko_row(2.33, 1.26, 1.3, -2979.229, 166, 288.8505, 2640.909, -10337160, 0, 748030423, -435.33, 2335, 6.516, 5, 0, 9.9404, 2.6774, -3.3071, 4.2587, 6.8482, -1.7420, -0.013093, 0.9623375),
    ko_row(1.18, 1.58, 1.13, -2456.7, 78, 254.2082, 1892.094, -8907106, 0, 721473646, -545.89, -8.383, 4.652, 4.1, 0, 5.6696, 1.1502, 8.5807, 1.0285, 2.319, 9.04, -0.003966, 0.7496),
    ko_row(10, 3, 2, -9648.078, 255, 827.28, 7125, -33587400, 0, 2748529727, -1036.5, 19059.8, 18.615, 5.8, 0, 6.5213, 1.6086, 3.4921, 2.1713, 1.8609, 0, 0, 0),
    ko_row(1.88, 1, 0.12, -2121.097, 84, 172.9, 1129, -1077100, 0, -218700000, 0, -13477, 4.104, 4.2, 0, 6.0432, 1.3283, 7.9979, 1.4188, 1.8609, 0, 0, 0),
    ko_row(1.91, 0.99, 0.11, -2120.883, 73.2, 158.56, 1220, -12297000, 0, 1484100000, 0, 7971.9, 3.991, 5, 0, 5.7947, 1.6369, -5.0022, 2.8422, 1.6149, 0, 0, 0),
    ko_row(1, 0, 1, -925.5, 64.4, 103.01, 1653.3, -2746200, 0, -136.5, 0, -114.1, 2.463, 4.712, 10, 17.5100, -0.0549, 204.4500, 11.0500, 8.7540, 18.0900, -0.02218, 0.5843),
    ko_row(2, 1, 0, -2174.14, 95.6, 165.8, 1855, -3971000, 0, 286100000, 0, -5610, 4.367, 5.2, 2, 7.4270, 1.2400, 6.9000, 1.7020, 2.0100, 13.9000, 0.001627, -0.338),
    ko_row(1, 1, 0, -1544.823, 65.3, 144.45, 188.2, -1350000, 0, 461200000, 0, -19380, 3.0918, 3.6, 1.5, 8.8958, 1.3547, 31.7162, 1.1524, 3.8770, 4.5002, 0.0003439, -1.723621),
    ko_row(2, 1, 0, -2140.5, 95, 172.9, 1129, -1077100, 0, -218700000, 0, -13477, 4.054, 4.3, 3, 5.5128, 0.9202, 8.8490, 1.1529, 2.3190, 9.04, -0.003966, 0.7496),
    ko_row(2, 1, 0, -2132.1, 90.6, 158.56, 1220, -12297000, 0, 1484100000, 0, 7971.9, 3.965, 4.3, 6, 5.07778, 1.3371, -29.854, 2.7822, 1.225, 11.04, 0.002496, -0.511),
    ko_row(1, 1, 0, -1494.1, 58.3, 125.5, 530.3, -6282000, 0, 766500000, 0, -3040, 2.635, 4, 3, 4.504, 0.82987, -3.8777, 1.2095, 2.2704, 6.8209, -0.001809, 0.3851),
    ko_row(1, 1, 0, -1449.5, 63.6, 137.9, 156.4, -3191000, 0, 719200000, 0, -14590, 2.4447, 4.1, 0, 3.7184, 0.23496, 22.508, 0, 4.903, -2.307, -0.01812, 2.543),
    ko_row(1, 0, 0, -601.49, 26.94, 45.49, 477.3, -2160000, 0, 174400000, 0, 2411, 1.125, 4.17, 0.4, 5.875, 1.101, 13.7, 0.48, 3.64, 8.35, 0.00085, -0.95),
    ko_row(0, 1, 0, -874, 24.642, 85.78, 0, -3605000, 0, 451100000, -345.5, 0, 1.401, 6, 10, 2.954, 0.8961, -32.9, 2.331, 0.23, 12, 0.0062, -1.13))
const KO_W    = ko_row(0, 0, 1, -241.814, 188.83, 7.368, 2746.8, -223160, -0.48117, 0, 361.74, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)  # Robie et al. 1978
const KO_DHMS = (true, true, true, true, true, false, false, true, false, false, false, false, false, false, false, false, true)
const KO_LAW  = (:none, :none, :none, :none, :none, :wa, :ri, :none, :ol, :hen, :wa, :ri, :ak, :pv, :mw, :st, :none)  # NAM law per phase
const MSH_SET = Set(["ol", "wa", "ri", "opx", "hpcpx", "ak", "pv", "mw", "st", "coe", "qtz"])                    # sb24 phases re-equilibrated in MSH

# Apparent Gibbs energy (J/mol) at P (GPa), T (K): closed-form Cp integrals from 298.15 K, Murnaghan ∫V dP from 1 bar
function g_ko(p, P::Float64, T::Float64)
    T0 = 298.15
    H = p.H + p.a*(T - T0) + p.b/2*(T^2 - T0^2) - p.c*(1/T - 1/T0) + p.d/3*(T^3 - T0^3) - p.e/2*(T^-2 - T0^-2) +
        2p.f*(sqrt(T) - sqrt(T0)) + p.g*log(T/T0)
    S = p.S + p.a*log(T/T0) + p.b*(T - T0) - p.c/2*(T^-2 - T0^-2) + p.d/2*(T^2 - T0^2) - p.e/3*(T^-3 - T0^-3) -
        2p.f*(T^-0.5 - T0^-0.5) - p.g*(1/T - 1/T0)
    p.V0 == 0.0 && return H - T*S
    V1 = p.V0*exp(p.α0*(T - T0) + p.α1/2*(T^2 - T0^2) + p.α2*log(T/T0) - p.α3*(1/T - 1/T0))
    K  = 1/(p.β0 + p.β1*T + p.β2*T^2 + p.β3*T^3)                                         # bar
    Kp = p.Kp + p.KpT*(T - T0)*log(T/T0)
    return H - T*S + V1*K/(Kp - 1)*((1 + Kp*(1e4P - 1)/K)^(1 - 1/Kp) - 1)
end

# Brodholt & Wood (1993) H2O, P in bar, V in cm3/mol: V(P) by safeguarded Newton; ∫P dV closed form (b is linear in V)
const BW = (b0=-3.90463e-2, b1=-0.991078, c=3.64905e4, d=-1.02451e7, e=-1.79681e8, f=2.18437e9, R=83.1446)
bw_a(T) = -582468 - 3038.79T - 9.24574e-3T^2 + 3.02674e9/T^2
bw_P(V, T) = (q = BW; q.R*T/((1 - q.b1)V - q.b0) - bw_a(T)/(sqrt(T)*V*((1 + q.b1)V + q.b0)) + q.c/V + q.d/V^2 + q.e/V^3 + q.f/V^4)
bw_F(V, T) = (q = BW; q.R*T/(1 - q.b1)*log((1 - q.b1)V - q.b0) - bw_a(T)/(sqrt(T)*q.b0)*log(V/((1 + q.b1)V + q.b0)) +
              q.c*log(V) - q.d/V - q.e/(2V^2) - q.f/(3V^3))
function bw_V(P::Float64, T::Float64)
    V, lb, ub = 12.0, 4.5, 80.0
    for _ in 1:200
        Pg = bw_P(V, T); abs(Pg - P) < 1e-9P && break
        Vn = V - (Pg - P)/ForwardDiff.derivative(v -> bw_P(v, T), V)
        Pg > P ? (lb = V) : (ub = V)
        V = (!isfinite(Vn) || Vn <= lb || Vn >= ub) ? 0.5(lb + ub) : Vn
    end
    return V
end

# μ_H2O (J/mol) of the MSH minimisation at a = 1: Robie 1-bar G, Pitzer & Sterner to 1 GPa, Brodholt & Wood above
function mu_h2o(P::Float64, T::Float64)
    g = g_ko(KO_W, 0.0, T) + Rgas*T*log(1e4fH2O(clamp(P, 1e-4, 1.0), T))
    P <= 1.0 && return g
    V1, V2 = bw_V(1e4, T), bw_V(1e4P, T)
    return g + 0.1*(1e4P*V2 - 1e4V1 - (bw_F(V2, T) - bw_F(V1, T)))                    # cm3 bar -> J
end

# ln a_H2O: thin-plate spline through the activities K&O fitted to fluid-present experiments (P GPa, T °C, a), a = 1 anchors
#  below 700 °C and to 9 GPa (Komabayashi et al. 2005), plus rows X; held flat outside 5-Pmax GPa, 500-1600 °C
const KO_AH2O = [5.0 550 1; 6.0 620 1; 7.1 690 1; 9.0 800 1; 11.7 580 1; 14.6 550 1; 21.55 728 1; 25.0 750 1;
                 12.2 1024 0.68; 12.8 1156 0.68; 13.5 1200 0.60; 14.5 1232 0.60; 15.0 900 0.85; 16.6 988 0.85;
                 16.1 1176 0.54; 17.0 1100 0.54; 19.1 1180 0.13; 19.9 1140 0.13; 20.0 1200 0.13; 23.3 1250 0.08;
                 25.0 1200 0.10; 26.3 1224 0.10; 28.0 1116 0.10; 26.0 1400 0.05; 31.0 1040 0.03; 24.8 1560 0.012]
function ah2o_fit(X=zeros(0, 3); Pmax=35.0)
    Q = vcat(KO_AH2O, X); x = Q[:, 1]./5; y = Q[:, 2]./200; n = length(x)
    φ(r) = r > 0 ? r^2*log(r) : 0.0
    A = [[φ(hypot(x[i] - x[j], y[i] - y[j])) for i in 1:n, j in 1:n] [ones(n) x y]; [ones(n) x y]' zeros(3, 3)]
    c = A \ [log.(Q[:, 3]); zeros(3)]
    return (x=x, y=y, w=c[1:n], c=c[n+1:end], Pmax=Pmax, Q=Q)
end
function ln_ah2o(P::Float64, T::Float64, A)
    x = clamp(P, 5.0, A.Pmax)/5; y = clamp(T - 273.15, 500.0, 1600.0)/200
    s = A.c[1] + A.c[2]*x + A.c[3]*y
    for i in eachindex(A.w); r = hypot(x - A.x[i], y - A.y[i]); r > 0 && (s += A.w[i]*r^2*log(r)); end
    return min(s, 0.0)
end

# Phase H MgSiO2(OH)2: Cp = pv + br - per; phase D's α and K_T(T) shape at K0 = 196 GPa (Nishi 2014, Ohtani 2014 cells).
#  ΔfH, S exact on Nishi's phD = phH + St line (Fig 2a); effective V0 and K' fitted to his phH-in curve (Fig 2b)
function ko_phH(A; p0=[3.35, 4.0])
    D, br, pv, per, st = KO[4], KO[8], KO[14], KO[15], KO[16]
    s = (D.β0 + D.β1*298.15 + D.β2*298.15^2 + D.β3*298.15^3)^-1/1.9607e6
    row(H, S, V0, Kp) = ko_row(1, 1, 1, H, S, pv.a + br.a - per.a, 1e5(pv.b + br.b - per.b), pv.c + br.c - per.c, 1e5(pv.d + br.d - per.d),
                               pv.e + br.e - per.e, pv.f + br.f - per.f, pv.g + br.g - per.g, V0, Kp, 0, 1e7s*D.β0, 1e10s*D.β1,
                               1e15s*D.β2, 1e17s*D.β3, 1e5D.α0, 1e9D.α1, D.α2, D.α3)
    mu(P, T) = mu_h2o(P, T) + Rgas*T*ln_ah2o(P, T, A)
    dGa(h, P, T) = 1.18g_ko(h, P, T) + 0.40g_ko(st, P, T) - g_ko(D, P, T) - 0.05mu(P, T)   # 0.05 H2O balances K&O's phD
    function hs(V0, Kp)                                                 # G is linear in (ΔfH, S): exact on the two Fig 2a points
        g1 = dGa(row(0, 0, V0, Kp), 45.4, 873.15)/1.18; g2 = dGa(row(0, 0, V0, Kp), 48.6, 1308.15)/1.18
        S = (g2 - g1)/435.0; return row((873.15S - g1)/1e3, S, V0, Kp)
    end
    function dGc(h, P, T)                                               # phH against the phH-free assemblage of its own bulk
        m = mu(P, T); G = ntuple(k -> g_ko(KO[k], P, T) - KO[k].h*m, Val(length(KO)))
        i, j, ni, nj = msh_min(1.0, 1.0, G, KO)
        return g_ko(h, P, T) - m - ni*G[i] - nj*G[j]
    end
    tg = ((28.0, 903.15), (31.5, 1073.15), (34.0, 1173.15))
    p = curve_fit((x, q) -> [dGc(hs(q[1], q[2]), P, T)/5e3 for (P, T) in tg], zeros(3), zeros(3), p0).param
    return hs(p[1], p[2])
end
# DHMS inputs: the a_H2O map extended with the activities that put Nishi et al. (2014) Fig 2 dehydration lines on this model
#  (phase D, K&O R88; phase H -> pv + H2O) and the 17-phase table; iterated, since both calibrations see the water potential
function msh_setup(; Pmax=60.0, n=3)
    D, pv, st = KO[4], KO[14], KO[16]
    aD(h, P, T) = exp(-((236g_ko(pv, P, T) + 80g_ko(st, P, T) - 200g_ko(D, P, T))/226 + mu_h2o(P, T))/(Rgas*T))
    aH(h, P, T) = exp((g_ko(h, P, T) - g_ko(pv, P, T) - mu_h2o(P, T))/(Rgas*T))
    pts = ((30.0, 1187.0, aD), (40.0, 1105.0, aD), (48.6, 1035.0, aD), (40.0, 997.0, aH), (50.0, 1030.0, aH), (57.9, 1056.0, aH))
    A = ah2o_fit(); h = ko_phH(A)
    for _ in 1:n
        X = [q[1] for q in pts]; Y = [q[2] for q in pts]; Z = [q[3](h, q[1], q[2] + 273.15) for q in pts]
        A = ah2o_fit([X Y Z]; Pmax=Pmax); h = ko_phH(A)
    end
    return A, (KO..., h)
end

# Exact minimum of Σ n G* over MgO-SiO2 with H2O at fixed μ: at most two solids, so every pair is solved
function msh_min(NM::Float64, NS::Float64, G, K)
    best = Inf; bi = 0; bj = 0; bni = 0.0; bnj = 0.0
    for i in 1:length(K), j in i+1:length(K)
        mi, si, mj, sj = K[i].m, K[i].s, K[j].m, K[j].s
        det = mi*sj - mj*si; abs(det) < 1e-12 && continue
        ni = (NM*sj - NS*mj)/det; nj = (mi*NS - si*NM)/det
        (ni < -1e-12 || nj < -1e-12) && continue
        g = ni*G[i] + nj*G[j]
        g < best && (best = g; bi = i; bj = j; bni = max(ni, 0.0); bnj = max(nj, 0.0))
    end
    return bi, bj, bni, bnj
end

# DHMS term on sb24 nodes: the (Mg+Fe, Si) of the MSH-type phases (sub, from assemble!) is re-equilibrated with H2O.
#  Where a DHMS forms the node is rebuilt (DHMS water + NAM laws on the MSH anhydrous phases); elsewhere it is untouched.
#  a: cell water activity, added to the map's effective one in μ; the anhydrous member's NAM law scales as a^n
function dhms!(cap, hyd, host, sub, P, T, L, A, K; Pmax=60.0, Tmax=1873.0, PDmax=50.0, a=1.0)
    Threads.@threads for k in eachindex(cap)
        (P[k] > Pmax || T[k] > Tmax || sub[1, k] <= 0.0) && continue
        mu = mu_h2o(P[k], T[k]) + Rgas*T[k]*(ln_ah2o(P[k], T[k], A) + log(a))
        G  = ntuple(i -> (i == 4 && P[k] > PDmax) ? 1e12 : g_ko(K[i], P[k], T[k]) - K[i].h*mu, Val(length(K)))  # no phD beyond PDmax
        i, j, ni, nj = msh_min(sub[1, k], sub[2, k], G, K)
        (KO_DHMS[i] && ni > 0) || (KO_DHMS[j] && nj > 0) || continue
        Mb = (sub[3, k] - 60.08sub[2, k])/sub[1, k]                     # (Mg,Fe)O-part molar mass: carries Fe and minor oxides
        W = 0.0; wN = 0.0; best = 0.0; bname = ""
        for (q, n) in ((i, ni), (j, nj))
            p = K[q]
            if KO_DHMS[q]
                W += n*p.h*18.015; w = 100n*p.h*18.015
            else
                law = KO_LAW[q]; l = L[k]; law === :hen && (law = P[k] < 8.0 ? :opx : :hpx)
                w = n*(p.m*Mb + 60.08p.s)*getfield(l, law)*a^NAM_N[law]; wN += w
            end
            w > best && (best = w; bname = KO_NAME[q])
        end
        cap[k] = (100W + wN + cap[k] - sub[4, k])/(1 + W); hyd[k] = 100W/(1 + W); host[k] = bname
    end
end

# ---- Ported from watercarb/Perple_X_MAGEMin/MAGEMin/s_fO2_dens_metal.jl (2026-09-19 production generator); arithmetic
#  kept operation-for-operation so the fO2, metal, step and density tables stay byte-identical
using Interpolations, StaticArrays
const Nₐ = 6.02214076e23
const MM = (SiO2=60.08, MgO=40.30, FeO=71.85, Fe=55.85, CaO=56.08, Al2O3=101.96, Na2O=61.98, Cr2O3=151.99,
            Fe2O3=159.69, O=16.0, TiO2=79.88, K2O=94.2, H2O=18.015, CO2=44.009, C=12.011)   # reference molar masses
get_Xoxmm(Xox) = [getfield(MM, Symbol(ox)) for ox in Xox]

# ---- Signed Fe3+/ΣFe bulk writer (reference oxidize_bulk, vector output only)
function oxidize_bulk(X, Rv, Xdummy; wt_in=false, wt_out=false, frac=false, FeFormat="Fe_O", oldXox=nothing)
    @assert -2 <= Rv <= 1 "R = Fe³⁺/ΣFe must lie in [-2, 1]"
    (Rv < 0.0) && (FeFormat = "Fe_O")
    entryXox = copy(Xdummy.Xox)
    if !isnothing(oldXox)
        for (i, ox) in enumerate(Xdummy.Xox)
            Xdummy.X[i] = ((ox=="O" && "FeO" ∈ oldXox) || ox=="Fe2O3" || !(ox ∈ oldXox)) ? 0.0 : X[oldXox.==ox][1]
        end
    else
        Xdummy.X .= Vector{Float64}(X)
    end
    Xd = Xdummy.X; Xoxd = Xdummy.Xox; mmd = Xdummy.mm
    wt_in && (Xd ./= mmd)
    n = sum(Xd); Xd ./= n
    if ("Fe" ∈ Xoxd) && ("O" ∈ Xoxd);        idxF, idxO = findfirst(Xoxd.=="Fe"), findfirst(Xoxd.=="O")
    elseif ("Fe2O3" ∈ Xoxd) && ("FeO" ∈ Xoxd); idxF, idxO = findfirst(Xoxd.=="FeO"), findfirst(Xoxd.=="Fe2O3")
    elseif ("FeO" ∈ Xoxd) && ("O" ∈ Xoxd);   idxF, idxO = findfirst(Xoxd.=="FeO"), findfirst(Xoxd.=="O")
    else error("oxidize_bulk: no iron pair in $(Xoxd)")
    end
    nFe = Xd[idxF]
    if FeFormat == "Fe_O";          Xd[idxF] = nFe;          Xd[idxO] = nFe*(1 + 0.5Rv);  Xoxd[idxF] = "Fe";  Xoxd[idxO] = "O"
    elseif FeFormat == "FeO_O";     Xd[idxF] = nFe;          Xd[idxO] = 0.5Rv*nFe;        Xoxd[idxF] = "FeO"; Xoxd[idxO] = "O"
    elseif FeFormat == "Fe2O3_FeO"; Xd[idxF] = (1 - Rv)*nFe; Xd[idxO] = 0.5Rv*nFe;        Xoxd[idxF] = "FeO"; Xoxd[idxO] = "Fe2O3"
    else error("oxidize_bulk: unknown FeFormat \"$FeFormat\"")
    end
    if wt_out
        for i in eachindex(Xd); Xd[i] *= (i==idxF || i==idxO) ? getfield(MM, Symbol(Xoxd[i])) : mmd[i]; end
    end
    n = sum(Xd); Xd ./= n
    !frac && (Xd .*= 1e2)
    Xdummy.Xox .= entryXox
    return copy(Xd)
end

function oxidized_Xox(Xox, FeFormat)
    out = Vector{String}(Xox)
    idxF = findfirst(x -> x in ("FeO", "Fe"), out); idxO = findfirst(x -> x in ("O", "Fe2O3"), out)
    (isnothing(idxF) || isnothing(idxO)) && return out
    if FeFormat == "Fe_O"; out[idxF], out[idxO] = "Fe", "O"
    elseif FeFormat == "FeO_O"; out[idxF], out[idxO] = "FeO", "O"
    elseif FeFormat == "Fe2O3_FeO"; out[idxF], out[idxO] = "FeO", "Fe2O3"
    else error("oxidized_Xox: unknown FeFormat \"$FeFormat\"")
    end
    return out
end

function filter_ox_list(Xin, Xox, newXox; wt_in=false, wt_out=false)
    Oflag1 = ("O"∈Xox)||("Fe2O3"∈Xox); Oflag2 = ("O"∈newXox)||("Fe2O3"∈newXox)
    FeFormat1, FeFormat2 = "", ""
    Oflag1 && (("O"∈Xox)&&("Fe"∈Xox) ? FeFormat1="Fe_O" : ("O"∈Xox)&&("FeO"∈Xox) ? FeFormat1="FeO_O" : "FeO_Fe2O3")
    Oflag2 && (("O"∈newXox)&&("Fe"∈newXox) ? FeFormat2="Fe_O" : ("O"∈newXox)&&("FeO"∈newXox) ? FeFormat2="FeO_O" : "FeO_Fe2O3")
    X = copy(Xin)
    wt_in && (X ./= get_Xoxmm(Xox))
    newX = zeros(length(newXox))
    for (i, ox) in enumerate(newXox); newX[i] = ox ∈ Xox ? X[findfirst(Xox.==ox)] : 0.0; end
    if (FeFormat1=="" && FeFormat2=="") || (FeFormat1==FeFormat2)
    elseif FeFormat1=="Fe_O" && (FeFormat2=="FeO_O" || FeFormat2=="FeO_Fe2O3")
        idxF1 = findfirst(Xox.=="Fe"); idxO1 = findfirst(Xox.=="O"); idxF2 = findfirst(newXox.=="FeO"); idxO2 = findfirst(newXox.=="O")
        newX[idxO2] = X[idxO1] - X[idxF1]; newX[idxF2] = 3X[idxF1] - 2X[idxO1]
    elseif FeFormat1=="FeO_O" && FeFormat2=="FeO_Fe2O3"
        idxF1 = findfirst(Xox.=="FeO"); idxO1 = findfirst(Xox.=="O"); idxF2 = findfirst(newXox.=="FeO"); idxO2 = findfirst(newXox.=="Fe2O3")
        newX[idxO2] = X[idxO1]; newX[idxF2] = X[idxF1] - 2X[idxO1]
    elseif FeFormat1=="FeO_O" && FeFormat2=="Fe_O"
        idxF1 = findfirst(Xox.=="FeO"); idxO1 = findfirst(Xox.=="O"); idxF2 = findfirst(newXox.=="Fe"); idxO2 = findfirst(newXox.=="O")
        newX[idxF2] = X[idxF1]; newX[idxO2] = X[idxF1] + X[idxO1]
    elseif FeFormat1=="FeO_Fe2O3" && FeFormat2=="Fe_O"
        idxF1 = findfirst(Xox.=="FeO"); idxO1 = findfirst(Xox.=="Fe2O3"); idxF2 = findfirst(newXox.=="Fe"); idxO2 = findfirst(newXox.=="O")
        newX[idxF2] = X[idxF1] + 2X[idxO1]; newX[idxO2] = X[idxF1] + 3X[idxO1]
    elseif FeFormat1=="FeO_Fe2O3" && FeFormat2=="FeO_O"
        idxO1 = findfirst(Xox.=="Fe2O3"); idxO2 = findfirst(newXox.=="O"); newX[idxO2] = X[idxO1]
    end
    wt_out && (newX .*= get_Xoxmm(newXox))
    newX ./= sum(newX)
    return newX
end

# ---- Fe3+/ΣFe axis (hybrid asinh + logit) and the isotonic repair of fO2 along it
function Rlogspace(Rlo::Real, Rhi::Real, n::Integer)
    (0 < Rlo < 1 && 0 < Rhi < 1 && Rlo < Rhi) || error("Rlogspace: bad endpoints ($Rlo, $Rhi)")
    u = LinRange(log(Rlo/(1-Rlo)), log(Rhi/(1-Rhi)), n)
    return [exp(ui)/(1+exp(ui)) for ui in u]
end
function Rasinhspace(Rlo::Real, Rhi::Real, n::Integer; R0::Real=2e-3)
    (Rlo < Rhi && R0 > 0 && Rlo >= -2) || error("Rasinhspace: bad arguments ($Rlo, $Rhi, $R0)")
    u = LinRange(asinh(Rlo/R0), asinh(Rhi/R0), n)
    return [R0*sinh(ui) for ui in u]
end
function Rhybridspace(Rlo::Real, Rhi::Real, n::Integer; R0::Real=2e-3, Rknee::Real=0.4, nhi::Integer=6)
    Rhi > Rknee || return Rasinhspace(Rlo, Rhi, n; R0=R0)
    n - nhi >= 2 || error("Rhybridspace: need n - nhi >= 2")
    return vcat(Rasinhspace(Rlo, Rknee, n - nhi; R0=R0), Rlogspace(Rknee, Rhi, nhi + 1)[2:end])
end

function _isotonic(y)
    v = Float64[]; w = Int[]
    for yi in y
        push!(v, yi); push!(w, 1)
        while length(v) > 1 && v[end-1] > v[end]
            nv = (w[end-1]*v[end-1] + w[end]*v[end])/(w[end-1] + w[end]); nw = w[end-1] + w[end]
            pop!(v); pop!(w); v[end] = nv; w[end] = nw
        end
    end
    out = similar(y, Float64); k = 1
    for (vi, wi) in zip(v, w), _ in 1:wi; out[k] = vi; k += 1; end
    n = length(out)
    if n > 1
        span = max(out[end] - out[1], 1.0)
        for i in 1:n; out[i] += 1e-6*span*(i - 1)/(n - 1); end
    end
    return out
end

# fO2 array f[P, T, R, comp]: every finite column forced non-decreasing along R
function enforce_R_monotone!(f)
    fixed = 0
    for k in axes(f, 4), j in axes(f, 2), i in axes(f, 1)
        v = @view f[i, j, :, k]
        (all(isfinite, v) && !issorted(v)) || continue
        v .= _isotonic(collect(v)); fixed += 1
    end
    return fixed
end

@inline isnativeFe(ph::String) = lowercase(ph) in ("feg", "fea", "fee")
@inline function metal_wt_frac(o)
    m = 0.0
    for i in eachindex(o.ph); isnativeFe(o.ph[i]) && (m += o.ph_frac_wt[i]); end
    return m
end

# ---- Solid and melt H2O density corrections (Gerya et al. 2004 fluid; Jing & Karato hard-sphere melt)
xlogx(x) = x > 0 ? x*log(x) : zero(x)
function Gₛ(P, T, Xliq)
    H298 = -286831.56; S298 = 65.188; Vₛ = 1.71382; ϕ = 6209; c₁, c₂ = 7.23576, 0.31482
    ΔHₛ₁ = 4586.46; ΔVₛ₁ = 0.04310884; ΔH⁰ₒᵣ = -44838.80; ΔS⁰ₒᵣ = -122.397; ΔCₚ⁰ₒᵣ = 21.486
    Wh₁ = -28793.19; Ws₁ = -11.704; Wcp₁ = 5.086; T₀ = 298.15; P₀ = 1
    Ψ  = (5/4)*(P₀ + ϕ)^(1/5)*((P + ϕ)^(4/5) - (1 + ϕ)^(4/5))
    e₁ = exp(-(ΔHₛ₁ + ΔVₛ₁*Ψ)/Rgas/T); e₂ = exp(-ΔHₛ₁/Rgas/T); e₀ = exp(-ΔHₛ₁/Rgas/T₀)
    return H298 - T*S298 + Vₛ*Ψ + Rgas*T*(c₁*log(1-e₁)+c₂*log(1-e₂)) -
           (c₁ + c₂)*(ΔHₛ₁*(1-T/T₀)*e₀/(1-e₀) + Rgas*T*log(1-e₀)) +
           Rgas*T*(xlogx(1-Xliq)+xlogx(Xliq)) +
           (1-Xliq)*Rgas*T*log(ϕ*(Xliq^2)+P) -
           (1-Xliq)*(ΔH⁰ₒᵣ - T*ΔS⁰ₒᵣ + ΔCₚ⁰ₒᵣ*(T - T₀ - T*log(T/T₀))) +
           (Wh₁ - T*Ws₁ + Wcp₁*(T - T₀ - T*log(T/T₀)))*Xliq*(1-Xliq)
end
∂Gₛ∂P(P, T, Xliq)    = ForwardDiff.derivative(P -> Gₛ(P, T, Xliq), P)
∂Gₛ∂Xliq(P, T, Xliq) = ForwardDiff.derivative(Xliq -> Gₛ(P, T, Xliq), Xliq)
function Gₛ_minXliq(P, T)
    xs = LinRange(1e-4, 1-1e-4, 200); i = argmin(Gₛ.(P, T, xs))
    lo, hi = xs[max(i-1,1)], xs[min(i+1,length(xs))]
    for _ in 1:50; m = 0.5(lo+hi); (∂Gₛ∂Xliq(P, T, m) > 0) ? (hi = m) : (lo = m); end
    return 0.5(lo+hi)
end
function Gerya_H2O_volume(P_GPa, T; branch=:condensed)
    Pbar = 1e4P_GPa
    Xl = branch == :condensed ? 1.0 : Gₛ_minXliq(Pbar, T)
    return 10∂Gₛ∂P(Pbar, T, Xl)
end

function ev_P(σᵢ₀, T, Tᵣ, ηᵢ, ξᵢ, α₀, V, V₀, X)
    σᵢ   = @. σᵢ₀ * (T/Tᵣ)^(ηᵢ) * exp( ξᵢ/3 * α₀ * (T - Tᵣ) ) * (V/V₀)^(ξᵢ/3)
    Vₘᵢ  = Nₐ * π / 6 * σᵢ.^3 .* X
    Vₘᵢ₀ = Nₐ * π / 6 * σᵢ₀.^3 .* X
    Vₘ = sum(Vₘᵢ); Vₘ₀ = sum(Vₘᵢ₀); ξ = sum(ξᵢ.*Vₘᵢ)/Vₘ
    f = Vₘ/V; f₀ = Vₘ₀/V₀
    Φ = (1 + f + f^2)/(1-f)^3; Φ₀ = (1 + f₀ + f₀^2)/(1-f₀)^3
    return 1e-3Rgas*T/V * ( (1-ξ)*Φ - Φ₀*(V₀/V)^(4/3 - 1) + ξ*Φ₀*(Vₘ₀/Vₘ)^(8/3 - 1) )
end
ev_∂P∂V(σᵢ₀, T, Tᵣ, ηᵢ, ξᵢ, α₀, V, V₀, X) = ForwardDiff.derivative(V -> ev_P(σᵢ₀, T, Tᵣ, ηᵢ, ξᵢ, α₀, V, V₀, X), V)

function Jing_Karato_silicate_melt_density(P, T, Xin, Xoxin; reg=1, niter=50, α₀=2.5e-4)
    Xox = ["SiO2", "Al2O3", "FeO", "MgO", "CaO", "H2O"]
    X = zeros(length(Xox)); X .= filter_ox_list(Xin, Xoxin, Xox)
    mma = X[1]*MM.SiO2 + X[2]*MM.Al2O3 + X[3]*MM.FeO + X[4]*MM.MgO + X[5]*MM.CaO + X[6]*MM.H2O
    Tᵣ = 1673.; V₀ᵢ = SA[26.86, 37.11, 13.65, 11.69, 16.57, 22.90]; V₀ = sum(V₀ᵢ.*X)
    if reg == 1
        σᵢ₀ = 1e-7SA[0.369, 0.326, 0.286, 0.272, 0.315, 0.172]; ηᵢ = SA[-0.04, -0.03, 0.01, 0.03, -0.07, 0.06]; ξᵢ = SA[0.70, 0.58, 0.17, 0.12, 0.09, -0.35]
    else
        σᵢ₀ = 1e-7SA[0.369, 0.326, 0.285, 0.272, 0.315, 0.180]; ηᵢ = SA[-0.04, -0.03, 0.01, 0.04, -0.06, 0.13]; ξᵢ = SA[0.69, 0.58, 0.15, 0.14, 0.07, 0.0]
    end
    Pf(V) = ev_P(σᵢ₀, T, Tᵣ, ηᵢ, ξᵢ, α₀, V, V₀, X)
    hi = 1.5V₀; lo = hi; plast = -Inf
    while true
        lo = 0.99lo; p = Pf(lo)
        (isfinite(p) && p >= P) && break
        (!isfinite(p) || (plast > 0 && p < plast) || lo < 1e-2V₀) && error("Jing_Karato_silicate_melt_density: no bracket at P=$P GPa, T=$T K")
        plast = p
    end
    hi = lo/0.99
    for _ in 1:80
        m = 0.5(lo+hi); (Pf(m) > P) ? (lo = m) : (hi = m)
        (hi - lo) < 1e-9hi && break
    end
    V = 0.5(lo+hi)
    for _ in 1:niter
        dP = Pf(V) - P
        (abs(dP) < 1e-4) && break
        Vn = V - dP/ev_∂P∂V(σᵢ₀, T, Tᵣ, ηᵢ, ξᵢ, α₀, V, V₀, X)
        (isfinite(Vn) && lo < Vn < hi) || break
        V = Vn
    end
    return mma/V
end

function bm3_volume(P, V₀, K₀, K′; P₀=0.0)
    f(V) = 1.5K₀*((V₀/V)^(7/3) - (V₀/V)^(5/3))*(1 + 0.75(K′-4)*((V₀/V)^(2/3) - 1)) - (P - P₀)
    lo, hi = 0.05V₀, 1.5V₀
    for _ in 1:200; m = 0.5(lo+hi); (f(m)*f(lo) < 0) ? (hi = m) : (lo = m); end
    return 0.5(lo+hi)
end
melt_H2O_partial_volume(P, V₁) = V₁*bm3_volume(P, 25.0, 1.91, 5.80)/25.0

# Anhydrous sb24 density + expansivity and the H2O volume on a (P,T) grid, as linear interpolants
function Gerya_solid_interpolator(Xin, Xox; nP=100, nT=100, Prange=(0.1, 75.), Trange=(500., 3000.), H2O_branch=:condensed)
    X = sum(Xin) > 1.0 ? 1e-2Xin : Xin; X = X./sum(X); noWX = X[1:end-1]
    nPnT = nP*nT; P = LinRange(Prange[1], Prange[2], nP); T = LinRange(Trange[1], Trange[2], nT)
    vP, vT = repeat(P, outer=nT), repeat(T, inner=nP)
    Vh2o = zeros(nP, nT); ρanh = zeros(nPnT); α = zeros(nPnT)
    Threads.@threads for it in 1:nT
        for ip in 1:nP; Vh2o[ip, it] = Gerya_H2O_volume(P[ip], T[it]; branch=H2O_branch); end
    end
    data = Initialize_MAGEMin("sb24", verbose=false)
    out  = multi_point_minimization(10vP, vT.-273.15, data, X=[noWX for _ in 1:nPnT], Xoxides=Xox[1:end-1], sys_in="mol", name_solvus=true, progressbar=false)
    Finalize_MAGEMin(data)
    for i in 1:nPnT; ρanh[i] = out[i].rho/1e3; α[i] = out[i].alpha[1]; end
    nbad = count(r -> r < 2.4 || r > 6.0, ρanh)
    nbad > 0 && @warn "Gerya_solid_interpolator: $nbad of $nPnT MAGEMin points have rho outside 2.4-6.0 g/cm³"
    ρanh = reshape(ρanh, nP, nT); α = reshape(α, nP, nT)
    return [extrapolate(interpolate((P,T), ρanh, Gridded(Linear())), Flat()),
            extrapolate(interpolate((P,T), Vh2o, Gridded(Linear())), Flat()),
            extrapolate(interpolate((P,T), α, Gridded(Linear())), Flat())]
end

function Gerya_solid_H2O_density_correction(P, T, Xin, VVρ, MMv)
    ρA = VVρ[1](P, T); VH = VVρ[2](P, T)
    X = sum(Xin) > 1.0 ? 1e-2Xin : Xin
    XH = X./sum(X); XA = @views X[1:end-1]./sum(X[1:end-1])
    massA = sum(XA.*(@view MMv[1:end-1])); massH = sum(XH.*MMv)
    rV = (1-XH[end])*(massA/ρA) + XH[end]*VH
    return 1e2(1 - (massH/rV)/ρA)
end

# S/M density corrections [P, T, H, slot], slot 1 basalt, 2 harzburgite (reference PT_H2O_ρ, default=true)
function PT_H2O_ρ(Pi, Pf, Ti, Tf; nH=50, nP=50, nT=50, Hmin_wt=0.0, Hmax_wt=8.0, melt_H2O_model=:JK, melt_floor=true, reg=1, H2O_branch=:condensed)
    Xox = ["SiO2", "MgO", "FeO", "CaO", "Al2O3", "Na2O", "Cr2O3", "O", "H2O"]
    XB  = [50.42, 9.77, 7.1, 12.54, 16.8, 2.23, 0.07, 0.0, 0.0]
    XH  = [43.43, 45.93, 8.34, 0.9, 1.0, 0.01, 0.3, 0.0, 0.0]
    comps = (filter_ox_list(XB, Xox, Xox, wt_in=true, wt_out=false), filter_ox_list(XH, Xox, Xox, wt_in=true, wt_out=false))
    melt_H2O_model in (:JK, :ideal) || error("PT_H2O_ρ: melt_H2O_model must be :JK or :ideal")
    P = LinRange(Pi, Pf, nP); T = LinRange(Ti, Tf, nT); Hwt = collect(LinRange(Hmin_wt, Hmax_wt, nH))
    SΔρ = zeros(nP, nT, nH, 2); MΔρ = zeros(nP, nT, nH, 2); anhM = zeros(nP, nT, 2)
    MMv = get_Xoxmm(Xox)
    for (is, Xc) in enumerate(comps)
        Xl = sum(Xc) > 1.0 ? 1e-2Xc : copy(Xc); Xl = Xl./sum(Xl)
        MA = sum(Xl.*MMv)
        VVρ = Gerya_solid_interpolator(Xl, Xox; nP=nP, nT=nT, Prange=(max(Pi,1e-3), Pf), Trange=(Ti, Tf), H2O_branch=H2O_branch)
        Threads.@threads for it in 1:nT
            for ip in 1:nP; anhM[ip, it, is] = Jing_Karato_silicate_melt_density(P[ip], T[it], Xl, Xox; reg=reg); end
        end
        for iH in 1:nH
            w = Hwt[iH]
            x = (w/MMv[end])/((w/MMv[end]) + (100-w)/MA)
            X = (1-x).*Xl; X[end] = x
            MH = (1-x)*MA + x*MMv[end]
            V₁ = zeros(nT)
            if melt_H2O_model == :ideal && x > 0.0
                for it in 1:nT
                    ρd1 = Jing_Karato_silicate_melt_density(1e-4, T[it], Xl, Xox; reg=reg)
                    ρh1 = Jing_Karato_silicate_melt_density(1e-4, T[it], X,  Xox; reg=reg)
                    V₁[it] = (MH/ρh1 - (1-x)*MA/ρd1)/x
                end
            end
            x == 0.0 && continue
            Threads.@threads for it in 1:nT
                for ip in 1:nP
                    SΔρ[ip, it, iH, is] = Gerya_solid_H2O_density_correction(P[ip], T[it], X, VVρ, MMv)
                    ρd = anhM[ip, it, is]
                    ρh = melt_H2O_model == :JK ? Jing_Karato_silicate_melt_density(P[ip], T[it], X, Xox; reg=reg) :
                                                 MH/((1-x)*MA/ρd + x*melt_H2O_partial_volume(P[ip], V₁[it]))
                    d = 1e2*(1 - ρh/ρd)
                    MΔρ[ip, it, iH, is] = melt_floor ? max(d, 0.0) : d
                end
            end
        end
    end
    return SΔρ, MΔρ, collect(P), collect(T), Hwt
end

# ---- Reference table-writer number format: 4 (or 6) decimals, NaN/Inf verbatim
@inline fmtv(v, d) = isfinite(v) ? round(v, digits=d) : v
