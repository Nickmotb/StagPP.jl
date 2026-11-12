# Calibration of the high-pressure extension of the Katz et al. (2003) water solidus depression: data, model, fit, 1-sigma band.
# Included by methods_figs.jl; grid fit, profiles and AICc ranking are cached in scripts/Paper1/out/katz_fit.jls (delete it to refit).
using Printf, Serialization

# Litasov & Ohtani (2002, PEPI 134:105) Table 2 runs: (P GPa, T °C, ± °C, liquid present). K-34 lists a melt-free low-T end at 2110 °C.
LO_dry = [(10.0, 2000, 20, true), (13.5, 2000, 50, false), (13.5, 2050, 40, false), (13.5, 2100, 50, true), (13.5, 2150, 40, true),
          (16.0, 2080, 40, false), (16.5, 2120, 20, true), (17.0, 2150, 40, true), (17.0, 2220, 70, true), (19.0, 1400, 10, false),
          (19.5, 1950, 40, false), (19.5, 2110, 50, false), (19.5, 2160, 50, true), (19.5, 2180, 50, true), (19.5, 2220, 40, true),
          (21.0, 2220, 70, true), (21.5, 2200, 70, true), (22.0, 2000, 20, false), (22.0, 2100, 50, false), (22.0, 2300, 70, true),
          (24.0, 2220, 20, false), (25.0, 2350, 50, true), (25.0, 2450, 100, true)]
LO_wet = [(10.0, 1500, 10, true), (10.0, 1650, 10, true), (10.0, 1750, 10, true), (10.0, 1820, 30, true), (13.5, 1400, 10, false),
          (13.5, 1550, 10, false), (13.5, 1600, 10, true), (13.5, 1650, 10, true), (13.5, 1930, 40, true), (14.0, 1850, 40, true),
          (16.0, 1400, 10, false), (16.0, 1600, 15, false), (16.0, 2030, 50, true), (16.5, 1900, 40, true), (16.5, 1950, 50, true),
          (16.5, 2070, 30, true), (16.5, 2180, 20, true), (17.0, 1750, 20, false), (17.0, 2000, 50, true), (18.0, 1900, 10, false),
          (18.5, 1400, 10, false), (18.5, 1600, 20, false), (18.5, 2020, 50, true), (19.0, 1800, 20, false), (19.5, 2050, 20, true),
          (19.5, 2070, 60, true), (19.5, 2120, 20, true), (19.5, 2200, 50, true), (21.5, 2150, 50, true), (22.0, 2020, 20, true),
          (22.0, 2120, 50, true), (24.5, 1400, 10, false), (24.5, 2000, 20, false), (25.0, 2180, 30, true), (25.0, 2250, 70, true),
          (24.5, 2350, 100, true)]
# Solidus bounds at P, assuming it rises with pressure: above every melt-free run at P' <= P, below every melting run at P' >= P
sol_lo(R, P) = maximum((T - u for (p, T, u, m) in R if !m && p <= P); init=-Inf)
sol_hi(R, P) = minimum((T + u for (p, T, u, m) in R if m && p >= P); init=Inf)
# ΔT brackets (K) = dry minus wet; at 10 GPa no melt-free run exists: dry from the authors' 1950 °C estimate, wet above the saturated 1000 °C
lo02 = map([10.0, 13.5, 16.25, 18.5, 19.5, 22.0, 24.75]) do P
    dlo = P == 10.0 ? 1900.0 : sol_lo(LO_dry, P); wlo = P == 10.0 ? 1000.0 : sol_lo(LO_wet, P)
    (P, max(dlo - sol_hi(LO_wet, P), 0.0), sol_hi(LO_dry, P) - wlo, dlo, sol_hi(LO_dry, P), wlo, sol_hi(LO_wet, P))
end
# Authors' apparent-solidus estimates, placed on the run groups they summarise: 200 K at 16 and 250 K at 16.5 GPa (one determination,
#  abstract and section 3.2) -> 225 K at 16.25; '1900 C at 20 GPa' (180 K) from the 19.5 GPa runs; '25 GPa' (120 K) from the 24.5-25 GPa runs
lo02e = [(16.25, 225.0), (19.5, 180.0), (24.75, 120.0)]
# H2O-saturated, against the Katz et al. (2003) dry solidus (as Katz did for these data), ±50 K beyond its 8 GPa calibration
Tdry_katz(P) = 1085.7 + 132.9P - 5.1P^2
kh97 = [(P, Tdry_katz(P) - 50 - 1050, Tdry_katz(P) + 50 - 1000) for P in (10.0, 11.0)]   # Kawamoto & Holloway (1997): ~1000-1050 °C
fr17 = [(12.0, Tdry_katz(12.0) - 50 - (1500 - 273.15), Tdry_katz(12.0) + 50 - (1450 - 273.15))]  # Freitas et al. (2017): melt onset 1450-1500 K
brk = [[(x[1], x[2], x[3]) for x in lo02]; kh97; fr17]

# Model: saturated Katz depression held at Psat, times a logistic step centred on the olivine-wadsleyite transition, renormalised to 1 at 7 GPa
Lg(P, Pc, w) = 1 / (1 + exp(min((P - Pc) / w, 60.0)))
gstep(P, Pc, w, g) = P <= 7 ? 1.0 : g + (1 - g) * Lg(P, Pc, w) / Lg(7.0, Pc, w)
katz_sat(P, Psat) = (Ph = min(P, Psat); 43 * (12 * Ph^0.6 + Ph)^0.75)
model(P, θ) = katz_sat(P, θ[4]) * gstep(P, θ[1], θ[2], θ[3])               # θ = (Pc, w, ginf, Psat)

function Phi(x)   # standard normal CDF (Abramowitz-Stegun 26.2.17)
    t = 1 / (1 + 0.2316419 * abs(x)); p = t * (0.319381530 + t * (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))))
    q = exp(-x^2 / 2) / sqrt(2π) * p
    return x >= 0 ? 1 - q : q
end
# brackets: probability that the model lies inside, blurred by a 25 K model error; authors' estimates: Gaussian, 50 K
nll(f) = -sum(log(max(Phi((hi - f(P)) / 25) - Phi((lo - f(P)) / 25), 1e-300)) for (P, lo, hi) in brk) + sum(((f(P) - d) / 50)^2 / 2 for (P, d) in lo02e)

# Hill damping of the by-eye calibration, ranked against the logistic step by AICc
hill(P, Ps, n, g) = P <= 7 ? 1.0 : g + (1 - g) / (1 + ((P - 7) / Ps)^n)

FIT = joinpath(@__DIR__, "out", "katz_fit.jls")
key = hash((brk, lo02e, :prof_aicc))
fitres = isfile(FIT) ? deserialize(FIT) : nothing
if fitres === nothing || fitres.key != key
    ax4 = (11.0:0.05:16.0, 0.05:0.05:2.5, 0.05:0.01:0.45, 8.0:0.25:13.0)    # Pc, w, ginf, Psat; Psat = 8 is the first slice
    V = [nll(P -> model(P, (a, b, c, d))) for a in ax4[1], b in ax4[2], c in ax4[3], d in ax4[4]]
    k = argmin(V); ok = [getindex.(ax4, Tuple(i)) for i in CartesianIndices(V) if V[i] <= V[k] + 0.5]   # profile-likelihood 1-sigma set
    Pb = range(0.01, 25.6, 1025); C = [model(P, θ) for P in Pb, θ in ok]
    axh = (3.0:0.1:12.0, 1.0:0.5:30.0, 0.05:0.01:0.45, 8.0:0.25:13.0)       # Ps, n, ginf, Psat
    Vh = [nll(P -> katz_sat(P, d) * hill(P, a, b, c)) for a in axh[1], b in axh[2], c in axh[3], d in axh[4]]
    N = length(brk) + length(lo02e); kh = argmin(Vh); kl = argmin(V[:, :, :, 1]); kh8 = argmin(Vh[:, :, :, 1])
    cmp = [(form, Psat, kp, v, 2kp + 2v + 2kp * (kp + 1) / (N - kp - 1), θ) for (form, Psat, kp, v, θ) in
           (("logistic", :free, 4, V[k], getindex.(ax4, Tuple(k))), ("Hill", :free, 4, Vh[kh], getindex.(axh, Tuple(kh))),
            ("logistic", 8.0, 3, V[kl, 1], (getindex.(ax4[1:3], Tuple(kl))..., 8.0)), ("Hill", 8.0, 3, Vh[kh8, 1], (getindex.(axh[1:3], Tuple(kh8))..., 8.0)))]
    fitres = (key=key, θ=getindex.(ax4, Tuple(k)), nll=V[k], rng=[extrema(getindex.(ok, i)) for i in 1:4], P=Pb, lo=vec(minimum(C, dims=2)), hi=vec(maximum(C, dims=2)),
              ax=ax4, prof=[vec(minimum(V, dims=Tuple(filter(!=(i), 1:4)))) .- V[k] for i in 1:4], N=N, cmp=cmp)   # prof: ΔNLL of each parameter
    mkpath(dirname(FIT)); serialize(FIT, fitres)
end
@printf("Katz fit: Pc %.2f [%.2f-%.2f]  w %.2f [%.2f-%.2f]  ginf %.2f [%.2f-%.2f]  Psat %.2f [%.2f-%.2f]  (NLL %.2f)\n",
        fitres.θ[1], fitres.rng[1]..., fitres.θ[2], fitres.rng[2]..., fitres.θ[3], fitres.rng[3]..., fitres.θ[4], fitres.rng[4]..., fitres.nll)
aicw = (a = [c[5] for c in fitres.cmp]; e = exp.(-(a .- minimum(a)) ./ 2); e ./ sum(e))   # Akaike weights
for (c, wt) in zip(fitres.cmp, aicw)
    @printf("  %-8s Psat %-4s k=%d  NLL %5.2f  AICc %6.2f  dAICc %5.2f  weight %.2f  θ = %s\n", c[1], c[2], c[3], c[4], c[5], c[5] - minimum(x[5] for x in fitres.cmp), wt, join(map(x -> @sprintf("%.3g", x), c[6]), ", "))
end
