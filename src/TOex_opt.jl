# TOex_opt.jl — the TOex engine of TOex.jl with its hot paths rewritten for speed; same physics, same answers
# to the solver tolerance. Entry point `partition_Oₑₓ_opt`, same arguments and `out` tuple as `partition_Oₑₓ`.
# Included after TOex.jl: every law, constant and inner solver not listed below is TOex.jl's own.
#
# What changed, and why (profile of the cached solve, 2026-11-19):
#  1. Solid curve: MAGEMin's points are processed ONCE per (P, T, bulk) and kept; a cell only rescales them by Ms.
#     TOex.jl rebuilt 115 bulks, the metal and QFM columns and an Interpolations object on every call (~30 %).
#  2. Solid and melt demands are inverted in closed form (piecewise-linear and log-ratio laws) instead of an
#     80-step bisection per sink and per residual (~13 % of the residual).
#  3. ZD09 fugacities on static vectors with the mixing matrices precomputed: no allocation, ε and σ built once.
#  4. The ideality loop keeps φ and the composition in static vectors: no heap work per pass.
# Outer solver (Illinois on [-120, 40]), inner Newton solves and their tolerances are untouched.

# ---- 1. solid curve: processed once, rescaled per cell ----

# MAGEMin's solid oxybarometer per unit solid mass: excess-O mass fraction `w` (ascending), fO2 (monotone,
#  repaired), metal mass fraction, the (P,T) QFM and the (negative) metal-floor value per unit mass
struct SolidCurve
    w::Vector{Float64}
    f::Vector{Float64}
    metal::Vector{Float64}
    qfm::Float64
    floor_w::Float64
    monotonic::Bool
end

function solid_curve(P, T, X, Xox; nr=100, Rmax=0.4, nrm=15, R0=2e-3, data=nothing, mag_cache=nothing)
    key  = (P, T, Tuple(X), nr, Rmax, nrm, R0)
    ckey = (key, :curve)
    mag_cache !== nothing && haskey(mag_cache, ckey) && return mag_cache[ckey]::SolidCurve
    Xoxv   = Vector{String}(Xox)
    iO     = findfirst(==("O"), Xoxv)
    Xred   = oxidized_Xox(Xoxv, "Fe_O")
    Xdummy = (X=zeros(length(Xoxv)), Xox=copy(Xoxv), mm=get_Xoxmm(Xoxv))
    oxb(r) = oxidize_bulk(X, r, Xdummy; wt_in=true, wt_out=true, frac=true, FeFormat="Fe_O", vector=true)
    rlist = Rasinhspace(-1.96, Rmax, nr + nrm; R0=R0)
    Xlist = [oxb(r) for r in rlist]
    w     = [Xlist[i][iO]*(0.5rlist[i])/(1 + 0.5rlist[i]) for i in eachindex(rlist)]
    out = mag_cache === nothing ? nothing : get(mag_cache, key, nothing)
    if out === nothing
        own = data === nothing
        own && (data = Initialize_MAGEMin("sb24", verbose=false))
        out = multi_point_minimization(10P*ones(length(rlist)), k2c(T)*ones(length(rlist)), data, X=Xlist,
                                       Xoxides=Xred, name_solvus=true, sys_in="wt", progressbar=false)
        own && Finalize_MAGEMin(data)
        mag_cache === nothing || (mag_cache[key] = out)
    end
    f     = [o.fO2 for o in out]
    metal = [sum(o.ph_frac_wt[i] for i in eachindex(o.ph) if lowercase(o.ph[i]) in ("feg", "fea", "fee"); init=0.0) for o in out]
    qv    = sort!(filter(isfinite, [o.fO2 - o.dQFM for o in out]))
    qfm   = isempty(qv) ? NaN : qv[(length(qv) + 1) ÷ 2]
    good  = isfinite.(f)
    w, f, metal = w[good], f[good], metal[good]
    mono  = issorted(f)
    mono || (f = _isotonic(f))
    iFe   = findfirst(==("Fe"), Xred)
    c = SolidCurve(w, f, metal, qfm, -oxb(-2.0)[iFe]*(mm.O/mm.Fe), mono)   # negative, as solid_fO2_interpolant's s_floor
    mag_cache === nothing || (mag_cache[ckey] = c)
    return c
end

# The cell's view of the curve: knots in grams (w·Ms), the log-ratio continuation above, flat below
struct SolidSink
    c::SolidCurve
    Ms::Float64
    s_stoich::Float64          # Fe³⁺/ΣFe = 1 ceiling [g O]
    upper::Bool                # room left above the table
    lnN::Float64               # log10(R/(1-R)) at the table's top
end
const SOLID_SLOPE = 1/0.1917   # log-linear continuation in R/(1-R), as TOex.jl

function SolidSink(c::SolidCurve, Ms, s_stoich)
    RN = c.w[end]*Ms/s_stoich
    up = RN < 1 - 1e-9
    return SolidSink(c, Ms, s_stoich, up, up ? log10(RN/(1 - RN)) : 0.0)
end

# Linear interpolation of ys on the ascending knots xs, flat outside: Interpolations' Gridded(Linear()) + Flat()
@inline function _lin(xs, ys, x)
    n = length(xs)
    x <= xs[1] && return ys[1]
    x >= xs[n] && return ys[n]
    k = searchsortedlast(xs, x)
    return ys[k] + (ys[k+1] - ys[k])*(x - xs[k])/(xs[k+1] - xs[k])
end

# Solid fO2 at an excess-O load s [g]: TOex.jl's sfO2
@inline function sfO2(ss::SolidSink, s)
    w = s/ss.Ms; c = ss.c
    w < c.w[1]   && return c.f[1]
    w <= c.w[end] && return _lin(c.w, c.f, w)
    ss.upper || return c.f[end]
    R = min(s/ss.s_stoich, 1 - 1e-12)
    return c.f[end] + SOLID_SLOPE*(log10(R/(1 - R)) - ss.lnN)
end

# The load reading fO2 = f, saturating at [lo, hi]: the limit of TOex.jl's 80-step bisection, in closed form
function solid_demand(ss::SolidSink, f, lo, hi)
    hi <= lo && return lo
    sfO2(ss, lo) >= f && return lo
    sfO2(ss, hi) <= f && return hi
    c = ss.c
    if f <= c.f[end]                                 # inside the table: invert the segment holding f
        k = searchsortedfirst(c.f, f)
        k == 1 && return clamp(c.w[1]*ss.Ms, lo, hi)
        w = c.w[k-1] + (c.w[k] - c.w[k-1])*(f - c.f[k-1])/(c.f[k] - c.f[k-1])
        return clamp(w*ss.Ms, lo, hi)
    end
    L = ss.lnN + (f - c.f[end])/SOLID_SLOPE          # above it: invert the log-ratio continuation
    return clamp(ss.s_stoich/(1 + 10.0^(-L)), lo, hi)
end

# ---- 2. melt demand in closed form ----

# melt_fO2 is (log10 ratio + K)/a with ratio = Fe³⁺/Fe²⁺; g1 = its value at ratio 1 fixes K
struct MeltSink
    nFe::Float64               # moles of Fe in the melt
    g1::Float64                # melt fO2 at Fe³⁺/Fe²⁺ = 1
    args::NTuple{6,Float64}    # T, TOₑₓ, Mf, nB_per_mass, IDV, (unused)
end
const MELT_A = 0.1917          # Hirschmann (2022) slope, as melt_fO2

function MeltSink(T, TOₑₓ, Mf, molfracXB, nB_per_mass, IDV)
    nFe = molfracXB[3]*Mf*nB_per_mass
    g1  = nFe > 0 ? melt_fO2(0.25nFe*mm.O, T, TOₑₓ, Mf, molfracXB, nB_per_mass, IDV) : 14.0
    return MeltSink(nFe, g1, (T, TOₑₓ, Mf, nB_per_mass, IDV, 0.0))
end

function melt_demand(ms::MeltSink, f, lo, hi, molfracXB)
    hi <= lo && return lo
    T, TOₑₓ, Mf, nBm, IDV, _ = ms.args
    melt_fO2(lo, T, TOₑₓ, Mf, molfracXB, nBm, IDV) >= f && return lo
    melt_fO2(hi, T, TOₑₓ, Mf, molfracXB, nBm, IDV) <= f && return hi
    ratio = 10.0^(MELT_A*(f - ms.g1))                # Fe³⁺/Fe²⁺ that reads f
    return clamp(ratio*ms.nFe/(2(1 + ratio))*mm.O, lo, hi)
end

# ---- 3. ZD09 on static vectors ----

# ε_ij and σ_ij of the mixing rules, built once (TOex.jl rebuilt them twice per call)
const ZD_E = SMatrix{7,7}([ZD_K1[i,j]*sqrt(ZD_eps[i]*ZD_eps[j]) for i in 1:7, j in 1:7])
const ZD_S = SMatrix{7,7}([ZD_K2[i,j]*(ZD_sig[i] + ZD_sig[j])/2 for i in 1:7, j in 1:7])

# (volume [cc/mol], fugacities [bar]) at composition vx: zd_fugacity without heap work
function zd_fugacity_s(vx::SVector{7,Float64}, Pbar, TK)
    εi = ZD_E*vx;  σi = ZD_S*vx
    ε  = vx'*εi;   σ  = vx'*σi
    nanf = SVector{7,Float64}(ntuple(_ -> NaN, 7))
    (ε <= 0 || σ <= 0) && return NaN, nanf
    Tm = 154.0*TK/ε
    Vr = zd_findroot(Tm, 3.06258806161*σ^3*Pbar/ε)
    isnan(Vr) && return NaN, nanf
    vol = Vr*(σ/3.691)^3                             # L/mol
    Z  = Pbar*vol/ZD_Rgas/TK
    Vm = Vr
    c  = ZD_c
    T2, T3 = Tm^2, Tm^3
    S1  = (c[1] + c[2]/T2 + c[3]/T3)/Vm
    S1 += (c[4] + c[5]/T2 + c[6]/T3)/(2Vm^2)
    S1 += (c[7] + c[8]/T2 + c[9]/T3)/(4Vm^4)
    S1 += (c[10] + c[11]/T2 + c[12]/T3)/(5Vm^5)
    S1 += c[13]/(2c[15]*T3) * (c[14] + 1 - (c[14] + 1 + c[15]/Vm^2)*exp(-c[15]/Vm^2))
    S2  = (2c[2]/T2 + 3c[3]/T3)/Vm
    S2 += (2c[5]/T2 + 3c[6]/T3)/(2Vm^2)
    S2 += (2c[8]/T2 + 3c[9]/T3)/(4Vm^4)
    S2 += (2c[11]/T2 + 3c[12]/T3)/(5Vm^5)
    S2 += 3c[13]/(2c[15]*T3) * (c[14] + 1 - (c[14] + 1 - c[15]/Vm^2)*exp(-c[15]/Vm^2))
    base = Z - 1 - log(Z) + S1
    fug  = SVector{7,Float64}(ntuple(i -> Pbar*exp(base - 2S2*(1 - εi[i]/ε) + 6(1 - Z)*(1 - σi[i]/σ))*vx[i], 7))
    return vol*1000.0, fug
end

# fluid_fugacity on static vectors: ZD09 above FL_T_HI, TOex.jl's hybrid/blend below (rare in the mantle)
function fluid_fugacity_s(vx::SVector{7,Float64}, Pbar, TK; pure=nothing)
    TK >= FL_T_HI && return zd_fugacity_s(vx, Pbar, TK)
    v, f = fluid_fugacity(Vector(vx), Pbar, TK; pure=pure)
    return v, SVector{7,Float64}(f)
end

# ---- 4. C-O-H state at one trial fO2: coh_state_u with φ and the composition in static vectors ----

function coh_state_opt(f, T, P, Mh, Mc, nOxideTotal, TOₑₓ, φ0::SVector{7,Float64}, sd::RootSeed, pt, sol;
                       nit_min=2, nit_max=30, id_tol=1e-8, wEMOG=0.05, XCO2max=0.5, melt_H2=true, melt_H2O=true,
                       h2_comp=:basalt, XH2max=0.3, Mf=0.0, h2o_law=:sunyao, fluid_Pmax=50.0, wGs=0.01)
    Pbar, fO₂ = 1e4*P, 10.0^f
    K  = pt.K
    nH, nC = 2*Mh/mm.H2O, Mc/mm.C
    fl_ok = P <= fluid_Pmax
    φ  = fl_ok ? φ0 : SVector{7,Float64}(ntuple(_ -> 1.0, 7))
    Ms, la, lC, φw0 = Float64(sol.Ms), sol.C.la, sol.C.lC, Float64(sol.φw0)
    y_ref = fl_ok ? log(Pbar/(1/φ[4] + K[3]*sqrt(fO₂)/φ[2])) : 0.0
    ylo, yhi = y_ref - 140.0, y_ref + 20.0
    nothing_there = nC <= 0 && nH <= 0
    mk(φ) = u_trial(f, T, P, Pbar, fO₂, K, φ, nC, nH, Float64(nOxideTotal), Float64(Mf), Ms, fl_ok, true,
                    la, lC, φw0, ylo, yhi, melt_H2, melt_H2O, h2_comp, XH2max, XCO2max, h2o_law, wGs; pt=pt)
    t  = mk(φ)
    n, x, y, hfloor, it, vol_fluid = 0.0, -Inf, -Inf, false, 0, NaN
    npass = nothing_there ? 0 : (fl_ok ? max(1, nit_max) : 1)
    while it < npass
        it += 1
        r = u_solve_n(t, sd)
        n, x, y, hfloor = r.n, r.xy.x, r.xy.y, r.xy.hfloor
        fl_ok || break
        sp = u_species(t, x, y)
        S  = sum(sp)
        vx = SVector{7,Float64}(sp[3], sp[1], sp[2], sp[5], sp[4], 0.0, 0.0)   # ZD09 order: CH4 H2O CO2 H2 CO O2 C2H6
        S > 0 && (vx = vx/S)
        vol_fluid, fug = fluid_fugacity_s(vx, Pbar, T; pure=pt.pure)
        dφ = 0.0
        φn = MVector(φ)
        for i in 1:7
            vx[i] > 0 || continue
            nφ = fug[i]/(vx[i]*Pbar)
            (isfinite(nφ) && nφ > 0) || continue
            dφ = max(dφ, abs(nφ - φn[i])/φn[i])
            φn[i] = it > 4 ? sqrt(φn[i]*nφ) : nφ
        end
        φ = SVector(φn)
        it >= nit_min && dφ < id_tol && break
        t = mk(φ)
    end
    # From here on, coh_state_u's bookkeeping verbatim (no heap work in it except the rare Pmax branch)
    sp  = nothing_there ? SVector(0.0, 0.0, 0.0, 0.0, 0.0) : u_species(t, x, y)
    S   = sum(sp)
    hs  = nothing_there ? (nH2 = 0.0, nHw = 0.0, nHs = 0.0, aw = 0.0, fH2 = 0.0) : u_hsinks(t, y, sp)
    XH2O, XCO2, XCH4, XCO, XH2 = S > 0 ? (sp[1]/S, sp[2]/S, sp[3]/S, sp[4]/S, sp[5]/S) : (0.0, 0.0, 0.0, 0.0, 0.0)
    aC  = isfinite(x) ? exp(x) : 0.0
    fE  = pt.fE
    aC_cap   = exp(t.x_cap)
    molH_fl  = n*(2sp[1] + 4sp[3] + 2sp[5])
    nC_fluid = n*(sp[2] + sp[3] + sp[4])
    XCO2_melt = nOxideTotal > 0 ? softcap(aC*t.Bmelt, t.XCO2cap) : 0.0
    nC_melt  = XCO2_melt*nOxideTotal/(1 - XCO2_melt)
    nC_solid = nC > 0 ? nC*exp(min(x - t.x_cap, 0.0)/t.wC) : 0.0
    x_carb   = 1/(1 + exp(-(f - fE)/wEMOG))
    nC_carb  = x_carb*nC_solid
    nC_gph   = nC_solid - nC_carb
    carbO    = 2*nC_carb*mm.O
    netO  = n > 0 ? n*(2sp[2] + sp[4] - 2sp[3] - sp[5])*mm.O : 0.0
    fₑₓ   = isfinite(netO) ? clamp(netO, -1e3, 1e3) : NaN
    O_up  = n > 0 ? n*(2sp[2] + sp[4])*mm.O : 0.0
    O_rel = n > 0 ? n*(2sp[3] + sp[5])*mm.O : 0.0
    if !fl_ok && n > 0
        vpw, _ = fluid_fugacity_s(SVector{7,Float64}(0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0), Pbar, T; pure=pt.pure)
        isfinite(vpw) && (vol_fluid = vpw)
    end
    nH2_melt, nH_H2O_melt, fH2 = hs.nH2, hs.nHw, hs.fH2
    PH2O_melt = hs.aw*Pbar
    meltH2ex  = -nH2_melt*mm.O
    XH2_melt  = (melt_H2 && nOxideTotal > 0) ? softcap(t.KH2*fH2, XH2max) : 0.0
    SH2Om, SOH = (melt_H2O && Mf > 0) ? melt_h2o_wt(P, max(T, SY_TMIN), PH2O_melt/Pbar; law=h2o_law) : (0.0, 0.0)
    nH2O_melt  = nH_H2O_melt/2
    XH2O_melt  = nOxideTotal > 0 ? nH2O_melt/(nH2O_melt + nOxideTotal) : 0.0
    nHs, a_solid = hs.nHs, hs.aw
    h_floor_melt = hfloor && Mf > 0
    if hfloor && !h_floor_melt
        nHs     = max(nH - molH_fl - 2nH2_melt - nH_H2O_melt, 0.0)
        a_solid = exp(la[end])
    end
    rH = nH > 0 ? abs(molH_fl + 2nH2_melt + nH_H2O_melt + nHs - nH)/nH : 0.0
    rC = nC > 0 ? abs(nC_melt + nC_fluid + nC_solid - nC)/nC : 0.0
    ctol = 1e-6*max(nC, eps())
    co_fluid = n > 0 && molH_fl <= 1e-9*max(nH, 1e-30) + 1e-300
    return (; aC, aC_cap, fE, vol_fluid, it_ideality = it, res_ideality = max(rC, rH),
              graphite = nC_gph > ctol, carbonate = nC_carb > ctol, valid = true,
              XH2O, XCO2, XCH4, XCO, XH2, n_fluid = n, nC_fluid,
              XCO2_melt, nC_melt, nC_solid, nC_carb, nC_gph, carbO, fₑₓ, O_up, O_rel,
              nH2_melt, molH_fl, molH_tot = nH, fH2, XH2_melt, meltH2ex,
              nH2O_melt, nH_H2O_melt, XH2O_melt, SH2Om, SOH, PH2O_melt, lam = min(1.0, S),
              nH_solid = nHs, a_solid, S_fluid = S, h_floor = hfloor, h_floor_melt,
              fluid_absent = n <= 0, fluid_restricted = !fl_ok, co_fluid,
              φ_out = φ, φ_valid = fl_ok && !nothing_there)
end

# ---- 5. the cell: residual and driver ----

# Everything a trial fO2 needs, concretely typed; φ seed and root seeds are carried across trials as in TOex.jl
mutable struct OptCell{PT,SOL}
    ss::SolidSink
    ms::MeltSink
    molfracXB::SVector{8,Float64}
    T::Float64; P::Float64; TOₑₓ::Float64; Mf::Float64; Mc::Float64; Mh::Float64
    nOxideTotal::Float64; Φ::Float64
    maxsOₑₓ::Float64; minsOₑₓ::Float64; maxmOₑₓ::Float64; lowclip::Float64
    pt::PT; sol::SOL
    φseed::SVector{7,Float64}
    rseed::RootSeed
    kw::NamedTuple{(:melt_H2, :melt_H2O, :h2_comp, :XH2max, :XCO2max, :h2o_law, :fluid_Pmax),
                   Tuple{Bool,Bool,Symbol,Float64,Float64,Symbol,Float64}}
end

function opt_derive(f, c::OptCell)
    sOₑₓ = solid_demand(c.ss, f, -c.minsOₑₓ + c.lowclip, c.maxsOₑₓ - c.lowclip)
    mOₑₓ = melt_demand(c.ms, f, c.lowclip, c.maxmOₑₓ - c.lowclip, c.molfracXB)
    k = c.kw
    cs = coh_state_opt(f, c.T, c.P, c.Mh, c.Mc, c.nOxideTotal, c.TOₑₓ, c.φseed, c.rseed, c.pt, c.sol;
                       melt_H2=k.melt_H2, melt_H2O=k.melt_H2O, h2_comp=k.h2_comp, XH2max=k.XH2max,
                       XCO2max=k.XCO2max, Mf=c.Mf, h2o_law=k.h2o_law, fluid_Pmax=k.fluid_Pmax)
    cs.φ_valid && (c.φseed = cs.φ_out)
    XCO₂    = cs.XCO2_melt
    carbonO = c.Φ*XCO₂/(1 - XCO₂)
    fₑₓ   = nancheck("coh_state.fₑₓ", cs.fₑₓ, f, c.Mh)
    carbO = nancheck("coh_state.carbO", cs.carbO, f, c.Mc)
    fₛ    = nancheck("sfO2(sOₑₓ)", sfO2(c.ss, sOₑₓ), sOₑₓ)
    ms    = c.ms
    fₗ    = nancheck("melt_fO2", melt_fO2(mOₑₓ, ms.args[1], ms.args[2], ms.args[3], c.molfracXB, ms.args[4], ms.args[5]), mOₑₓ)
    mH2ex = nancheck("coh_state.meltH2ex", cs.meltH2ex, f, c.Mh)
    return (; f, sOₑₓ, mOₑₓ, XCO₂, carbonO, carbO, fₑₓ, mH2ex, fₛ, fₗ, cs,
              total = sOₑₓ + mOₑₓ + carbonO + carbO + fₑₓ + mH2ex)
end

# Same arguments and `out` tuple as partition_Oₑₓ (no printed report: verbose=true calls TOex.jl's driver)
function partition_Oₑₓ_opt(P::Real, T::Real, p::Real, ϕ::Real, TOex::Real, TC::Real, H::Real;
                           nr=100, Rmax=0.4, nrm=15, R0=2e-3, maxiter=500, tol=1e-8, tol_floor=1e-12,
                           verbose=false, data=nothing, out=true, mag_cache=nothing,
                           melt_H2=true, melt_H2O=true, h2_comp=:basalt, XH2max=0.3, XCO2max=0.5,
                           h2o_law=:sunyao, fluid_Pmax=50.0, sH2Oa=nothing, T_quench=298.15, closure_tol=1e-6)
    verbose && return partition_Oₑₓ(P, T, p, ϕ, TOex, TC, H; nr, Rmax, nrm, R0, maxiter, tol, tol_floor, verbose,
                                    data, out, mag_cache, melt_H2, melt_H2O, h2_comp, XH2max, XCO2max, h2o_law,
                                    fluid_Pmax, sH2Oa, T_quench, closure_tol)
    T = max(Float64(T), T_quench); P = Float64(P)
    (; X, XB, Xox, mmXox) = endmember_composition(p)
    Mt = 1e3
    (; Mf, Ms, Mc, Mh, TOₑₓ, molfracXB, nB_per_mass, Φ) = mass_budget(Mt, ϕ, TC, TOex, H, XB, mmXox)
    tbl  = sH2Oa === nothing ? load_StagH2Oa() : sH2Oa
    Csat = solid_h2o_curve(tbl, P, T, p)
    pt   = cell_consts(T, P, 1e4*P, h2_comp; pure = T < FL_T_HI ? hyb_pure(1e4*P, T) : nothing)
    sol  = (; Ms, C = (la = log.(tbl.aw), lC = log.(max.(Csat, 1e-12))), φw0 = P <= fluid_Pmax ? pure_water_φ(P, T; pure=pt.pure) : 1.0)
    IDV  = nancheck("solve_∫ΔVdP", solve_∫ΔVdP(P, T), P, T)
    maxsOₑₓ = 0.5*(Ms*X[3]/mmXox[3])*mm.O
    nOxideTotal = Mf*nB_per_mass
    maxmOₑₓ = 0.5*molfracXB[3]*nOxideTotal*mm.O
    lowclip = 1e-10
    # The solid sink: absent at ϕ = 1, as solid_fO2_interpolant's placeholder
    if Ms > 0
        curve = solid_curve(P, T, X, Xox; nr, Rmax, nrm, R0, data, mag_cache)
        ss = SolidSink(curve, Ms, maxsOₑₓ);  s_floor = curve.floor_w*Ms;  qfm_out = curve.qfm
    else
        curve = SolidCurve([0.0, 1.0], [0.0, 0.0], [0.0, 0.0], NaN, 0.0, true)
        ss = SolidSink(curve, 1.0, 1.0);  s_floor = 0.0;  qfm_out = NaN
    end
    minsOₑₓ = -s_floor
    ms  = MeltSink(T, TOₑₓ, Mf, molfracXB, nB_per_mass, IDV)
    kw  = (; melt_H2, melt_H2O, h2_comp, XH2max=Float64(XH2max), XCO2max=Float64(XCO2max), h2o_law, fluid_Pmax=Float64(fluid_Pmax))
    cell = OptCell(ss, ms, SVector{8,Float64}(molfracXB), T, P, TOₑₓ, Mf, Mc, Mh, nOxideTotal, Φ,
                   maxsOₑₓ, minsOₑₓ, maxmOₑₓ, lowclip, pt, sol, SVector{7,Float64}(ntuple(_ -> 1.0, 7)), RootSeed(), kw)
    atol_g = tol*max(abs(TOₑₓ), 1e-6*max(maxsOₑₓ + s_floor, maxmOₑₓ, floatmin(Float64)))
    atol_g = max(atol_g, tol_floor*(Mc + Mh + max(maxsOₑₓ + s_floor, maxmOₑₓ)))
    res  = bisect_solve(f -> opt_derive(f, cell).total - TOₑₓ, -120.0, 40.0; maxiter=maxiter, atol=atol_g)
    st   = opt_derive(res.x, cell)
    conv = res.converged && st.cs.res_ideality <= closure_tol && !st.cs.h_floor_melt
    sOₑₓ, mOₑₓ, XCO₂ = st.sOₑₓ, st.mOₑₓ, st.XCO₂

    # Iron bookkeeping and the output tuple: partition_Oₑₓ's, line for line
    nFe_solid = Ms > 0 ? Ms*X[3]/mmXox[3] : 0.0
    nFe_melt  = molfracXB[3]*nOxideTotal
    nO_solid  = sOₑₓ/mm.O
    FeR_melt  = nFe_melt > 0 ? 2*(mOₑₓ/mm.O)/nFe_melt : 0.0
    metal_wt   = Ms > 0 ? max(0.0, _lin(curve.w, curve.metal, clamp(sOₑₓ/Ms, curve.w[1], curve.w[end]))) : 0.0
    nFe0_solid = metal_wt*Ms/mm.Fe
    FeM_solid  = nFe_solid > 0 ? nFe0_solid/nFe_solid : 0.0
    FeM_net    = nFe_solid > 0 ? max(0.0, -nO_solid/nFe_solid) : 0.0
    FeR_net    = nFe_solid > 0 ? max(0.0, 2nO_solid/nFe_solid) : 0.0
    FeM_diss   = min(FeM_net, FeM_solid)
    FeM_disp   = max(0.0, FeM_solid - FeM_diss)
    R_net      = nFe_solid > 0 ? 2nO_solid/nFe_solid : 0.0
    disp_cap   = max(0.0, (1 - max(R_net, 0.0) - FeM_diss)/3)
    fe_clipped = FeM_disp - min(FeM_disp, disp_cap) > 1e-9
    FeM_disp   = min(FeM_disp, disp_cap)
    FeM_solid  = FeM_diss + FeM_disp
    FeR_solid  = min(1.0, max(0.0, R_net + 2FeM_solid))
    eqfO2, cohsink, cs = st.f, st.fₑₓ, st.cs
    fmq, iw = reference_buffers(P, T; qfm=qfm_out)
    out || return nothing
    oS(x) = TOₑₓ != 0 ? x/TOₑₓ : 0.0
    fC(x) = Mc > 0 ? x*mm.C/Mc : 0.0
    fH(x) = cs.molH_tot > 0 ? x/cs.molH_tot : 0.0
    return (converged=conv, fO2=eqfO2, FMQ=eqfO2-fmq, IW=eqfO2-iw, aC=cs.aC, sOₑₓ=oS(sOₑₓ), mOₑₓ=oS(mOₑₓ), XCO₂=XCO₂,
            carbonO=oS(st.carbonO), carbO=oS(st.carbO), cohsink=oS(cohsink),
            fluid_O_uptake=oS(cs.O_up), fluid_O_release=oS(cs.O_rel),
            Oex_added=max(0.0, -oS(cohsink)), rock_pool=1.0-oS(cohsink),
            FeR_solid=FeR_solid, FeR_net=FeR_net, FeM_solid=FeM_solid, FeM_net=FeM_net,
            FeM_diss=FeM_diss, FeM_disp=FeM_disp, metal_wt=metal_wt, fe_clipped=fe_clipped,
            nFe_solid=nFe_solid, FeR_melt=FeR_melt,
            fC_gph=fC(cs.nC_gph), fC_carb=fC(cs.nC_carb), fC_melt=fC(cs.nC_melt), fC_fluid=fC(cs.nC_fluid),
            fH_fluid=fH(cs.molH_fl), fH_melt=fH(cs.molH_tot - cs.molH_fl - cs.nH_solid),
            fH_solid=fH(cs.nH_solid), aH2O_solid=cs.a_solid, C_solid=solid_h2o_wt(cs.a_solid, sol.C),
            C_solid_sat=Csat[1], φ_H2O_pure=sol.φw0,
            iters=res.iters, resnorm=res.resnorm, res_id=cs.res_ideality, it_id=cs.it_ideality,
            h_floor=cs.h_floor, h_floor_melt=cs.h_floor_melt,
            XH2O=cs.XH2O, XCO2=cs.XCO2, XCH4=cs.XCH4, XCO=cs.XCO, XH2=cs.XH2,
            meltH2ex=oS(st.mH2ex), nH2_melt=cs.nH2_melt, XH2_melt=cs.XH2_melt, fH2=cs.fH2,
            molH_fluid=cs.molH_fl, molH_total=cs.molH_tot,
            XH2O_melt=cs.XH2O_melt, nH2O_melt=cs.nH2O_melt, SH2O_melt=cs.SH2Om+cs.SOH,
            SH2Om=cs.SH2Om, SOH=cs.SOH, PH2O_melt=cs.PH2O_melt, fluid_absent=cs.fluid_absent,
            fluid_restricted=cs.fluid_restricted, vol_fluid=cs.vol_fluid, lam=cs.lam,
            TOₑₓ_g=TOₑₓ, sOₑₓ_g=sOₑₓ, mOₑₓ_g=mOₑₓ, carbonO_g=st.carbonO, carbO_g=st.carbO,
            cohsink_g=cohsink, meltH2ex_g=st.mH2ex, atol_g=atol_g,
            mC_gph=cs.nC_gph*mm.C, mC_carb=cs.nC_carb*mm.C, mC_melt=cs.nC_melt*mm.C, mC_fluid=cs.nC_fluid*mm.C,
            mH_fluid=cs.molH_fl*0.5*mm.H2O, mH_melt=(cs.molH_tot - cs.molH_fl - cs.nH_solid)*0.5*mm.H2O,
            mH_solid=cs.nH_solid*0.5*mm.H2O)
end
