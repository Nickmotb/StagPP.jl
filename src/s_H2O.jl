# =============================================
# ==== Pitzer and Sterner (1994) solver =======
# =============================================

    # Pitzer and Sterner (1994) coefficients and assemblers
    const cᵢⱼ = [
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
    get_cᵢ(cᵢⱼ, Tᵢ) = cᵢⱼ * [Tᵢ^(-4) Tᵢ^(-2) Tᵢ^(-1) 1.0 Tᵢ Tᵢ^2]' # cᵢ = cᵢ₁T⁻⁴ + cᵢ₂T⁻² + cᵢ₃T⁻¹ + cᵢ₄ + cᵢ₅T + cᵢ₆T²
    Pitzer_and_Sterner_EoS(cᵢ, ρ, Tᵢ) = (ρ + cᵢ[1]*ρ^2 - ρ^2*((cᵢ[3] + 2cᵢ[4]*ρ + 3cᵢ[5]*ρ^2 + 4cᵢ[6]*ρ^3) / (cᵢ[2] 
                                            + cᵢ[3]*ρ + cᵢ[4]*ρ^2 + cᵢ[5]*ρ^3 + cᵢ[6]*ρ^4)^2) + cᵢ[7]*ρ^2*exp(-cᵢ[8]*ρ) 
                                                + cᵢ[9]*ρ^2*exp(-cᵢ[10]*ρ))*(1e6R*Tᵢ)
    # Newton solver P-V-T
    function Pitzer_and_Sterner_solver(P_target, Tᵢ, cᵢⱼ; verbose=false, volume=false, ρ₀=0.1, Lbound=1e-3, Ubound=1e3, tol=1e-6, maxiter=300)
        # Initial guess
        ρ, lb, lu, ΔP, iter = ρ₀, Lbound, Ubound, Inf, 0 # mol/cm³
        cᵢ = get_cᵢ(cᵢⱼ, Tᵢ)
        ∂P∂ρ(ρ, cᵢ, Tᵢ) = ForwardDiff.derivative(ρ -> Pitzer_and_Sterner_EoS(cᵢ, ρ, Tᵢ), ρ) # ∂P/∂ρ at constant T
        while (ΔP > tol && iter < maxiter)
            # Initial guess
            P_guess = Pitzer_and_Sterner_EoS(cᵢ, ρ, Tᵢ) # Pa
            # Error assessment
            ΔP = P_guess - P_target;
            # Step estimation using derivatives
            ρn = ρ - ΔP / ∂P∂ρ(ρ, cᵢ, Tᵢ) # Newton-Raphson step
            # If exiting bounds, use bisection to fall back in domain
            (!isfinite(ρn) || ρn < lb || ρn > lu) && (ρn = 0.5(lb + lu)) # Bisection fallback
            # Update bounds depending on pressure guess
            (P_guess > P_target) ? (lu = ρ) : (lb = ρ)
            ρ = ρn; iter+=1
            (verbose) && println(" Pitzer and Sterner (1994) at $Tᵢ K and $(round(1e-5P_target, digits=4)) bar ($(round(1e-9P_target, digits=4)) GPa) \t:\t ρ = $(H₂O_mm*ρ) g/cm³ \t|\t V = $(1/ρ) cm³/mol \t (P error = $ΔP Pa) in $iter iterations")
        end
        # (verbose) && println(" Pitzer and Sterner (1994) at $Tᵢ K and $(round(1e-5P_target, digits=4)) bar ($(round(1e-9P_target, digits=4)) GPa) \t:\t ρ = $(H₂O_mm*ρ) g/cm³ \t|\t V = $(1/ρ) cm³/mol \t (P error = $ΔP Pa) in $iter iterations")
        volume ? (return 1/ρ) : (return ρ)# in cm³/mol or mol/cm³
    end

    # H₂O fugacity solver -> yields fugacity map in GPa
    function fH₂O(P::AbstractVector, T::AbstractVector)
        fmap = zeros(Float64, length(P), length(T))
        Threads.@threads for t in eachindex(T)
            Tᵢ = T[t]
            cᵢ = get_cᵢ(cᵢⱼ, Tᵢ)
            for p in eachindex(P)
                ρ = Pitzer_and_Sterner_solver(1e9min(P[p], 100.), Tᵢ, cᵢⱼ) # P -> from GPa to Pa; T in K
                # Fugacity calculation
                fmap[p,t] = exp(log(ρ) + cᵢ[1]*ρ + (1/(cᵢ[2] + cᵢ[3]*ρ + cᵢ[4]*ρ^2 + cᵢ[5]*ρ^3 + cᵢ[6]*ρ^4) - 1/cᵢ[2]) 
                                - cᵢ[7]/cᵢ[8] * (exp(-cᵢ[8]*ρ) - 1) - cᵢ[9]/cᵢ[10] * (exp(-cᵢ[10]*ρ) - 1) + 1e9min(P[p], 100.)/(ρ*1e6R*Tᵢ) + log(1e6R*Tᵢ) - 1)/1e9
            end
        end
        return fmap # Convert back to GPa
    end

    function fH₂O(P::Float64, T::Float64)
        cᵢ = get_cᵢ(cᵢⱼ, T)
        ρ = Pitzer_and_Sterner_solver(1e9min(P, 100.), T, cᵢⱼ) # P -> from GPa to Pa; T in K
        # Fugacity calculation
        return exp(log(ρ) + cᵢ[1]*ρ + (1/(cᵢ[2] + cᵢ[3]*ρ + cᵢ[4]*ρ^2 + cᵢ[5]*ρ^3 + cᵢ[6]*ρ^4) - 1/cᵢ[2])
                        - cᵢ[7]/cᵢ[8] * (exp(-cᵢ[8]*ρ) - 1) - cᵢ[9]/cᵢ[10] * (exp(-cᵢ[10]*ρ) - 1) + 1e9min(P, 100.)/(ρ*1e6R*T) + log(1e6R*T) - 1)/1e9
    end

# ========================
# ==== Integration =======
# ========================

    function ∫sᴴ²ᴼ!(fmap, Pv, Tv, min_s, outH, outB)

        # Post-stishovite boundary
        CaCl₂_st, α_PbO₂_st = CaCl₂_α_PbO₂_boundary()

        # Dry phase checker
        @inline dry_phase(ph::String) = ph ∈ ["q", "nal", "crn", "plg"]

        # Phase iterator call
        function s_phase_sum!(fmap, i, min_s, phase, ph, out, phwt, p, t, slot)
                # Transition zone database (mtl)
                (phase=="ol")      && (fmap[i, slot] += phwt * (min_s.ol(p, t) + min(out[i].SS_vec[ph].Comp[2]/0.1005, 1.0)*min_s.opx_al(p, t)))
                (phase=="wad")     && (fmap[i, slot] += phwt * min_s.wad(p, t))
                (phase=="ring")    && (fmap[i, slot] += phwt * min_s.rw(p, t))
                (phase=="opx")     && (fmap[i, slot] += phwt * min_s.opx(p, t))
                (phase=="coe")     && (fmap[i, slot] += phwt * min_s.coe(p, t))
                (phase=="crst")    && (fmap[i, slot] += phwt * min_s.crst)
                (phase=="fp")      && (fmap[i, slot] += phwt * min_s.fp)
                ((phase=="hpx"))   && (fmap[i, slot] += phwt * min_s.cpx_hp(p, t))
                (phase=="cpx")     && (fmap[i, slot] += out[i].SS_vec[ph].emFrac_wt[1] * min_s.cpx_lp_di(p, t))
                (phase=="cpx")     && (fmap[i, slot] += out[i].SS_vec[ph].emFrac_wt[4] * min_s.cpx_lp_jd(p, t))
                # Lower mantle database (Stx)
                (phase=="ak")  && (fmap[i, slot] += phwt * min_s.rw(p,t)/min_s.D_rw_aki)
                (phase=="ppv") && (fmap[i, slot] += phwt * min_s.pv * (out[i].SS_vec[ph].Comp[3]*min_s.Dppv_al(p,t) + (phwt-out[i].SS_vec[ph].Comp[3])*min_s.Dppv_noal(p,t)))
                # Found in both
                (phase=="cf")  && (fmap[i, slot] += phwt * min_s.cf)
                (phase=="cpv"   || phase=="capv")  && (fmap[i, slot] += phwt * min_s.cpv)
                (phase=="pv"    || phase=="mpv")   && (fmap[i, slot] += phwt * min_s.pv)
                (phase=="stv"   || phase=="st")    && (fmap[i, slot] += phwt * (p<=CaCl₂_st(t) ? min_s.st(p, t) : (p<=α_PbO₂_st(t) ? min_s.CaCl₂_st(p, t) : min_s.α_PbO₂_st(p, t))))
                (phase=="g"     || phase=="gtmj")  && (fmap[i, slot] += phwt * min_s.grt(p, t))
        end

        # Vectorized mesh iterator
        Threads.@threads for i in eachindex(Pv)
            # Harzburgite
            for ph in eachindex(outH[i].ph)
                # Skip if dry phases
                phase = outH[i].ph[ph]
                dry_phase(phase) && continue
                # Fraction of current phase
                phwt = outH[i].ph_frac_wt[ph]
                # Contribute to sum
                s_phase_sum!(fmap, i, min_s, phase, ph, outH, phwt, Pv[i], Tv[i], 1)
            end

            # Basalt
            for ph in eachindex(outB[i].ph)
                # Skip if dry phase
                phase = outB[i].ph[ph]
                dry_phase(phase) && continue
                # Fraction of current phase
                phwt = outB[i].ph_frac_wt[ph]
                # Contribute to sum
                s_phase_sum!(fmap, i, min_s, phase, ph, outB, phwt, Pv[i], Tv[i], 2)
            end
        end

    end

    function CaCl₂_α_PbO₂_boundary()

        # CaCl₂-type post-stishovite Transition (Umemoto 2006)
        function umemoto(x,p)
            c = zeros(length(x)); mask = x .< 1400
            @. c[mask] = p[1]*x[mask]^4 + p[2]*x[mask]^3 + p[3]*x[mask]^2 + p[4]*x[mask]
            @. c[!mask] = p[5]*x[!mask]^2 + p[6]*x[!mask] + p[7]
            return c
        end
        PSP = [19.00, 22.80, 26.60, 30.40, 34.20, 38.00, 41.80, 45.60, 49.40, 53.20, 57.00, 60.80, 64.60, 68.40]
        PST = [563.2, 706.0, 833.9, 945.3, 1055, 1141, 1212, 1282, 1347, 1400, 1476, 1544, 1652, 1894]
        fT = LinRange(300, 1894, 100)
        fit = curve_fit(umemoto, PST, PSP, ones(7))
        fP = umemoto(fT, fit.param)
        CaCl₂ = extrapolate(Interpolations.interpolate((PST,), PSP, Gridded(Linear())), Line())

        # α-PbO₂-type post-stishovite Transition (Murakami et al. 2007)
        PST = ([2000., 4050.],)
        PSP =  [110., 128.]
        α_PbO₂ = extrapolate(Interpolations.interpolate(PST, PSP, Gridded(Linear())), Line())

        return CaCl₂, α_PbO₂

    end

# =======================
# ==== Capacities =======
# =======================

    function min_sᴴ²ᴼ_assembler(ns)
        # Interpolator fields
        ol = sᴴ²ᴼ_olivine(;ns=ns)
        wad = sᴴ²ᴼ_wadsleyite(;ns=ns)
        rw = sᴴ²ᴼ_ringwoodite(;ns=ns)
        grt = sᴴ²ᴼ_garnet(;ns=ns)
        opx, opx_al = sᴴ²ᴼ_opx(;ns=ns)
        cpx_lp_di, cpx_lp_jd, cpx_hp = sᴴ²ᴼ_cpx(;ns=ns)
        st = sᴴ²ᴼ_stishovite(;ns=ns)
        CaCl₂_st = sᴴ²ᴼ_CaCl₂_stishovite(;ns=ns)
        α_PbO₂_st = sᴴ²ᴼ_α_PbO₂_stishovite(;ns=ns)
        coe = sᴴ²ᴼ_coesite(;ns=ns)
        Dppv_al, Dppv_noal = sᴴ²ᴼ_post_perovskite(;ns=ns)

        # Constants
        pv = 3e-3
        cpv = 0.0
        crst = 1e-4
        cf = pv
        fp = 5e-3
        D_rw_aki = 21.
        cor = 0.0

        return min_sᴴ²ᴼ(ol, wad, rw, grt, opx, opx_al, cpx_lp_di, cpx_lp_jd, cpx_hp, st, CaCl₂_st, α_PbO₂_st, coe, Dppv_al, 
                            Dppv_noal, pv, cpv, cf, crst, fp, D_rw_aki, cor)
    end

    function fit_Keppler(PT, W, fug, n; unit="GPa", p=ones(3))
        if unit=="Pa"
            model = function (x, p)
                P = x[:,1]; T = x[:,2]
                A, ΔH, ΔV = p
                res = @. A * fug(P,T)^n * exp(-(ΔH + ΔV*P*1e9)/(R*T))
            end
        elseif unit=="GPa"
            model = function (x, p)
                P = x[:,1]; T = x[:,2]
                A, ΔH, ΔV = p
                res = @. A * fug(P,T)^n * exp(-(ΔH + ΔV*P)/(R*T))
            end
        end
        fit = curve_fit(model, PT, W, p)
        # A, ΔH, ΔV
        return fit.param
    end

    function sᴴ²ᴼ_olivine(; ns=100)
        # Generate curve (Dong 2020)
        b = 4905.5403; n = 0.6447
        P, T = LinRange(0.01, 14.2, ns), LinRange(1200, 2400, ns)
        ol = zeros(Float64, ns, ns)
        fug = local_fug(P, T, ns)

        # Model
        for j in eachindex(T)
            for i in eachindex(P)
                ol[i,j] = 1e-4exp(0.5n*log(fug(P[i], T[j])) + b/(T[j]))
            end
        end
        average_spikes!(ol, 3.7)

        # Interpolation
        return extrapolate(Interpolations.interpolate((P,T), ol, Gridded(Linear())), Flat())
    end

    function sᴴ²ᴼ_wadsleyite(; ns=100)

        # Data Dong + Litasov
        P = [14.5, 15.2, 15.2, 15.5, 16, 16, 17, 18.5, 19, 19, 20, 20, 20, 14.0, 14.0, 14.0, 15.0, 15.5, 16.0, 18.5, 20.0, 20.0]
        T = [1873, 1673, 1873, 1673, 1773, 1773, 1673, 1673, 1723, 1923, 1773, 1873, 2173, 1473, 1573, 1673, 1573, 1673, 1773, 1573, 1773, 1873]
        W = [0.6, 2.49, 0.56, 1.06, 0.47, 1., 1.15, 1.31, 0.85, 0.35, 0.95, 0.8, 0.4, 2.07, 1.83, 0.81, 1.02, 0.58, 0.55, 1.13, 0.52, 0.44]
        
        # Computation
        fug = local_fug(P, T, ns)
        n = 1.0
        A, ΔH, ΔV = fit_Keppler([P T], W, fug, n; unit="GPa", p=[1e-8, -1e4, 1e4])

        # Generate maps
        P, T = LinRange(minimum(P), maximum(P), ns), LinRange(minimum(T), maximum(T), ns)
        wad = zeros(Float64, ns, ns)
        for j in eachindex(T)
            for i in eachindex(P)
                wad[i,j] = pmodel_GPa(P[i], T[j], A, ΔH, ΔV, n, fug)
            end
        end
        average_spikes!(wad, 10.0)

        # Unit conversion
        ΔH, ΔV = ΔH*1e-3, ΔV*1e6*1e-9

        return extrapolate(Interpolations.interpolate((P,T), wad, Gridded(Linear())), Flat())

    end

    function sᴴ²ᴼ_ringwoodite(; ns=100, verbose=false)

        # Generate curve (Dong 2020)
        n = 1.0
        P = [19., 21., 21., 23., 23., 23.]
        T = [1823., 1823., 1973., 1723., 1800., 2000.]
        W = [1.16, 1.04, 0.89, 1.60, 1.74, 0.87]
        fug = local_fug(P, T, ns)
        A, ΔH, ΔV = fit_Keppler([P T], W, fug, n; unit="GPa", p=[1e-2, -1e4, 1e4])
        
        P, T = LinRange(minimum(P), maximum(P), ns), LinRange(minimum(T), maximum(T), ns)
        rw = zeros(Float64, ns, ns)
        for j in eachindex(T)
            for i in eachindex(P)
                rw[i,j] = pmodel_GPa(P[i], T[j], A, ΔH, ΔV, n, fug)
            end
        end

        # Unit conversion
        ΔH, ΔV = ΔH*1e-3, ΔV*1e6/1e9
        verbose && println("A = $(round(A, digits=7)) ppm/bar , ΔH = $(round(ΔH, digits=3)) kJ/mol, ΔV = $(round(ΔV,digits=3)) cm³/mol, n = $n")

        return extrapolate(Interpolations.interpolate((P,T), rw, Gridded(Linear())), Flat())
    end

    function sᴴ²ᴼ_garnet(; ns=100)
        P, T = LinRange(0.01, 10.0, ns), LinRange(700., 1500., ns)
        py = zeros(Float64, ns, ns) # Lu and Keppler (1997)

        # Data
        A, ΔH, ΔV, n = 0.866, 2580, 5.71e-6, 0.5
        fug = local_fug(P, T, ns)

        # Generate map
        for j in eachindex(T)
            for i in eachindex(P)
                py[i,j] = 1e-4pmodel_GPa(P[i], T[j], A, ΔH, ΔV, n, fug)
            end
        end
        average_spikes!(py, 8.e-1)

        return extrapolate(Interpolations.interpolate((P,T), py, Gridded(Linear())), Flat())

    end

    function sᴴ²ᴼ_opx(; ns=100)
        P, T = LinRange(0.01, 14.0, ns), LinRange(700., 1500., ns)
        opx, opx_al = zeros(Float64, ns, ns), zeros(Float64, ns, ns)

        # Data
        A, ΔH, ΔV, n = 0.01354, -4563, 12.1*1e-6, 1
        A_al, ΔH_al, ΔV_al, n_al = 0.042, -79.685, 11.3*1e-6, 0.5
        fug = local_fug(P, T, ns)

        # Generate map
        for j in eachindex(T)
            for i in eachindex(P)
                opx[i,j]    = 1e-4pmodel_bar_Pa(P[i], T[j], A, ΔH, ΔV, n, fug)
                opx_al[i,j] = 1e-4pmodel_bar_Pa(P[i], T[j], A_al, ΔH_al, ΔV_al, n_al, fug)
            end
        end
        average_spikes!(opx, 2e-1)

        # Interpolations
        opx    = extrapolate(Interpolations.interpolate((P,T), opx, Gridded(Linear())), Flat())
        opx_al = extrapolate(Interpolations.interpolate((P,T), opx_al, Gridded(Linear())), Flat())

        return opx, opx_al

    end

    function sᴴ²ᴼ_cpx(; ns=100, verbose=false)

        # Low pressure cpx (Diopside + Jadeite)
            P, T, n = LinRange(1e-4, 10., ns), LinRange(873, 973, ns), 0.5
            fug = local_fug(P, T, ns)

            # Jadeite -> Bromiley and Keppler 2004
            A, ΔH, ΔV = 7.144, 0, 8.019*1e-6 # ppm/bar, J/mol, m³/mol
            # Diopside -> Bromiley et al. 2004
            A_di, ΔH_di, ΔV_di = 2.15, 0, 7.43*1e-6 # ppm/bar, J/mol, m³/mol

            # Output
            cpx_lp_di, cpx_lp_jd = zeros(Float64, ns, ns), zeros(Float64, ns, ns)
            for j in eachindex(T)
                for i in eachindex(P)
                    cpx_lp_di[i,j] = 1e-4pmodel_bar_Pa(P[i], T[j], A_di, ΔH_di, ΔV_di, n, fug)
                    cpx_lp_jd[i,j] = 1e-4pmodel_bar_Pa(P[i], T[j], A, ΔH, ΔV, n, fug)
                end
            end
            average_spikes!(cpx_lp_di, 1.5e-2)
            average_spikes!(cpx_lp_jd, 5e-2)
            cpx_lp_di = extrapolate(Interpolations.interpolate((P,T), cpx_lp_di, Gridded(Linear())), Flat())
            cpx_lp_jd = extrapolate(Interpolations.interpolate((P,T), cpx_lp_jd, Gridded(Linear())), Flat())

        # High pressure cpx (Clino-enstatite)
            n = 0.5
            P = [8., 10, 10.2, 11.6, 13.4, 13.4]
            T = [1300, 1100, 1400, 1400, 1300, 1400].+273.15;
            W = [0.12, 0.10, 0.19, 0.16, 0.49, 0.35].*1e4;
            fug = local_fug(P, T, ns)
            A, ΔH, ΔV = fit_Keppler([P T], W, fug, n; unit="Pa", p=[1e-2, -1e2, 1e-9])

            # Output
            P, T = LinRange(minimum(P), maximum(P), ns), LinRange(minimum(T), maximum(T), ns)
            cpx_hp = zeros(Float64, ns, ns)
            for j in eachindex(T)
                for i in eachindex(P)
                    cpx_hp[i,j] = 1e-4pmodel_Pa(P[i], T[j], A, ΔH, ΔV, n, fug)
                end
            end
            average_spikes!(cpx_hp, 1.5)
            cpx_hp = extrapolate(Interpolations.interpolate((P,T), cpx_hp, Gridded(Linear())), Flat())

            # unit conversion
            A = A*1e4/sqrt(1e5) # wt%/GPa⁰⁵ → ppm/bar⁰⁵
            ΔH *=1e-3; ΔV = ΔV*1e6/1e9
            verbose && println("A = $(round(A, digits=7)) ppm/bar⁰⁵ , ΔH = $(round(ΔH, digits=3)) kJ/mol, ΔV = $(round(ΔV,digits=3)) cm³/mol, n = $n")

        return cpx_lp_di, cpx_lp_jd, cpx_hp

    end

    function sᴴ²ᴼ_stishovite(; ns=100)

        # Generate curve (Chen et al. 2021)
        chen_model(p, t, A, B, ΔH, ΔV, n, Alwt) = A * fug(p,t)^n * exp(-(ΔH + ΔV*p*1e9)/(R*t)) * exp(B*Alwt/R/t) * 1e-4

        # Vectors
        P, T = LinRange(0.01, 60.01, ns), LinRange(800, 1500, ns)
        n = 0.5;        # Due to hydrogarnet substitution
        A, ΔH, ΔV, B, Alwt = 24, -3065, 4.29*1e-6, 7690, 0.0 # ppm/GPa⁰⁵, J/mol, m³/mol, J/mol, wt%
        fug = local_fug(P, T, ns)

        # Generate map
        st = zeros(Float64, ns, ns)
        for j in eachindex(T)
            for i in eachindex(P)
                st[i,j] = chen_model(P[i], T[j], A, B, ΔH, ΔV, n, Alwt)
            end
        end
        average_spikes!(st, 0.1)

        return extrapolate(Interpolations.interpolate((P,T), st, Gridded(Linear())), Flat())

    end

    function sᴴ²ᴼ_CaCl₂_stishovite(; ns=100, verbose=false)
        # --- data
        T = [2900, 2950, 3000, 3050, 3250, 3300, 3350, 3400, 3400, 3450, 3500]# K
        P = [25., 42., 41., 29, 58., 31., 51., 42., 59., 50., 100]# GPa
        W = [0.51, 0.79, 1.27, 0.82, 3.38, 0.76, 1.18, 1.51, 1.85, 1.13, 2.58] # wt%
        # avoided (3450K, 121GPa, 3.55wt%) and (3500K, 103GPa, 1.63wt%) due pressures exceeding the fugacity stability region.

        # Computation
        fug = local_fug(P, T, ns)
        n = 0.5
        A, ΔH, ΔV = fit_Keppler([P T], W, fug, n; unit="GPa")

        # Generate map
        CaCl₂_st = zeros(Float64, ns, ns)
        P, T = LinRange(minimum(P), maximum(P), ns), LinRange(minimum(T), maximum(T), ns)
        for j in eachindex(T)
            for i in eachindex(P)
                CaCl₂_st[i,j] = pmodel_GPa(P[i], T[j], A, ΔH, ΔV, n, fug)
            end
        end

         # Unit conversion
        A = A*1e4/sqrt(1e5) # wt%/GPa⁰⁵ → ppm/bar⁰⁵
        ΔH, ΔV = ΔH*1e-3, ΔV*1e6*1e-9
        verbose && println("A = $(round(A, digits=7)) ppm/bar⁰⁵ , ΔH = $(round(ΔH, digits=3)) kJ/mol, ΔV = $(round(ΔV,digits=3)) cm³/mol, n = $n")

        # Interpolation
        return extrapolate(Interpolations.interpolate((P,T), CaCl₂_st, Gridded(Linear())), Flat())
    end

    function sᴴ²ᴼ_α_PbO₂_stishovite(; ns=100, verbose=false)
        # --- data
        T = [3650., 3800., 3850., 3850., 4100.] # K
        P = [143., 144., 141., 128., 142.] # GPa
        W = [1.99, 1.97, 1.84, 1.94, 1.74] # wt%
        # avoided (3450K, 121GPa, 3.55wt%) and (3500K, 103GPa, 1.63wt%) due pressures exceeding the fugacity stability region.

        # Computation
        fug = local_fug(P, T, ns)
        n = 0.5
        A, ΔH, ΔV = fit_Keppler([P T], W, fug, n; unit="GPa", p = [1e-6, -1e4, 1e3])

        # Generate map
        α_PbO₂_st = zeros(Float64, ns, ns)
        P, T = LinRange(minimum(P), maximum(P), ns), LinRange(minimum(T), maximum(T), ns)
        for j in eachindex(T)
            for i in eachindex(P)
                α_PbO₂_st[i,j] = pmodel_GPa(P[i], T[j], A, ΔH, ΔV, n, fug)
            end
        end

         # Unit conversion
        A = A*1e4/sqrt(1e5) # wt%/GPa⁰⁵ → ppm/bar⁰⁵
        ΔH, ΔV = ΔH*1e-3, ΔV*1e6*1e-9
        verbose && println("A = $(round(A, digits=7)) ppm/bar⁰⁵ , ΔH = $(round(ΔH, digits=3)) kJ/mol, ΔV = $(round(ΔV,digits=3)) cm³/mol, n = $n")

        return extrapolate(Interpolations.interpolate((P,T), α_PbO₂_st, Gridded(Linear())), Flat())
    end

    function sᴴ²ᴼ_coesite(; ns=100)
        P, T = LinRange(1e-4, 10., ns), LinRange(800, 1073, ns)
        coe = zeros(Float64, ns, ns)
        for j in eachindex(T)
            @. coe[:,j] = max((-49 + 6.0P + 0.06T[j]), 0.0) # ppm
        end
        return extrapolate(Interpolations.interpolate((P,T), 1e-4coe, Gridded(Linear())), Flat())
    end

    function sᴴ²ᴼ_post_perovskite(; ns=100)

        function fitmodel(x, p)
            fP = x[:,1]
            fT = x[:,2]
            return @. exp( (p[1] + p[2]*fP - fT*(p[3] + p[4]*fP + p[5]*fT + p[6])) / fT )
        end

        # Vectors
        P, T = LinRange(90, 135., ns), LinRange(1000., 4000., ns)
        Dppv_al, Dppv_noal = zeros(Float64, ns, ns), zeros(Float64, ns, ns)

        # Data
        T_al = [1800, 2750, 3550, 1850, 2620, 3300, 3900, 1890, 2500, 3100, 3650, 1950, 2420, 2800, 3350, 3750, 1910, 2380, 2750, 3200, 3510, 4000]
        P_al = [90, 90, 90, 100, 100, 100, 100, 110, 110, 110, 110, 120, 120, 120, 120, 120, 130, 130, 130, 130, 130, 130]
        D_al = [0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.4, 0.1, 0.2, 0.3, 0.4, 0.1, 0.2, 0.3, 0.4, 0.5, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
        T_noal = [1080, 1200, 1310, 1670, 2820, 1075, 1130, 1230, 1380, 1700, 2900, 1030, 1100, 1170, 1240, 1420, 1770, 3000, 
                    1050, 1120, 1190, 1250, 1450, 1800, 3030, 1060, 1125, 1200, 1260, 1470, 1850, 3200,]
        P_noal = [90, 90, 90, 90, 90, 100, 100, 100, 100, 100, 100, 110, 110, 110, 110, 110, 110, 110, 120, 120, 120, 120, 120, 120, 120, 130, 130, 130, 130, 130, 130, 130]
        D_noal = [12.5, 10.0, 7.5, 5.0, 2.5, 15.0, 12.5, 10.0, 7.5, 5.0, 2.5, 17.5, 15.0, 12.5, 10.0, 7.5, 5.0, 2.5, 17.5, 15.0, 12.5, 10.0, 7.5, 5.0, 2.5, 17.5, 15.0, 12.5, 10.0, 7.5, 5.0, 2.5,]
        
        # fit
        fit_al = curve_fit(fitmodel, hcat(P_al, T_al), D_al, [1e-6; 1e-6; 1e-6; 1e-2; 1e-5; 1e-2])
        fit_noal = curve_fit(fitmodel, hcat(P_noal, T_noal), D_noal, [1e-6; 1e-6; 1e-6; 1e-2; 1e-5; 1e-2])
        p_al, p_noal = fit_al.param, fit_noal.param

        # Generate maps
        Pv = repeat(P, outer=length(T)); Tv = repeat(T, inner=length(P))
        Dppv_al .= reshape(fitmodel([Pv Tv], p_al), ns, ns)
        Dppv_noal .= reshape(fitmodel([Pv Tv], p_noal), ns, ns)

        # Interpolations
        Dppv_al = extrapolate(Interpolations.interpolate((P,T), Dppv_al, Gridded(Linear())), Flat())
        Dppv_noal = extrapolate(Interpolations.interpolate((P,T), Dppv_noal, Gridded(Linear())), Flat())
        return Dppv_al, Dppv_noal

    end

# ======================
# ======= Plots ========
# ======================

    function plot_sᴴ²ᴼ(s::sᴴ²ᴼ; cmap=:vik100, interp=false, cmap_reverse=false, logscale=true)

        # Inputs
        xlabsz, ylabsz, titlesz, xticklabsz, yticklabsz, xticksz, yticksz = 20, 20, 22, 16, 16, 12, 12

        fig = Figure(size = (1000, 1000))
        cmap_reverse && (cmap = Reverse(cmap))
        # Upper Mantle
            ax = Axis(fig[1, 1], xlabel=L"Pressure\;[\mathrm{GPa}]", ylabel=L"Temperature\;[\mathrm{K}]", title=L"Upper\;Mantle\;(Depleted, wt%)", yreversed=true,
                        xlabelsize=xlabsz, ylabelsize=ylabsz, titlesize=titlesz, xticklabelsize=xticklabsz, yticklabelsize=yticklabsz, xticksize=xticksz, yticksize=yticksz)
            hm = heatmap!(ax, s.Tum, s.Pum, logscale ? log10.(s.um[:,:,1]') : s.um[:,:,1]'; colormap=cmap, interpolate=interp); Colorbar(fig[1, 2], hm)
            ax = Axis(fig[1, 3], xlabel=L"Pressure\;[\mathrm{GPa}]", ylabel=L"Temperature\;[\mathrm{K}]", title=L"Upper\;Mantle\;(Enriched, wt%)", yreversed=true,
                        xlabelsize=xlabsz, ylabelsize=ylabsz, titlesize=titlesz, xticklabelsize=xticklabsz, yticklabelsize=yticklabsz, xticksize=xticksz, yticksize=yticksz)
            hm = heatmap!(ax, s.Tum, s.Pum, logscale ? log10.(s.um[:,:,2]') : s.um[:,:,2]'; colormap=cmap, interpolate=interp); Colorbar(fig[1, 4], hm)
        # Transition zone
            ax = Axis(fig[2, 1], xlabel=L"Pressure\;[\mathrm{GPa}]", ylabel=L"Temperature\;[\mathrm{K}]", title=L"Transition\;Zone\;(Depleted, wt%)", yreversed=true,
                        xlabelsize=xlabsz, ylabelsize=ylabsz, titlesize=titlesz, xticklabelsize=xticklabsz, yticklabelsize=yticklabsz, xticksize=xticksz, yticksize=yticksz)
            hm = heatmap!(ax, s.Ttz, s.Ptz, logscale ? log10.(s.tz[:,:,1]') : s.tz[:,:,1]'; colormap=cmap, interpolate=interp); Colorbar(fig[2, 2], hm)
            ax = Axis(fig[2, 3], xlabel=L"Pressure\;[\mathrm{GPa}]", ylabel=L"Temperature\;[\mathrm{K}]", title=L"Transition\;Zone\;(Enriched, wt%)", yreversed=true,
                        xlabelsize=xlabsz, ylabelsize=ylabsz, titlesize=titlesz, xticklabelsize=xticklabsz, yticklabelsize=yticklabsz, xticksize=xticksz, yticksize=yticksz)
            hm = heatmap!(ax, s.Ttz, s.Ptz, logscale ? log10.(s.tz[:,:,2]') : s.tz[:,:,2]'; colormap=cmap, interpolate=interp); Colorbar(fig[2, 4], hm)
        # Lower Mantle
            ax = Axis(fig[3, 1], xlabel=L"Pressure\;[\mathrm{GPa}]", ylabel=L"Temperature\;[\mathrm{K}]", title=L"Lower\;Mantle\;(Depleted, wt%)", yreversed=true,
                        xlabelsize=xlabsz, ylabelsize=ylabsz, titlesize=titlesz, xticklabelsize=xticklabsz, yticklabelsize=yticklabsz, xticksize=xticksz, yticksize=yticksz)
            hm = heatmap!(ax, s.Tlm, s.Plm, logscale ? log10.(s.lm[:,:,1]') : s.lm[:,:,1]'; colormap=cmap, interpolate=interp); Colorbar(fig[3, 2], hm)
            ax = Axis(fig[3, 3], xlabel=L"Pressure\;[\mathrm{GPa}]", ylabel=L"Temperature\;[\mathrm{K}]", title=L"Lower\;Mantle\;(Enriched, wt%)", yreversed=true,
                        xlabelsize=xlabsz, ylabelsize=ylabsz, titlesize=titlesz, xticklabelsize=xticklabsz, yticklabelsize=yticklabsz, xticksize=xticksz, yticksize=yticksz)
            hm = heatmap!(ax, s.Tlm, s.Plm, logscale ? log10.(s.lm[:,:,2]') : s.lm[:,:,2]'; colormap=cmap, interpolate=interp); Colorbar(fig[3, 4], hm)
        display(fig)
    end

# ======================
# ======= Others =======
# ======================

    function sᴴ²ᴼ_assembler!(map::Array{Float64, 2}, outH, outB, n::Int64)
        Threads.@threads for i in 1:n
            H₂O = 0.0
            # Depleted (Harzburgite)
            for j in eachindex(outH[i].SS_vec)
                ("H2O" in outH[i].SS_vec[j].emNames) && continue
                H₂O += sum(outH[i].SS_vec[j].Comp[end-1])
            end; map[i, 1] = sum(H₂O)

            # Enriched (Lherzolite)
            H₂O = 0.0
            for j in eachindex(outB[i].SS_vec)
                ("H2O" in outB[i].SS_vec[j].emNames) && continue
                H₂O += sum(outB[i].SS_vec[j].Comp[end-1])
            end; map[i, 2] = sum(H₂O)
        end
    end

    function sᴴ²ᴼ_assembler!(map::Array{Float64, 2}, out, n::Int64)
        Threads.@threads for i in 1:n
            H₂O = 0.0
            for j in eachindex(outH[i].SS_vec)
                ("H2O" in outH[i].SS_vec[j].emNames) && continue
                H₂O += sum(outH[i].SS_vec[j].Comp[end-1])
            end; map[i, 1] = sum(H₂O)
        end
    end

    function mesh_vectorization!(P, T, nP, nT, Pv, Tv)
        Pv .= repeat(P, outer=nT); Tv .= repeat(T, inner=nP)
    end

    function local_fug(P, T, ns)
        Pvec, Tvec = LinRange(minimum(P), maximum(P), ns), LinRange(minimum(T), maximum(T), ns)
        return Interpolations.interpolate((Pvec,Tvec), fH₂O(Pvec, Tvec), Gridded(Linear()))
    end

    # --- keppler NAM Functions
    pmodel_GPa(p, t, A, ΔH, ΔV, n, fug) = A * fug(p,t)^n * exp(-(ΔH + ΔV*p)/(R*t))
    pmodel_Pa(p, t, A, ΔH, ΔV, n, fug) = A * fug(p,t)^n * exp(-(ΔH + ΔV*p*1e9)/(R*t))
    pmodel_bar_Pa(p, t, A, ΔH, ΔV, n, fug) = A * (1e4fug(p,t))^n * exp(-(ΔH + ΔV*p*1e9)/(R*t))
    pmodel_bar_GPa(p, t, A, ΔH, ΔV, n, fug) = A * (1e4fug(p,t))^n * exp(-(ΔH + ΔV*p)/(R*t))

    # --- Field Smoother
    function average_spikes!(f, threshold)
        for j in axes(f, 2)
            for i in axes(f, 1)
                if f[i,j]>=threshold
                    f[i,j] = 0.0; cnt = 0
                    (i>1 && f[i-1,j] < threshold) && (f[i,j] += f[i-1,j]; cnt += 1)
                    (j>1 && f[i,j-1] < threshold) && (f[i,j] += f[i,j-1]; cnt += 1)
                    (i<size(f, 1) && f[i+1,j] < threshold) && (f[i,j] += f[i+1,j]; cnt += 1)
                    (j<size(f, 2) && f[i,j+1] < threshold) && (f[i,j] += f[i,j+1]; cnt += 1)
                    f[i,j] *= 1.0/cnt
                end
            end
        end
    end