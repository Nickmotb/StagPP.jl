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

    function ∫sᴴ²ᴼ!(fmap, Pv, Tv, min_s, outHB, DHMS, ppaths, paths, n)

        # Post-stishovite boundary
        CaCl₂_st, α_PbO₂_st = CaCl₂_α_PbO₂_boundary()

        # Dry phase checker
        @inline dry_phase(ph::String) = ph ∈ ["q", "nal", "crn", "plg"]

        # Phase iterator call
        function s_phase_sum!(fmap, i, ib, min_s, phase, ph, out, phwt, p, t, slot, nmol)

                # Transition zone database (mtl)
                (phase=="ol")      && (fmap[i, slot] += (phwt+nmol[8])/nmol[end] * min_s.ol(p, t))
                (phase=="wad")     && (fmap[i, slot] += (phwt+nmol[9])/nmol[end] * min_s.wad(p, t))
                (phase=="ring")    && (fmap[i, slot] += phwt/nmol[end] * min_s.rw(p, t))
                (phase=="opx")     && (fmap[i, slot] += (phwt+nmol[7])/nmol[end] * min_s.opx(p, t)+ min(out[ib].SS_vec[ph].Comp[2]/0.1005, 1.0)*min_s.opx_al(p, t))
                (phase=="coe")     && (fmap[i, slot] += phwt/nmol[end] * min_s.coe(p, t))
                (phase=="crst")    && (fmap[i, slot] += phwt/nmol[end] * min_s.crst)
                (phase=="fp")      && (fmap[i, slot] += phwt/nmol[end] * min_s.fp)
                ((phase=="hpx"))   && (fmap[i, slot] += phwt/nmol[end] * min_s.cpx_hp(p, t))
                (phase=="cpx")     && (fmap[i, slot] += phwt/nmol[end] * out[ib].SS_vec[ph].emFrac[1] * min_s.cpx_lp_di(p, t))
                (phase=="cpx")     && (fmap[i, slot] += phwt/nmol[end] * out[ib].SS_vec[ph].emFrac[4] * min_s.cpx_lp_jd(p, t))
                # Lower mantle database (Stx)
                (phase=="ak")  && (fmap[i, slot] += phwt/nmol[end] * min_s.rw(p,t)/min_s.D_rw_aki)
                (phase=="ppv") && (fmap[i, slot] += phwt/nmol[end] * min_s.pv * (out[ib].SS_vec[ph].Comp[3]*min_s.Dppv_al(p,t) + (phwt-out[ib].SS_vec[ph].Comp[3])*min_s.Dppv_noal(p,t)))
                # Found in both
                (phase=="cf")  && (fmap[i, slot] += phwt/nmol[end] * min_s.cf)
                (phase=="cpv"   || phase=="capv")  && (fmap[i, slot] += phwt/nmol[end] * min_s.cpv)
                (phase=="pv"    || phase=="mpv")   && (fmap[i, slot] += (phwt+nmol[11])/nmol[end] * min_s.pv)
                (phase=="stv"   || phase=="st")    && (fmap[i, slot] += (phwt+nmol[10])/nmol[end] * (p<=CaCl₂_st(t) ? min_s.st(p, t) : (p<=α_PbO₂_st(t) ? min_s.CaCl₂_st(p, t) : min_s.α_PbO₂_st(p, t))))
                (phase=="g"     || phase=="gtmj")  && (fmap[i, slot] += phwt/nmol[end] * min_s.grt(p, t))

        end

        # Vectorized mesh iterator
        Threads.@threads for i in 1:n

            # DHMS correction
            nmol = @MVector zeros(Float64, 12); nmol[end]=1.0

            # Basalt
            for ph in eachindex(outHB[i+n].ph)
                ib = i+n
                # Skip if dry phase
                phase = outHB[ib].ph[ph]
                dry_phase(phase) && continue
                # Fraction of current phase
                phwt = outHB[ib].ph_frac[ph]
                # Contribute to sum
                s_phase_sum!(fmap, i, ib, min_s, phase, ph, outHB, phwt, Pv[i], Tv[i], 2, nmol)
            end

            # Harzburgite
            DHMS && (nmol .= DHMS_solve(outHB[i], ppaths, paths))
            for ph in eachindex(outHB[i].ph)
                # Skip if dry phases
                phase = outHB[i].ph[ph]
                dry_phase(phase) && continue
                # Fraction of current phase
                phwt = outHB[i].ph_frac[ph]
                # Contribute to sum
                s_phase_sum!(fmap, i, i, min_s, phase, ph, outHB, phwt, Pv[i], Tv[i], 1, nmol)
            end
            # Dense hydrous magnesium silicates
            DHMS && (fmap[i, 1] += (nmol[1:5]./nmol[end] ⋅ [min_s.PhA, min_s.PhE, min_s.shB, min_s.PhD, min_s.PhH]))
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

# ==============================================
# ==== Dense Hydrous Magnesium Silicates =======
# ==============================================

    function build_reactions()
        # r1
        P = [0.1923, 0.4103, 0.7821, 1.167, 1.590, 1.974, 2.410, 2.872, 3.218, 3.551, 3.936, 4.218, 4.513, 4.795, 5.038, 5.282, 5.385, 5.436, 5.462]
        T = [886.8, 937.6, 980.3, 1001.0, 1019.0, 1025.0, 1025.0, 1013.0, 1001.0, 974.2, 931.5, 899.0, 854.2, 801.4, 746.4, 699.7, 659.0, 600.0, 543.1]
        r1 = extrapolate(Interpolations.interpolate((P,), T, Gridded(Linear())), Line())
        # r2
        P = [12.70, 12.80, 12.90, 13.00, 13.10, 13.20, 13.30, 13.40, 13.50, 13.70, 13.705]
        T = [1276., 1265., 1246., 1231., 1215., 1199., 1183., 1170., 1151., 1119., 1100.]
        r2 = extrapolate(Interpolations.interpolate((P,), T, Gridded(Linear())), Line())
        # r3
        P = [12.00, 12.20, 12.40, 12.60, 12.80, 13.00, 13.20, 13.40]
        T = [805.0, 838.6, 866.7, 894.3, 922.5, 950.1, 977.6, 1005.0]
        r3 = extrapolate(Interpolations.interpolate((P,), T, Gridded(Linear())), Line())
        # r4
        P = [13.90, 14.20, 14.50, 14.80, 15.10, 15.40, 15.70, 16.00, 16.30, 16.60, 16.90, 17.20]
        T = [1070., 1110., 1143., 1171., 1197., 1222., 1249., 1281., 1303., 1331., 1357., 1383.]
        r4 = extrapolate(Interpolations.interpolate((P,), T, Gridded(Linear())), Line())
        # r5
        P = [25.60, 25.80, 26.00, 26.20, 26.40, 26.60, 26.80, 27.00]
        T = [1435., 1396., 1363., 1325., 1296., 1252., 1220., 1183.]
        r5 = extrapolate(Interpolations.interpolate((P,), T, Gridded(Linear())), Line())
        # r6
        P = [36.90, 37.80, 38.70, 39.60, 40.50, 41.40, 42.30, 43.20, 44.10, 45.00, 45.90, 46.80]
        T = [1369., 1376., 1379., 1387., 1394., 1396., 1391., 1383., 1364., 1320., 1190., 1038.]
        r6 = extrapolate(Interpolations.interpolate((P,), T, Gridded(Linear())), Line())
        # exit
        P = [0.1923, 0.4103, 0.7821, 1.167, 1.590, 1.974, 2.410, 2.872, 3.218, 3.551, 3.936, 4.218, 4.513, 4.795, 5.038, 5.200, 5.600, 6.000, 6.400, 6.800, 7.200, 7.600, 8.000, 8.400, 8.800, 9.200, 9.600, 10.00, 
                10.40, 10.80, 11.20, 11.60, 12.00, 12.40, 12.80, 13.20, 13.60, 14.00, 14.40, 14.80, 15.20, 
                    15.60, 16.00, 16.40, 16.80, 17.20, 17.60, 18.00, 18.40, 18.80, 19.20, 19.60, 20.00, 20.40, 
                        20.80, 21.20, 21.60, 22.00, 22.40, 22.80, 23.20, 23.60, 24.00, 24.40, 24.80, 25.20, 25.60, 
                            26.00, 26.40, 26.80, 27.20, 27.60, 28.00, 28.40, 35.56, 36.40, 
                                37.10, 37.80, 38.50, 39.20, 39.90, 40.60, 41.30, 42.00, 42.70, 43.40, 44.10, 44.80, 45.50, 
                                    46.20, 46.90, 47.60, 48.30, 49.00, 49.70, 50.40, 51.10, 51.80, 52.50, 53.20, 53.90, 54.60, 
                                        55.30, 56.00, 56.70, 57.40, 58.10, 58.80, 59.50, 59.90, 60.00, 60.21, 60.42, 60.52, 60.84, 
                                            60.94, 61.26, 61.57, 61.78, 62.10, 62.31, 62.62, 62.83]
        T = [886.8, 937.6, 980.3, 1001.0, 1019.0, 1025.0, 1025.0, 1013.0, 1001.0, 974.2, 931.5, 899.0, 854.2, 801.4, 746.4, 816.9, 850.2, 881.6, 913.4, 944.6, 971.2, 1009.0, 1039.0, 1073.0, 1107.0, 1144.0, 1157.0, 1183.0, 
            1205.0, 1222.0, 1239.0, 1254.0, 1267.0, 1273.0, 1297.0, 1397.0, 1480.0, 1474.0, 1467.0, 1459.0, 1452.0, 
                1444.0, 1434.0, 1425.0, 1416.0, 1407.0, 1407.0, 1417.0, 1427.0, 1437.0, 1448.0, 1458.0, 1468.0, 1477.0, 
                    1484.0, 1490.0, 1500.0, 1502.0, 1510.0, 1518.0, 1525.0, 1524.0, 1516.0, 1507.0, 1494.0, 1478.0, 1462.0, 
                        1452.0, 1440.0, 1429.0, 1419.0, 1407.0, 1395.0, 1385.0, 1363.0, 1370.0, 
                            1384.0, 1398.0, 1413.0, 1427.0, 1441.0, 1455.0, 1469.0, 1483.0, 1496.0, 1510.0, 1524.0, 1539.0, 1553.0, 
                                1566.0, 1581.0, 1594.0, 1607.0, 1615.0, 1627.0, 1640.0, 1651.0, 1658.0, 1663.0, 1669.0, 1674.0, 1679.0, 
                                    1676.0, 1673.0, 1669.0, 1659.0, 1650.0, 1628.0, 1587.0, 1553.0, 1518.0, 1481.0, 1433.0, 1390.0, 1333.0, 
                                        1280.0, 1214.0, 1135.0, 1056.0, 984.2, 905.6, 839.7, 783.1]
        e = extrapolate(Interpolations.interpolate((P,), T, Gridded(Linear())), Line())
        return DHMS_boundaries(r1, r2, r3, r4, r5, r6, e)
    end

    function compute_path(P, s; 
                            T1=700.0, P1=5.0,
                            cold=(g1=10.0, g∞=5.0,  L=1.0, m=1.1),
                            warm=(g1=100.0, g∞=40.0, L=4.0, m=1.1),
                            bump=(A=20.0, Pb=5.0, w=4.0))

        g1 = (1-s)*cold.g1 + s*warm.g1
        g∞ = (1-s)*cold.g∞ + s*warm.g∞
        L  = (1-s)*cold.L  + s*warm.L
        m  = (1-s)*cold.m  + s*warm.m
        A, Pb, w = bump.A, bump.Pb, bump.w
        @inline ∂T∂P(P) = @. g∞ + (g1 - g∞) / (1 + ((P - P1)/L)^(m)) + A * exp(-0.5 * ((P - Pb)/w)^2)

        path = ∂T∂P(P) * step(P); path[1] += T1
        path .= cumsum(path)
        return path
    end

    function test_paths(paths, reactions)
        P = LinRange(5.0, 65.0, 500)
        fig = Figure(size=(800,600))
        ax = Axis(fig[1,1])
        for i in eachindex(paths)
            plot!(ax, P, paths[i].T(P), color = paths[i].Bᵢ==0.0 ? :gray : :red)
        end
        lines!(ax, P[1:20], reactions.r1(P[1:20]), color=:red, linewidth=2)
        lines!(ax, P[1:20], reactions.r2(P[1:20]), color=:green, linewidth=2)
        lines!(ax, P[1:20], reactions.r3(P[1:20]), color=:orange, linewidth=2)
        lines!(ax, P, reactions.r4(P), color=:purple, linewidth=2)
        lines!(ax, P, reactions.r5(P), color=:brown, linewidth=2)
        lines!(ax, P, reactions.e(P), color=:blue, linewidth=2)
        ylims!(ax, 0.0, 2500.)
        # display(fig)
        CairoMakie.save("path_test.png", fig)
    end

    # Assumes DHMS are entirely isolated within region of stability. And only DHMS
    # phase carries over across boundaries.

    function path_solve(XH, Clist, phase_out, Pvtz, Tvtz, DBswitchP, outH, test_path; npaths=50, ns=100, Pend=130.0)

        # Initialize variables
        P = LinRange(5.0, Pend, ns)
        blends = LinRange(0.0, 1.0, npaths)
        path_collection = PTpath[]
        reactions = build_reactions()
        sidx, eidx = findfirst(P.>=DBswitchP), findlast(P.<=25.0) # P range in which to look for seeds

        # Initialise MAGEMin_C
        data    = Initialize_MAGEMin("um", verbose=false, buffer="aH2O");

        # Generate paths
        for s in blends
            # Temperature path
            Tpath = interpolate((P,), compute_path(P, s), Gridded(Linear()))
            # Reaction history vector
            rh = fill("", ns)
            # Find first crossing
            idx = max(findfirst(Tpath(P) .>= reactions.r1(P))-1, 1)
            # Minimize crossing point
            rm_list = remove_phases(phase_out, "um")
            out = single_point_minimization(10*P[idx], Tpath(P[idx])-273.15, data, X=XH, Xoxides=Clist, B=1.0, sys_in="wt", name_solvus=true, rm_list=rm_list)
            # Extract atg mol%, if none get out.
            atg = 0.0
            any(out.ph .== "atg") ? (atg = out.ph_frac[findfirst(out.ph .== "atg")]; rh .*= "1") : (push!(path_collection, PTpath(Tpath, atg, rh, 0.0, 0.0, 0.0, 0.0)); continue)
            # Initialise seeds
            en_seed, fo_seed, wad_seed, st_seed = 0.0, 0.0, 0.0, 0.0

            # Now iterate on path
            for i in idx+1:ns
                # Assess closest minimized point and track seed
                if i >= sidx && i <= eidx
                    min_idx = argmin((Pvtz .- P[i]).^2 .+ (Tvtz .- Tpath(P[i])).^2)
                    ph2em_map = Dict(ph => n for (n, ph) in enumerate(outH[min_idx].ph))
                    ("ol" in outH[min_idx].ph) && (fo_seed = max(fo_seed, outH[min_idx].SS_vec[ph2em_map["ol"]].emFrac[1]))
                    ("opx" in outH[min_idx].ph) && (en_seed = max(en_seed, outH[min_idx].SS_vec[ph2em_map["opx"]].emFrac[1]))
                    ("wad" in outH[min_idx].ph) && (wad_seed = max(wad_seed, outH[min_idx].ph_frac[ph2em_map["wad"]]))
                    ("stv" in outH[min_idx].ph) && (st_seed = max(st_seed, outH[min_idx].ph_frac[ph2em_map["stv"]]))
                end

                # Current P-T
                p, t = P[i], Tpath(P[i])
                # Check boundaries and update reaction history
                (t >= reactions.e(p)) && (rh[i] = "e"; rh[i+1:end] .= ""; break)
                # PhA → PhE region
                (t >= reactions.r2(p) && t>= reactions.r4(p)) && (occursin("2", rh[i]) ? continue : (rh[i:end].*="2"; continue))
                # PhE → ShB region
                ((t <= reactions.r4(p) && t <= reactions.r5(p)) && occursin("2", rh[i-1])) && (occursin("4", rh[i]) ? continue : (rh[i:end].*="4"; continue))
                # PhA → ShB region
                ((t <= reactions.r3(p) && t <= reactions.r5(p)) && !occursin("2", rh[i-1]) && (occursin("3", rh[i]) ? continue : (rh[i:end].*="3"; continue)))
                # ShB → PhD region
                (t >= reactions.r5(p) && t <= reactions.r6(p)) && (occursin("5", rh[i]) ? continue : (rh[i:end].*="5"; continue))
                # PhD → PhH region
                (t >= reactions.r6(p)) && (occursin("6", rh[i]) ? continue : (rh[i:end].*="6"; continue))
            end

            # Ensure non-zero seeds
            en_seed==0.0 && (en_seed = 0.1)
            fo_seed==0.0 && (fo_seed = 0.1)
            wad_seed==0.0 && (wad_seed = 0.1)
            st_seed==0.0 && (st_seed = 0.1)
            
            # Store path
            push!(path_collection, PTpath(Tpath, atg, rh, en_seed, fo_seed, wad_seed, st_seed))
        end
        Finalize_MAGEMin(data);
        test_path && test_paths(path_collection, reactions)
        return P, path_collection
    end

    function DHMS_solve(out, ppaths, paths)
        # P and T
        P = 1e-1out.P_kbar; T = out.T_C+273.15
        # Assign closest path (requires grid search)
        # Δmap = zeros(Float64, length(paths), length(ppaths)) # Empty distance map
        # [Δmap[i,:] .= (paths[i].T(ppaths) .- T).^2 for i in eachindex(paths)] # Fill with Temperature Δ
        # [Δmap[i,:] .+= (ppaths .- P).^2 for i in eachindex(paths)] # Fill with Pressure Δ
        # closest_point = argmin(Δmap)
        # path_idx = closest_point[1]
        # point_rh = paths[path_idx].rh[closest_point[2]] # Path history
        path_idx = argmin([abs(paths[i].T(P) - T) for i in eachindex(paths)])
        point_rh = paths[path_idx].rh[argmin(abs.(ppaths.-P))]
        # Reaction/Output vector
        # Order: (PhA, PhE, shB, PhD, PhH, atg, en, fo, wad, st, pv, normalization)
        nmol = @MVector zeros(Float64, 12); nmol[end]=1.0
        (point_rh=="" || point_rh=="e") && (return nmol) # No DHMS present
        nmol[6] = paths[path_idx].Bᵢ
        # Assign seeds to static vector
        for (i, ph) in enumerate(out.ph)
            nmol[7] += ph=="opx" ? out.SS_vec[i].emFrac[1] : paths[path_idx].en_seed
            nmol[8] += ph=="ol"  ? out.SS_vec[i].emFrac[1] : paths[path_idx].fo_seed
            nmol[9] += ph=="wad" ? out.ph_frac[i] : paths[path_idx].wad_seed
            nmol[10] += (ph=="stv" || ph=="st") ? out.ph_frac[i] : paths[path_idx].st_seed
            occursin("pv", ph) && (nmol[11] += out.ph_frac[i])
        end
        nmol2 = copy(nmol) # stores present phases but doesn't change. Nmol only tracks transitions

        # Compute mass balance
        if occursin("1", point_rh)
            nmol[1] += 14/5*nmol[6] # 5 atg → 14 PhA
            nmol[7] += 142/5*nmol[6] # 5 atg → 142 en
            nmol[6] = 0.0
        end
        if occursin("2", point_rh)
            # Ensure only DHMS crosses boundary
            nmol[7] = nmol2[7] # en
            # Assess limiting term
            lim = (nmol[1]/5790)<(nmol[7]/21831) ? 1 : 2
            # Form products
            nmol[9] += 18757/(lim==1 ? 5790 : 21831)*nmol[lim==1 ? 1 : 7] # wad
            nmol[2] += 11630/(lim==1 ? 5790 : 21831)*nmol[lim==1 ? 1 : 7] # PhE
            # Remove reactants
            nmol[7] -= lim==1 ? 21831/5790*nmol[1] : 0.0; nmol[1] = 0.0 # en + PhA
        end
        if occursin("3", point_rh)
            # Ensure only DHMS crosses boundary
            nmol[7] = nmol2[7] # en
            # Assess limiting term
            lim = (nmol[1]/2)<(nmol[8]/11) ? 1 : 2
            # Form products
            nmol[3] += 3/(lim==1 ? 2 : 11)*nmol[lim==1 ? 1 : 8] # shB
            nmol[7] += 6/(lim==1 ? 2 : 11)*nmol[lim==1 ? 1 : 8] # en
            # Remove reactants
            nmol[8] -= lim==1 ? 11/2*nmol[1] : 0.0; nmol[1] = 0.0 # ol + PhA
        end
        if occursin("4", point_rh)
            # Ensure only DHMS crosses boundary
            nmol[9] = nmol2[9] # wad
            # Assess limiting term
            lim = (nmol[2]/19265)<(nmol[9]/62762) ? 1 : 2
            # Form products
            nmol[3]  += 16288/(lim==1 ? 19265 : 62762)*nmol[lim==1 ? 2 : 9] # shB
            nmol[10] += 38172/(lim==1 ? 19265 : 62762)*nmol[lim==1 ? 2 : 9] # st
            # Remove reactants
            nmol[9] -= lim==1 ? 62762/19265*nmol[2] : 0.0; nmol[2] = 0.0 # wad + PhE
        end
        if occursin("5", point_rh)
            # Ensure only DHMS crosses boundary
            nmol[7] = nmol2[7]
            nmol[10] = nmol2[10]
            # Assess limiting term
            lim = (nmol[3]/1017)<(nmol[10]/7839) ? 1 : 2
            # Form products
            nmol[4] += 1800/(lim==1 ? 1017 : 7839)*nmol[lim==1 ? 3 : 10] # PhD
            nmol[11] += 8046/(lim==1 ? 1017 : 7839)*nmol[lim==1 ? 3 : 10] # pv
            # Remove reactants
            nmol[10] -= lim==1 ? 7839/1017*nmol[3] : 0.0; nmol[3] = 0.0 # st + shB
        end
        if occursin("6", point_rh)
            # Ensure only DHMS crosses boundary
            nmol[11] = nmol2[11]
            # Form products
            nmol[5]  += nmol[4] # PhH
            nmol[10] += nmol[4] # st
            # Remove reactants
            nmol[4] = 0.0
        end
        # Remove already present phases base (because sumed during contribution)
        nmol[6:end-1] .= max.(nmol[6:end-1] .- nmol2[6:end-1], 0.0)
        # Total amount of moles to normalize on during contribution
        nmol[end] = 1.0+sum(nmol[1:end-1])
        return nmol
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

        # DHMS (temporary)
        PhA = 11.3 # Maurice et al. 2018
        PhE = 11.9 # Maurice et al. 2018
        shB = 0.5(11.1 + 11.1/1.9) # Kakizawa et al. 2018
        PhD = 0.5(6.7 + 11.2) # Bolfan Casanova et al. 2018
        PhH = 15.2 # Nishiet et al. 2014

        return min_sᴴ²ᴼ(ol, wad, rw, grt, opx, opx_al, cpx_lp_di, cpx_lp_jd, cpx_hp, st, CaCl₂_st, α_PbO₂_st, coe, Dppv_al, 
                            Dppv_noal, pv, cpv, cf, crst, fp, D_rw_aki, cor, PhA, PhE, shB, PhD, PhH)
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

    """
        Plot sᴴ²ᴼ results.

        \t Basic usage: \t `plot_sᴴ²ᴼ(s; kwargs...)`

        Optional arguments (kwargs):
            - `cmap::Symbol`: Colormap to use. Default is `:vik100`.
            - `interp::Bool`: Whether to interpolate the heatmap. Default is `false`.
            - `cmap_reverse::Bool`: Whether to reverse the colormap. Default is `false`.
            - `logscale::Bool`: Whether to use log scale for color mapping. Default is `true`.
            - `savein::String`: Path to save the figure. Default is `""` (does not save).
            - `bigpicture::Tuple{Bool, Int}`: Whether to create a big picture figure. Default is `(false, 1)`.
    """
    function plot_sᴴ²ᴼ(s; cmap=:vik100, interp=false, cmap_reverse=false, logscale=true, savein="", bigpicture=(false, 1) )

        # Inputs
        xlabsz, ylabsz, titlesz, xticklabsz, yticklabsz, xticksz, yticksz = 20, 20, 22, 16, 16, 12, 12

        fig = Figure(size = (1000, 1000))
        cmap_reverse && (cmap = Reverse(cmap))
        # Upper Mantle
            ax = Axis(fig[1, 1], ylabel=L"Pressure\;[\mathrm{GPa}]", xlabel=L"Temperature\;[\mathrm{K}]", title=L"Upper\;Mantle\;(Depleted, wt%)", yreversed=true,
                        xlabelsize=xlabsz, ylabelsize=ylabsz, titlesize=titlesz, xticklabelsize=xticklabsz, yticklabelsize=yticklabsz, xticksize=xticksz, yticksize=yticksz)
            hm = heatmap!(ax, s.Tum, s.Pum, logscale ? log10.(s.um[:,:,1]') : s.um[:,:,1]'; colormap=cmap, interpolate=interp, colorrange= logscale ? (-2., log10(maximum(s.um[:,:,1]))) : (minimum(s.um[:,:,1]), maximum(s.um[:,:,1]))); Colorbar(fig[1, 2], hm)
            ax = Axis(fig[1, 3], ylabel=L"Pressure\;[\mathrm{GPa}]", xlabel=L"Temperature\;[\mathrm{K}]", title=L"Upper\;Mantle\;(Enriched, wt%)", yreversed=true,
                        xlabelsize=xlabsz, ylabelsize=ylabsz, titlesize=titlesz, xticklabelsize=xticklabsz, yticklabelsize=yticklabsz, xticksize=xticksz, yticksize=yticksz)
            hm = heatmap!(ax, s.Tum, s.Pum, logscale ? log10.(s.um[:,:,2]') : s.um[:,:,2]'; colormap=cmap, interpolate=interp, colorrange= logscale ? (-2., log10(maximum(s.um[:,:,2]))) : (minimum(s.um[:,:,2]), maximum(s.um[:,:,2]))); Colorbar(fig[1, 4], hm)
        # Transition zone
            ax = Axis(fig[2, 1], ylabel=L"Pressure\;[\mathrm{GPa}]", xlabel=L"Temperature\;[\mathrm{K}]", title=L"Transition\;Zone\;(Depleted, wt%)", yreversed=true,
                        xlabelsize=xlabsz, ylabelsize=ylabsz, titlesize=titlesz, xticklabelsize=xticklabsz, yticklabelsize=yticklabsz, xticksize=xticksz, yticksize=yticksz)
            hm = heatmap!(ax, s.Ttz, s.Ptz, logscale ? log10.(s.tz[:,:,1]') : s.tz[:,:,1]'; colormap=cmap, interpolate=interp); Colorbar(fig[2, 2], hm)
            ax = Axis(fig[2, 3], ylabel=L"Pressure\;[\mathrm{GPa}]", xlabel=L"Temperature\;[\mathrm{K}]", title=L"Transition\;Zone\;(Enriched, wt%)", yreversed=true,
                        xlabelsize=xlabsz, ylabelsize=ylabsz, titlesize=titlesz, xticklabelsize=xticklabsz, yticklabelsize=yticklabsz, xticksize=xticksz, yticksize=yticksz)
            hm = heatmap!(ax, s.Ttz, s.Ptz, logscale ? log10.(s.tz[:,:,2]') : s.tz[:,:,2]'; colormap=cmap, interpolate=interp); Colorbar(fig[2, 4], hm)
        # Lower Mantle
            ax = Axis(fig[3, 1], ylabel=L"Pressure\;[\mathrm{GPa}]", xlabel=L"Temperature\;[\mathrm{K}]", title=L"Lower\;Mantle\;(Depleted, wt%)", yreversed=true,
                        xlabelsize=xlabsz, ylabelsize=ylabsz, titlesize=titlesz, xticklabelsize=xticklabsz, yticklabelsize=yticklabsz, xticksize=xticksz, yticksize=yticksz)
            hm = heatmap!(ax, s.Tlm, s.Plm, logscale ? log10.(s.lm[:,:,1]') : s.lm[:,:,1]'; colormap=cmap, interpolate=interp); Colorbar(fig[3, 2], hm)
            ax = Axis(fig[3, 3], ylabel=L"Pressure\;[\mathrm{GPa}]", xlabel=L"Temperature\;[\mathrm{K}]", title=L"Lower\;Mantle\;(Enriched, wt%)", yreversed=true,
                        xlabelsize=xlabsz, ylabelsize=ylabsz, titlesize=titlesz, xticklabelsize=xticklabsz, yticklabelsize=yticklabsz, xticksize=xticksz, yticksize=yticksz)
            hm = heatmap!(ax, s.Tlm, s.Plm, logscale ? log10.(s.lm[:,:,2]') : s.lm[:,:,2]'; colormap=cmap, interpolate=interp); Colorbar(fig[3, 4], hm)
        
        # Big Picture
        if bigpicture[1]
            bigP = vcat(s.Pum, s.Ptz, s.Plm)
            itp = extrapolate(interpolate((s.Pum, s.Tum), s.um[:,:,bigpicture[2]], Gridded(Linear())), 0.0)
            bigH = cat(itp(s.Pum, s.Ttz), s.tz[:,:,bigpicture[2]], s.lm[:,:,bigpicture[2]], dims=1)
            ax = Axis(fig[1:3, 5], ylabel=L"Pressure\;[\mathrm{GPa}]", xlabel=L"Temperature\;[\mathrm{K}]", title=L"Big\;Picture\;(Depleted, wt%)", yreversed=true,
                        xlabelsize=xlabsz, ylabelsize=ylabsz, titlesize=titlesz, xticklabelsize=xticklabsz, yticklabelsize=yticklabsz, xticksize=xticksz, yticksize=yticksz)
            hm = heatmap!(ax, s.Ttz, bigP, logscale ? log10.(bigH') : bigH'; colormap=cmap, interpolate=interp,colorrange= logscale ? (-2., log10(maximum(bigH))) : (minimum(bigH), maximum(bigH))); Colorbar(fig[1:3, 6], hm)
        end
        display(fig)
        savein!="" && CairoMakie.save(savein*".png", fig)
    end

# ======================
# ======= Others =======
# ======================

    function sᴴ²ᴼ_assembler!(map, outHB, n)
        Threads.@threads for i in 1:n
            H₂O = 0.0
            # Depleted (Harzburgite)
            for j in eachindex(outHB[i].SS_vec)
                ("H2O" in outHB[i].SS_vec[j].emNames) && continue
                H₂O += sum(outHB[i].SS_vec[j].Comp[end-1])
            end; map[i, 1] = sum(H₂O)

            # Enriched (Lherzolite)
            H₂O = 0.0
            for j in eachindex(outHB[i+n].SS_vec)
                ("H2O" in outHB[i+n].SS_vec[j].emNames) && continue
                H₂O += sum(outHB[i+n].SS_vec[j].Comp[end-1])
            end; map[i, 2] = sum(H₂O)
        end
    end

    function mesh_vectorization!(P, T, nP, nT, Pv, Tv)
        vP, vT = repeat(P, outer=nT), repeat(T, inner=nP)
        Pv .= vcat(vP, vP); Tv .= vcat(vT, vT)
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

# ======================
# ======= Export =======
# ======================

    function write_output(smap, fmap; fname_s ="StagH2O.dat", fname_fO2="StagfO2.dat", s=true, fO2=true)
        # Array dimensions
        if s
            nP, nT = length(smap.Pum), length(smap.Tum)
            pmap = zeros(Float64, nP, nT)
            open(fname_s, "w") do io
                # Header
                println(io, nP, " ", nT, "\n");
                println(io, smap.Pum[1], " ", smap.Pum[end], " ", smap.Tum[1], " ", smap.Tum[end])
                println(io, smap.Ptz[1], " ", smap.Ptz[end], " ", smap.Ttz[1], " ", smap.Ttz[end])
                println(io, smap.Plm[1], " ", smap.Plm[end], " ", smap.Tlm[1], " ", smap.Tlm[end], "\n")
                # Data
                for (n, slot) in enumerate([2, 1, 2, 1, 2, 1])
                    pmap .= n>4 ? smap.lm[:,:,slot] : n>2 ? smap.tz[:,:,slot] : smap.um[:,:,slot]
                    for i in 1:nP
                        for j in 1:nT
                            print(io, pmap[i,j], " ")
                        end; println(io, "")
                    end; println(io, "")
                end
            end
        end

        if fO2
            nP, nT = size(fmap, 1), size(fmap, 2)
            open(fname_fO2, "w") do io
                # Header
                println(io, nP, " ", nT, "\n");
                println(io, smap.Pum[1], " ", smap.Pum[end], " ", smap.Tum[1], " ", smap.Tum[end], "\n")
                # Data
                for (n, slot) in enumerate([1, 2])
                    for i in 1:nP
                        for j in 1:nT
                            print(io, fmap[i,j], " ")
                        end; println(io, "")
                    end; println(io, "")
                end
            end
        end
    end