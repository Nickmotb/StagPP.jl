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
    function Pitzer_and_Sterner_solver(P_target, Tᵢ, cᵢⱼ; verbose=false, volume=false, ρ₀=1.0, Lbound=1e-3, Ubound=1e3, tol=1e-6, maxiter=300)
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
    function fH₂O(P, Tᵢ, cᵢⱼ)
        fmap = zeros(Float64, length(P), length(T))
        Threads.@threads for t in eachindex(T)
            Tᵢ = T[t]
            cᵢ = get_cᵢ(cᵢⱼ, Tᵢ)
            for p in eachindex(P)
                ρ = Pitzer_and_Sterner_solver(1e9P[p], Tᵢ, cᵢⱼ) # P -> from GPa to bar; T in K
                # Fugacity calculation
                fmap[p,t] = exp(log(ρ) + cᵢ[1]*ρ + (1/(cᵢ[2] + cᵢ[3]*ρ + cᵢ[4]*ρ^2 + cᵢ[5]*ρ^3 + cᵢ[6]*ρ^4) - 1/cᵢ[2]) 
                                - cᵢ[7]/cᵢ[8] * (exp(-cᵢ[8]*ρ) - 1) - cᵢ[9]/cᵢ[10] * (exp(-cᵢ[10]*ρ) - 1) + P[p]/(ρ*1e6R*Tᵢ) + log(1e6R*Tᵢ) - 1)/1e5
            end
        end
        return fmap./1e-4 # Convert back to GPa
    end

# ======================
# ==== Functions =======
# ======================

    function Xᵢ_ig2stx(Xᵢ, Clist_ig, Clist_stx)
        # Map indices
        ig2stx_map = Dict{String, Int64}(); [ig2stx_map[oxide] = i for (i, oxide) in enumerate(Clist_ig)]
        return [Xᵢ[ig2stx_map["SiO2"]], Xᵢ[ig2stx_map["Al2O3"]], Xᵢ[ig2stx_map["CaO"]], Xᵢ[ig2stx_map["MgO"]], Xᵢ[ig2stx_map["FeO"]], Xᵢ[ig2stx_map["Na2O"]]]
    end