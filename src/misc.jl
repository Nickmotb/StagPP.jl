# Stagno and Frost XCO2 equilibrium (2 Oex per 1 CO₂) | cOₑₓ in mass fraction of TOₑₓ
function cOₑₓ_to_fO2(cOₑₓ, Pin, Tin, s, Φ; giveXCO₂=false)
    # Checks
    @assert cOₑₓ>=0.0 "Mass of carbon being oxidized for XCO₂ call must be positive! (Don't call when reducing)"
    # Mass of Oₑₓ used for oxidation as mass fraction of molten tracer
    XOₑₓ = Φ*cOₑₓ
    # Convert Oₑₓ to molar fraction, and assess amount of XCO₂ in mols [0.5XCO₂ for each XOₑₓ]
    XCO₂ = 0.5XOₑₓ/(s+XOₑₓ)
    # Limit to rexplored ranges
    P = Pin >= 11. ? 11. : Pin
    T = Tin >= c2k(1600.) ? c2k(1600.) : Tin
    # Compute logfO₂
    if giveXCO₂
        return XCO₂
    else
        return 5.44 - 21380/T + 0.078(1e5P-1)/T + log10(XCO₂)
    end
end

# Model:
# Solution space is [P-T-TOₑₓ] -> looking to find ∂Oₑₓ/∂Mmelt

# Inputs:
# - Parent source bulk composition [Xs] -> defined as Xs(p) = p*XB + (1-p)*XH
# - Local parcel [P] and [T]
# - Total oxygen budget [TOₑₓ] collected from solid (and molten if present) tracer

# Tools:
# - Hirschmann 2022 melt mapping from XFe₂O₃ (Oₑₓ) ↔ fO₂
# - Stixrude and Bertelloni 2024 (MAGEMin) solid mapping from XFe₂O₃ (Oₑₓ) ↔ fO₂
# - Stagno and Frost 2010 parameterization of melt EDDOG2 buffer fO₂ ↔ XCO₂
function partition_Oₑₓ(P::K, T::K; p::K=0.2, ϕ::K=0.01, Rs::K=0.02, Rf::K=0.0, Ctot::K=0.1, nr=50, niter=100, verbose=false, TOex=nothing, data=nothing, respace=(false, 20)) where {K <: Real}

    # Endmember bulks
    XH      = @SVector [0.4343, 0.4593, 0.0834, 0.0090, 0.0100, 0.0001, 0.0030, 0.0]
    XB      = @SVector [0.5042, 0.0977, 0.0710, 0.1254, 0.1680, 0.0223, 0.0007, 0.0]
    Xox     = @SVector ["SiO2", "MgO", "FeO", "CaO", "Al2O3", "Na2O", "Cr2O3", "O"]
    SymXox  = Tuple(Symbol.(Xox))
    mmXox   = @SVector [mm.SiO2, mm.MgO, mm.FeO, mm.CaO, mm.Al2O3, mm.Na2O, mm.Cr2O3, mm.O]
    X       = p*XB + (1-p)*XH
    molXB   = XB./mmXox

    # Hirschmann 2022 parameters
    a=0.1917; b=-1.961; c=4158.1; ΔCₚ=33.25; T₀=1673.15; y1=-520.46; y2=-185.37; y3=494.39; y4=1838.34; y5=2888.48; y6=3473.68; y7=-4473.6; y8=-1245.09; y9=-1156.86

    # Solid / Molten tracer mass [kg]
    Ms = 1.0e+17
    Mf=ϕ*Ms; Ms-=Mf; Mc=Ctot*Ms
    # molMf = sum(1e3Mf.*XB./mmXox) # Mf in mols

    # Mols of available carbon
    # molCav = 1e3Mc/mm.C

    # IDV = ∫ΔVdP for melts
    IDV = solve_∫ΔVdP([P-0.05P, P, P+0.05P],[T-0.05T, T, T+0.05T])[2,2,1]

    # Flag whether to manage MAGEMin initialisation and finalization
    flag = isnothing(data)

    # Index of oxygen component
    idxO = findfirst(Xox.=="O")
    
    # Allocate iterative memory and create bulk structures
    J           = @SMatrix zeros(2,2)   # Jacobian
    sol         = @SVector zeros(2)     # Solution Vector
    dummy       = zeros(length(XB))   # for Hirschmann calls
    Xdummy      = (X=Vector{Float64}(X), Xox=Vector{String}(Xox), mm=get_Xoxmm(Xox)) # For oxidizing calls
    if isnothing(TOex)
        Xs = oxidize_bulk(X, Rs, Xdummy, FeFormat="FeO_O", wt_out=true, frac=true, SymXox=SymXox); 
        Xmo = oxidize_bulk(XB, Rf, Xdummy, FeFormat="FeO_O", wt_out=true, frac=true, SymXox=SymXox);
    end
    Xm = oxidize_bulk(XB, 0.0, Xdummy, FeFormat="FeO_O", wt_out=true, frac=true);

    # Compute total oxygen budget
    TOₑₓ = isnothing(TOex) ? (Xs.O*Ms + Xmo.O*Mf) : TOex*Ms # kg
    if TOₑₓ==0.0 && verbose
        println("---- Solution (P=$(P)GPa | T=$(T)K | ϕ=$(ϕ) | Rs=$(Rs) | Rf=$(Rf) | Source mix = $p) ----")
        println("Total Oₑₓ budget = 0.0 kg")
        println("TOₑₓ partitioning → [0.0% solid, 0.0% melt]")
        return
    end

    # Generate solid fO₂ space
    Rlist = LinRange(0.00001, 0.05, nr)
    Xlist = Vector{Vector{Float64}}(undef, nr); 
    sOₑₓlist = zeros(nr)
    for i in 1:nr
        Xl = oxidize_bulk(X, Rlist[i], Xdummy, wt_out=true, frac=true, FeFormat="FeO_O", SymXox=SymXox)
        Xlist[i] = [getfield(Xl, f) for f in SymXox]
        sOₑₓlist[i] = Xl.O
    end
    # -- Minimizer call
        flag && (data = Initialize_MAGEMin("sb24", verbose=false))
        out = multi_point_minimization(10P*ones(nr), k2c(T)*ones(nr), data, X=Xlist, Xoxides=Vector{String}(Xox), name_solvus=true, sys_in="wt", progressbar=false)
        flag && Finalize_MAGEMin(data);
    # -- Create interpolation object
        sfO2 = extrapolate(interpolate((sOₑₓlist.*Ms./TOₑₓ,), [out[i].fO2 for i in eachindex(out)], Gridded(Linear())), Line())
        sample_sOlist = LinRange(1e-7, 0.995, 250)
        sample_sOlist05 = 0.5(sample_sOlist[1:end-1] + sample_sOlist[2:end])
        sampled_sfO2 = sfO2(sample_sOlist)
    # -- Solid partial derivative
        ∂Sᵢ = extrapolate(interpolate((sample_sOlist05,), ∂S∂sOₑₓ(sampled_sfO2, sample_sOlist), Gridded(Linear())), Line())

    # Compute Jacobian variables
    Φ = evΦ(TOₑₓ, Mf, mm.O)                                          # Proportion of TOₑₓ → normalized XOₑₓ component
    s = sum(molXB)                                                   # Sum of non-normalized molar components
    _ln10, _T = 1/log(10), 1/T                                       
    Ys1 = (y1*molXB[1] + y3*molXB[2] + y4*molXB[4] + y5*molXB[6])*_T  # Sum of linear parameterized molar components
    Ys2 = molXB[1]*(y8*molXB[5] + y9*molXB[2])*_T                     # Sum of non-linear parameterized molar components

    if !respace[1]
        # Newton Solver
        minsOₑₓ, maxsOₑₓ =  1e-7, 0.99
        mincOₑₓ, maxcOₑₓ = 1e-7, 0.99 #(2molCav*mm.O)/TOₑₓ
        x = sol + [0.7(minsOₑₓ+maxsOₑₓ), 0.1(mincOₑₓ+maxcOₑₓ)]
        etol, damp = 1e-1, 1.0
        XCO₂ = 0.0
        mat = zeros(niter, 2)
            for it in 1:niter
                # Evaluate current stage
                sOₑₓ, cOₑₓ = x
                # Reset dummy
                dummy .= molXB
                # Compute residual
                Fx = Rx(sfO2, P, T, sOₑₓ, cOₑₓ, TOₑₓ, Mf, T₀, ΔCₚ, a, b, c, y1, y3, y4, y5, y8, y9, IDV, SymXox, dummy, idxO, _ln10, _T, s, Φ, XCO₂)
                # Exit if below tolerance
                if sum(abs.(Fx))<=etol
                    if verbose
                        println("---- Solution (P=$(P)GPa | T=$(T)K | ϕ=$(ϕ) | Rs=$(Rs) | Rf=$(Rf) | Source mix = $p) ----")
                        println("Shared fO₂ = $(sfO2(sOₑₓ)) |  residual = $(sum(abs.(Fx)))")
                        println("Total Oₑₓ budget = $TOₑₓ kg ($(1e2TOₑₓ/Ms)% of solid tracer mass)")
                        println("TOₑₓ partitioning → [$(round(1e2x[1], digits=4))% solid, $(round(1e2(1 - x[1] - x[2]), digits=4))% melt], $(round(1e2x[2], digits=4))% used to oxidize C → CO₂")
                        println("Melt XCO₂ = $(cOₑₓ_to_fO2(cOₑₓ, P, T, s, Φ, giveXCO₂=true))")
                        println("Converged in $it iterations.")
                    end
                    break
                end
                mat[it,:] .= Fx;
                # Partial derivatives
                α = evα(sOₑₓ, cOₑₓ, Φ); θ = evθ(α, s); θ² = θ^2
                ∂S = ∂Sᵢ(sOₑₓ)
                ∂C = ∂C∂cOₑₓ(cOₑₓ, Φ, s, _ln10)
                ∂M∂iOₑₓ = ∂M∂xOₑₓ(Φ, Ys1, Ys2, α, θ, θ², _ln10, a, molXB[3])
                # Jacobian inverse
                J⁻¹ = J + inv([(∂S-∂M∂iOₑₓ) -∂M∂iOₑₓ; ∂S         -∂C])
                x = x - J⁻¹*Fx*damp
                # Clamp
                x = SA[ min(max(x[1], minsOₑₓ), maxsOₑₓ), min(max(x[2], mincOₑₓ), maxcOₑₓ)]
                # Output if not converged
                if verbose && it==niter
                    println("---- Solution (P=$(P)GPa | T=$(T)K | ϕ=$(ϕ) | Rs=$(Rs) | Rf=$(Rf) | Source mix = $p) ----")
                    println("Shared fO₂ = $(sfO2(sOₑₓ)) |  residual = $(sum(abs.(Fx)))")
                    println("Total Oₑₓ budget = $TOₑₓ kg ($(1e2TOₑₓ/Ms)% of solid tracer mass)")
                    println("TOₑₓ partitioning → [$(round(1e2x[1], digits=4))% solid, $(round(1e2(1 - x[1] - x[2]), digits=4))% melt], $(round(1e2x[2], digits=4))% used to oxidize C → CO₂")
                    println("Melt XCO₂ = $(cOₑₓ_to_fO2(cOₑₓ, P, T, s, Φ, giveXCO₂=true)))")
                    println("Did not converge in $it iterations.")
                end
            end
        mOₑₓ = 1 - x[1] - x[2]
        fig = Figure(size=(1000, 700))
        ax = Axis(fig[1,1]); plot!(ax, 1:niter, mat[:,1]); plot!(ax, 1:niter, mat[:,2], color=:red)
        display(fig)

        return x[1], mOₑₓ, x[2] 

    else

        # Export 1D sensitivity to P | T | ϕ | TOₑₓ
        v = zeros(respace[2])
        sOr = LinRange(0.005, 0.995, respace[2])
        for i in 1:respace[2]
            v[i] = ΔR(sfO2, T, Xm, TOₑₓ*sOr[i]/Ms, TOₑₓ, Ms, Mf, T₀, ΔCₚ, a, b, c, y1, y3, y4, y5, y8, y9, IDV, SymXox, dummyarray, mmXox)
        end     
        return v

    end

end

function Hirsch(T, sOₑₓ, cOₑₓ, TOₑₓ, Mf, T₀, ΔCₚ, a, b, c, y1, y3, y4, y5, y8, y9, IDV, SymXox, dummy, idxO, _ln10, _T, s, Φ)
    # Checks
    @assert (:O∈SymXox && :FeO∈SymXox) "This function requires FeO + O format!"
    # Extract melt fOₑₓ
    mfOₑₓ = Φ*(1-sOₑₓ-cOₑₓ)
    # Oxidize
    dummy[idxO] = mfOₑₓ
    # Normalize
    dummy./=(s+mfOₑₓ)
    # Construct bulk and assess hardlimit
    Xl = Cbulk((; zip(SymXox, dummy)...))
    (Xl.O>=0.5Xl.FeO || Xl.O == 0.0) && (return -Inf) # Too much oxygen!! Above hard-limit.
    # Compute mfO2
    mfO2 = (log10(Xl.O/(Xl.FeO-2Xl.O)) - b - c*_T + (ΔCₚ/R*_ln10 * (1 - T₀*_T - log(T/T₀))) + IDV/(1e-3R)*_T*_ln10 
                        - _T*(y1*Xl.SiO2 + y3*Xl.MgO + y4*Xl.CaO + y5*Xl.Na2O + y8*Xl.SiO2*Xl.Al2O3 + y9*Xl.SiO2*Xl.MgO))/a
    return mfO2
end

function Rx(sfO2, P, T, sOₑₓ, cOₑₓ, TOₑₓ, Mf, T₀, ΔCₚ, a, b, c, y1, y3, y4, y5, y8, y9, IDV, SymXox, dummy, idxO, _ln10, _T, s, Φ, XCO₂)
    Rₛ = sfO2(sOₑₓ)
    return SA[Rₛ - Hirsch(T, sOₑₓ, cOₑₓ, TOₑₓ, Mf, T₀, ΔCₚ, a, b, c, y1, y3, y4, y5, y8, y9, IDV, SymXox, dummy, idxO, _ln10, _T, s, Φ)
              Rₛ - cOₑₓ_to_fO2(cOₑₓ, P, T, s, Φ)]
end

function Rx1D(sfO2, T, sOₑₓ, cOₑₓ, TOₑₓ, Mf, T₀, ΔCₚ, a, b, c, y1, y3, y4, y5, y8, y9, IDV, SymXox, dummy, idxO, _ln10, _T, s)
    Rₛ = sfO2(sOₑₓ)
    return Rₛ - Hirsch(T, sOₑₓ, cOₑₓ, TOₑₓ, Mf, T₀, ΔCₚ, a, b, c, y1, y3, y4, y5, y8, y9, IDV, SymXox, dummy, idxO, _ln10, _T, s)
end

function P_T_ϕ_TOₑₓ_Sspace(Pr::K, Tr::K, ϕr::K, TOr::K, data; cutϕ::Float64=-1.0, cutTO::Float64=-1.0, p::Float64=0.2, niter::Int64=100, nr::Int64=50) where {K}
    # Resolutions
    nP, nT, nϕ, nTO = length(Pr), length(Tr), length(ϕr), length(TOr)
    iϕ05, iTO05 = Int(floor(0.5nϕ)), Int(floor(0.5nTO))
    # Correct cut values
    cutϕ = (cutϕ==-1) ? ϕr[iϕ05] : min(max(cutϕ, first(ϕr)), last(ϕr))
    cutTO = (cutTO==-1) ? TOr[iTO05] : min(max(cutTO, first(TOr)), last(TOr))
    # Construct map
    mapϕ, mapTO, nmax = zeros(nP, nT, nϕ), zeros(nP, nT, nTO), max(nϕ, nTO)
    for n in 1:nmax
        println("n = $n / $nmax")
        for ip in 1:nP
            for it in 1:nT
                (n<=nϕ)  && (mapϕ[ip, it, n]  = partition_Oₑₓ(Pr[ip], Tr[it]; p=p, ϕ=ϕr[n], Rs=0.0, Rf=0.0, nr=nr, niter=niter, verbose=false, TOex=cutTO, data=data))
                (n<=nTO) && (mapTO[ip, it, n] = partition_Oₑₓ(Pr[ip], Tr[it]; p=p, ϕ=cutϕ, Rs=0.0, Rf=0.0, nr=nr, niter=niter, verbose=false, TOex=TOr[n], data=data))
            end
        end
    end
    fig = Figure(size=(1200, 700))
    ax = Axis3(fig[1,1], xlabel="Pressure [GPa]", ylabel="Temperature [K]", zlabel="ϕ [%]")
    ax2 = Axis3(fig[1,2], xlabel="Pressure [GPa]", ylabel="Temperature [K]", zlabel="TOₑₓ [wt% Mₜ]")
    for i in 1:max(nϕ, nTO)
        (i<=nϕ) && surface!(ax, Pr, Tr, ϕr[i]*ones(nP,nT), color=1.0.-mapϕ[:,:,i], colormap=:vik100)
        (i<=nTO) && surface!(ax2, Pr, Tr, TOr[i]*ones(nP,nT), color=1.0.-mapTO[:,:,i], colormap=:vik100)
        (i==1) && Colorbar(fig[2,1:2], colorrange=(0.0, 100.0), colormap=:vik100, label=L"TO_{ex}\;in\;melt\;\mathrm{[\%]}", labelsize=20, vertical=false)
    end
    display(fig)
    return mapϕ, mapTO
end

function P_T_ϕ_TOₑₓ_Rspace(Pr::K, Tr::K, ϕr::K, TOr::K, data; p::Float64=0.2, nres=20) where {K}
    # Resolutions
    nP, nT, nϕ, nTO = length(Pr), length(Tr), length(ϕr), length(TOr)
    defP, defT, defϕ, defTO = 3.0, 1600., 0.02, 3e-4
    # Construct map
    mapP = zeros(nP, nres); mapT = zeros(nP, nres)
    mapϕ = zeros(nP, nres); mapTO = zeros(nP, nres)
    nmax = max(nP, nT, nϕ, nTO)
    sOr = LinRange(0.005, 0.995, nres)
    for n in 1:nmax
        println("n = $n / $nmax")
        (n<=nϕ)  && (mapϕ[n, :]  .= partition_Oₑₓ(defP, defT; p=p, ϕ=ϕr[n], verbose=false, TOex=defTO, data=data, respace=(true, nres)))
        (n<=nTO) && (mapTO[n, :] .= partition_Oₑₓ(defP, defT; p=p, ϕ=defϕ, verbose=false, TOex=TOr[n], data=data, respace=(true, nres)))
        (n<=nP) && (mapP[n, :] .= partition_Oₑₓ(Pr[n], defT; p=p, ϕ=defϕ, verbose=false, TOex=defTO, data=data, respace=(true, nres)))
        (n<=nT) && (mapT[n, :] .= partition_Oₑₓ(defP, Tr[n]; p=p, ϕ=defϕ, verbose=false, TOex=defTO, data=data, respace=(true, nres)))
    end
    # Set asbolutes
    mapP.=abs.(mapP)
    mapT.=abs.(mapT)
    mapϕ.=abs.(mapϕ)
    mapTO.=abs.(mapTO)

    fig = Figure(size=(1200, 700))
    ax = Axis3(fig[1,1], xlabel="Pressure [GPa]", ylabel="sOₑₓ [%TOₑₓ]", zlabel="ΔR")
        surface!(ax, Pr, sOr, mapP, colormap=:Purples, alpha=1.0)
        wireframe!(ax, Pr, sOr, mapP, color=:black)
    ax = Axis3(fig[1,2], xlabel="Temperature [GPa]", ylabel="sOₑₓ [%TOₑₓ]", zlabel="ΔR")
        surface!(ax, Tr, sOr, mapT, colormap=:Purples, alpha=1.0)
        wireframe!(ax, Tr, sOr, mapT, color=:black)
    ax = Axis3(fig[2,1], xlabel="ϕ [%]", ylabel="sOₑₓ [%TOₑₓ]", zlabel="ΔR")
        surface!(ax, ϕr, sOr, mapϕ, colormap=:Purples, alpha=1.0)
        wireframe!(ax, ϕr, sOr, mapϕ, color=:black)
    ax = Axis3(fig[2,2], xlabel="TOₑₓ [% of Mₜ]", ylabel="sOₑₓ [%TOₑₓ]", zlabel="ΔR")
        surface!(ax, 1e2TOr, sOr, mapTO, colormap=:Purples, alpha=1.0)
        wireframe!(ax, 1e2TOr, sOr, mapTO, color=:black)
    display(fig)
end

# Variables
@inline evΦ(TOₑₓ, Mf, mmO) = TOₑₓ/Mf/mmO            # Conversion factor
@inline evα(sOₑₓ, cOₑₓ, Φ) = Φ*(1 - sOₑₓ - cOₑₓ)    # mass fraction of TOₑₓ → non-normalized XFe₂O₃
@inline evθ(α, s) = s - α                           # 
@inline evαₙ(α, θ) = α/θ

# Partial derivatives
∂S∂sOₑₓ(sfO2, sOlist) = @views diff(sfO2)./diff(sOlist)
∂C∂cOₑₓ(cOₑₓ, Φ, s, _ln10) = s*_ln10/(cOₑₓ*(s + Φ*cOₑₓ))
∂M∂xOₑₓ(Φ, Ys1, Ys2, α, θ, θ², _ln10, a, XFeO) = Φ/θ²/a*(Ys1 + 2Ys2/θ - _ln10*(θ²/α + 2(θ+α)/(XFeO - 2α)))
