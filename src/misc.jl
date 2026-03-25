# Stagno and Frost XCO2 equilibrium (2 Oex per 1 CO₂) (Melt C⁰ -> C⁴⁺ ==> fO₂ to XCO₂ mapping)
function eq_XCO2(logfO2, Pin, Tin)
    P = Pin >= 11. ? 11. : Pin
    T = Tin >= c2k(1600.) ? c2k(1600.) : Tin
    return max(min(10^(logfO2 - 5.44 + 21380/T - 0.078(1e-4P-1)/T), 1.0), 0.0)
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
function partition_Oₑₓ(P::K, T::K; p::K=0.2, ϕ::K=0.01, Rs::K=0.02, Rf::K=0.0, nr=50, niter=100, verbose=false, TOex=nothing) where {K <: Real}

    # Hirschmann 2022 parameters
    a=0.1917; b=-1.961; c=4158.1; ΔCₚ=33.25; T₀=1673.15; y1=-520.46; y2=-185.37; y3=494.39; y4=1838.34; y5=2888.48; y6=3473.68; y7=-4473.6; y8=-1245.09; y9=-1156.86

    # Solid / Molten tracer mass [kg]
    Ms = 1.0e+17
    Mf=ϕ*Ms; Ms-=Mf

    # IDV = ∫ΔVdP for melts
    IDV = solve_∫ΔVdP([P-0.05P, P, P+0.05P],[T-0.05T, T, T+0.05T])[2,2,1]

    # Endmember bulks
    XH  = @SVector [0.4343, 0.4593, 0.0834, 0.0090, 0.0100, 0.0001, 0.0030, 0.0]
    XB  = @SVector [0.5042, 0.0977, 0.0710, 0.1254, 0.1680, 0.0223, 0.0007, 0.0]
    Xox = @SVector ["SiO2", "MgO", "FeO", "CaO", "Al2O3", "Na2O", "Cr2O3", "O"]
    X   = p*XB + (1-p)*XH
    
    # Allocate iterative memory and create bulk structures
    dummyarray = zeros(length(X)) # for Hirschmann calls
    Xdummy = (X=Vector{Float64}(X), Xox=Vector{String}(Xox), mm=get_Xoxmm(Xox)) # For oxidizing calls
    if isnothing(TOex)
        Xs = oxidize_bulk(X, Xox, Rs, Xdummy, FeFormat="FeO_O", wt_out=true, frac=true); 
        Xmo = oxidize_bulk(XB, Xox, Rf, Xdummy, FeFormat="FeO_O", wt_out=true, frac=true);
    end
    Xm = oxidize_bulk(XB, Xox, 0.0, Xdummy, FeFormat="FeO_O", wt_out=true, frac=true);

    # Compute total oxygen budget
    TOₑₓ = isnothing(TOex) ? (Xs.O*Ms + Xmo.O*Mf) : TOex*Ms # kg
    if TOₑₓ==0.0 && verbose
        println("---- Solution (P=$(P)GPa | T=$(T)K | ϕ=$(ϕ) | Rs=$(Rs) | Rf=$(Rf) | Source mix = $p) ----")
        println("Total Oₑₓ budget = 0.0 kg")
        println("TOₑₓ partitioning → [0.0% solid, 0.0% melt]")
    end

    # Generate solid fO₂ space
    Rlist = LinRange(0.00001, 0.5, nr)
    Xlist = Vector{Vector{Float64}}(); 
    for i in 1:nr
        Xl = oxidize_bulk(X, Xox, Rlist[i], Xdummy, wt_out=true, frac=true, FeFormat="FeO_O", onlyvals=true)
        push!(Xlist, [getfield(Xl, Symbol(f)) for f in Xox])
    end
    sOₑₓlist = [Xlist[i][end] for i in 1:nr] # mass fraction Oₑₓ in solid
    # -- Minimizer call
        data    = Initialize_MAGEMin("sb24", verbose=false);
        out = multi_point_minimization(10P*ones(nr), k2c(T)*ones(nr), data, X=Xlist, Xoxides=Vector{String}(Xox), name_solvus=true, sys_in="wt", progressbar=false)
        Finalize_MAGEMin(data);
    # -- Create interpolation object
        # fO2 = [out[i].fO2 for i in eachindex(out)]
        sfO2 = extrapolate(interpolate((sOₑₓlist,), [out[i].fO2 for i in eachindex(out)], Gridded(Linear())), Line())

    # Generate melt fO₂ space
    mfO2_sOₑₓlist = LinRange(0.005TOₑₓ/Ms, 0.995TOₑₓ/Ms, nr)
    # Bisection solver
    minsOₑₓ, maxsOₑₓ, residual, sharedfO2, sOₑₓ, etol, iter = 0.005TOₑₓ/Ms, 0.995TOₑₓ/Ms, 0.0, 0.0, 0.0, 1e-5, 0
    for it = 1:niter
        iter += 1
        sOₑₓ = 0.5(maxsOₑₓ+minsOₑₓ)
        residual = ΔR(sfO2, P, T, Xm, sOₑₓ, TOₑₓ, Ms, Mf, T₀, ΔCₚ, a, b, c, y1, y2, y3, y4, y5, y6, y7, y8, y9, IDV, Xox, dummyarray)
        if (abs(residual)<=etol)
            sharedfO2 = Hirsch(P, T, Xm, sOₑₓ, TOₑₓ, Ms, Mf, T₀, ΔCₚ, a, b, c, y1, y2, y3, y4, y5, y6, y7, y8, y9, IDV, Xox, dummyarray)
            break
        end
        if residual>0.0
            maxsOₑₓ = sOₑₓ
        elseif isnan(residual) || residual<0.0
            minsOₑₓ = sOₑₓ
        end
    end

    if verbose
        println("---- Solution (P=$(P)GPa | T=$(T)K | ϕ=$(ϕ) | Rs=$(Rs) | Rf=$(Rf) | Source mix = $p) ----")
        println("Shared fO₂ = $sharedfO2 |  residual = $(abs(residual))")
        println("Total Oₑₓ budget = $TOₑₓ kg ($(1e2TOₑₓ/Ms)% of solid tracer mass)")
        println("TOₑₓ partitioning → [$(round(1e2(sOₑₓ*Ms/TOₑₓ), digits=4))% solid, $(round(1e2(TOₑₓ-sOₑₓ*Ms)/TOₑₓ, digits=4))% melt]")
        println("Converged in $iter iterations.")
    end

    return 

end

function Hirsch(P, T, Xm, sOₑₓ, TOₑₓ, Ms, Mf, T₀, ΔCₚ, a, b, c, y1, y2, y3, y4, y5, y6, y7, y8, y9, IDV, Xox, dummy)
    # Checks
    @assert ("O"∈Xox && "FeO"∈Xox) "This function requires FeO + O format!"
    # Extract melt fOₑₓ
    mfOₑₓ = (TOₑₓ-sOₑₓ*Ms)/Mf
    # Fill dummy with current composition
    dummy .= [getfield(Xm, Symbol(f)) for f in Xox]
    # Oxidize
    idxO = findfirst(Xox.=="O"); dummy[idxO] = mfOₑₓ
    # mass fraction → molar fraction conversion
    for ox in eachindex(Xox)
        dummy[ox]/=getfield(mm, Symbol(Xox[ox]))
    end
    # Normalize
    dummy./=sum(dummy)
    Xl = Cbulk((; zip(Symbol.(Xox), dummy)...))
    
    (Xl.O>=0.5Xl.FeO) && (return NaN) # Too much oxygen!! Above hard-limit.
    
    # Compute mfO2
    _ln10 = 1/log(10)
    mfO2 = (log10(Xl.O/(Xl.FeO-2Xl.O)) - b - c/T + (ΔCₚ/R*_ln10 * (1 - T₀/T - log(T/T₀))) + IDV/(1e-3R)/T*_ln10 
                        - (1/T)*(y1*Xl.SiO2 + y3*Xl.MgO + y4*Xl.CaO + y5*Xl.Na2O + y8*Xl.SiO2*Xl.Al2O3 + y9*Xl.SiO2*Xl.MgO))/a
    return mfO2
end

function ΔR(sfO2, P, T, Xm, sOₑₓ, TOₑₓ, Ms, Mf, T₀, ΔCₚ, a, b, c, y1, y2, y3, y4, y5, y6, y7, y8, y9, IDV, Xox, dummyarray)
    return sfO2(sOₑₓ) - Hirsch(P, T, Xm, sOₑₓ, TOₑₓ, Ms, Mf, T₀, ΔCₚ, a, b, c, y1, y2, y3, y4, y5, y6, y7, y8, y9, IDV, Xox, dummyarray)

end