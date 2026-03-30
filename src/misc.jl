# Model:
# Assess partitioning of TOₑₓ into solid and molten tracer + oxidation reduction of parcel carbon through equilibrium XCO₂
#
# Inputs:
# - Pressure (P)
# - Temperature (T)
# - Total oxygen budget either as wt% of solid (TOₑₓ) or solid and melt Fe³⁺/Feᵀ ratios (Rs, Rf)
# - Melt tracer initial XCO₂ (iXCO₂)
# - Available reduced carbon in the parcel as wt% of solid tracer (avRC)
#
# Independent variables:
# - Equilibrium mass fraction of TOₑₓ in the solid (sOₑₓ)
# - Equilibrium mass fraction of TOₑₓ in the melt  (mOₑₓ)
# - Equilibrium molar XCO₂ in the melt (XCO₂)
#
# Equations:
# (1) Equilibrium constraint    :  solid fO₂(P,T,sOₑₓ) = melt fO₂(P,T,mOₑₓ)
# (2) Equilibrium constraint    :  solid fO₂(P,T,sOₑₓ) = EDDOG fO₂(P,T,XCO₂) <-- Maybe exclude residual if XCO₂<0.0 or XCO₂>1.0
# (3) Mass conservation         :  1 = sOₑₓ + mfOₑₓ + Oₑₓ_in_CO₂(XCO₂)  →  1 = sOₑₓ + mfOₑₓ + Φ*cα/(sw+cα)
#
# Jacobian : [∂(1)∂sOₑₓ ∂(1)∂mOₑₓ ∂(1)∂XCO₂         [∂S∂sOₑₓ     -∂M∂mOₑₓ        0
#             ∂(2)∂sOₑₓ ∂(2)∂mOₑₓ ∂(2)∂XCO₂    =     ∂S∂sOₑₓ        0       -∂C∂XCO₂
#             ∂(3)∂sOₑₓ ∂(3)∂mOₑₓ ∂(3)∂XCO₂]           -1           -1      -∂[Φ*cα/(sw+cα)]∂XCO₂]
#
# Variable extentions:
# cα = 2XCO₂*mm.O  # non-normalized mass fraction Oₑₓ in melt
# Φ  = Mf/TOₑₓ
# sw, s = ∑(oxᵢ*mmᵢ), ∑(oxᵢ/mmᵢ)
# Oₑₓ_in_CO₂ = Φ*cα/(sw+cα)
#
# Tools:
# - Hirschmann 2022 melt mapping from XFe₂O₃ (Oₑₓ) ↔ fO₂
# - Stixrude and Bertelloni 2024 (MAGEMin) solid mapping from XFe₂O₃ (Oₑₓ) ↔ fO₂
# - Stagno and Frost 2010 parameterization of melt EDDOG2 buffer fO₂ ↔ XCO₂
function partition_Oₑₓ(P::K, T::K; p::K=0.2, ϕ::K=0.01, Rs::K=0.02, Rf::K=0.0, Ctot::K=0.1, iXCO₂::K=0.01, nr=50, niter=100, verbose=false, TOex=nothing, data=nothing, respace=(false, 20), plotevo=false) where {K <: Real}

    # Endmember bulks
    XH      = @SVector [0.4343, 0.4593, 0.0834, 0.0090, 0.0100, 0.0001, 0.0030, 0.0]
    XB      = @SVector [0.5042, 0.0977, 0.0710, 0.1254, 0.1680, 0.0223, 0.0007, 0.0]
    Xox     = @SVector ["SiO2", "MgO", "FeO", "CaO", "Al2O3", "Na2O", "Cr2O3", "O"]
    SymXox  = Tuple(Symbol.(Xox))
    mmXox   = @SVector [mm.SiO2, mm.MgO, mm.FeO, mm.CaO, mm.Al2O3, mm.Na2O, mm.Cr2O3, mm.O]
    X       = p*XB + (1-p)*XH
    molXB   = XB./mmXox; molXB = molXB ./ sum(molXB)
    molX    = X./mmXox;  molX  = molX ./sum(molX)

    # Hirschmann 2022 parameters
    a=0.1917; b=-1.961; c=4158.1; ΔCₚ=33.25; T₀=1673.15; y1=-520.46; y2=-185.37; y3=494.39; y4=1838.34; y5=2888.48; y8=-1245.09; y9=-1156.86

    # Solid / Molten tracer mass [kg]
    Ms = 1.0e+17
    Mf=ϕ*Ms; Ms-=Mf; Mc=Ctot*Ms
    molMf = sum(1e3Mf.*XB./mmXox) # Mf in mols

    # Mols of available carbon
    molCav = 1e3Mc/mm.C

    # IDV = ∫ΔVdP for melts
    IDV = solve_∫ΔVdP([P-0.05P, P, P+0.05P],[T-0.05T, T, T+0.05T])[2,2,1]

    # Flag whether to manage MAGEMin initialisation and finalization
    flag = isnothing(data)
    converged = false

    # Index of oxygen component
    idxO = findfirst(Xox.=="O")
    
    # === Allocate iterative memory and create bulk structures
        J           = @SMatrix zeros(3,3)   # Jacobian
        sol         = @SVector zeros(3)     # Solution Vector
        dummy       = zeros(length(XB))   # for Hirschmann calls
        Xdummy      = (X=Vector{Float64}(X), Xox=Vector{String}(Xox), mm=get_Xoxmm(Xox)) # For oxidizing calls
        if isnothing(TOex)
            Xs = oxidize_bulk(X, Rs, Xdummy, FeFormat="FeO_O", wt_out=true, frac=true, SymXox=SymXox); 
            Xmo = oxidize_bulk(XB, Rf, Xdummy, FeFormat="FeO_O", wt_out=true, frac=true, SymXox=SymXox);
        end
        Xm = oxidize_bulk(XB, 0.0, Xdummy, FeFormat="FeO_O", wt_out=true, frac=true);

    # Sum TOₑₓ contributions from solid and molten tracer Oₑₓ
    verb_flag = isnothing(TOex) ? -1 : TOex
    TOₑₓ = isnothing(TOex) ? (Xs.O*Ms + Xmo.O*Mf) : TOex*Ms # kg
    if TOₑₓ==0.0 && verbose
        println("---- Solution (P=$(P)GPa | T=$(T)K | ϕ=$(ϕ) | Rs=$(Rs) | Rf=$(Rf) | Source mix = $p) ----")
        println("Total Oₑₓ budget = 0.0 kg")
        println("TOₑₓ partitioning → [0.0% solid, 0.0% melt]")
        return
    end

    # Compute Jacobian variables
    s, sw = sum(molXB), sum(XB)           
    αᵢ = evcα(iXCO₂); TOₑₓ += αᵢ/(sw+αᵢ)*Mf # Add XCO₂ contribution to TOₑₓ
    Φ, Φₘ = evΦ(TOₑₓ, Mf), evΦₘ(TOₑₓ, Mf)                                     
    _ln10, _T = 1/log(10), 1/T                                       
    Ys1 = (y1*molXB[1] + y3*molXB[2] + y4*molXB[4] + y5*molXB[6])*_T  # Sum of linear parameterized molar components
    Ys2 = molXB[1]*(y8*molXB[5] + y9*molXB[2])*_T                     # Sum of non-linear parameterized molar components

    # sOₑₓ and mOₑₓ caps
    maxsOₑₓ = min((0.5molX[3]*mm.O)/(sum(X) + 0.5molX[3]*mm.O)*(Mf/TOₑₓ), 1.0)
    maxmOₑₓ = min((0.5molXB[3]*mm.O)/(sum(XB) + 0.5molXB[3]*mm.O)*(Mf/TOₑₓ), 1.0)

    # === Generate solid fO₂ space
        Rlist = LinRange(0.00001, 0.15, nr)
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
            sample_sOlist = LinRange(1e-9, maxsOₑₓ, 100)
            sample_sOlist05 = 0.5(sample_sOlist[1:end-1] + sample_sOlist[2:end])
            sampled_sfO2 = sfO2(sample_sOlist)
        # -- Solid partial derivative
            ∂Sᵢ = extrapolate(interpolate((sample_sOlist05,), ∂S∂sOₑₓ(sampled_sfO2, sample_sOlist), Gridded(Linear())), Line())

    if !respace[1]
        # Newton Solver
        minsOₑₓ =  1e-7
        minmOₑₓ = 1e-7
        minXCO₂ = 1e-7; maxXCO₂ = max(min(1.0, iXCO₂+molCav/molMf), minXCO₂)
        x = sol + [0.7(minsOₑₓ+maxsOₑₓ), 0.1(minmOₑₓ+maxmOₑₓ), 1e-4]
        etol, damp, aR = 1e-2, 0.35, Inf
        plotevo && (mat = zeros(niter, 3, 3); Rs = zeros(3))
        for it in 1:niter
            # Evaluate current stage
            sOₑₓ, mOₑₓ, XCO₂ = x
            # Iteration variables
            α   = evα(mOₑₓ, Φₘ)
            θₘ  = evθₘ(α, s)
            cα  = evcα(XCO₂)
            θ   = evθ(cα, sw)
            # Compute residual
            Fx = Rx(sfO2, P, T, sOₑₓ, mOₑₓ, XCO₂, T₀, ΔCₚ, a, b, c, y1, y3, y4, y5, y8, y9, IDV, SymXox, dummy, idxO, _ln10, _T, s, Φ, sw, Φₘ, cα, molXB)
            if plotevo
                mat[it,:,1] .= Fx;
                mat[it,:,3] .= x
                mat[it,1,2] = sfO2(sOₑₓ)
                mat[it,2,2] = Hirsch(T, mOₑₓ, T₀, ΔCₚ, a, b, c, y1, y3, y4, y5, y8, y9, IDV, SymXox, dummy, idxO, _ln10, _T, s, Φₘ, molXB)
                mat[it,3,2] = XCO₂_to_fO2(XCO₂, P, T)
            end
            # Exit if below tolerance
            aR = (x[3]==maxXCO₂ || x[3]==minXCO₂) ? max(abs(Fx[1]), abs(Fx[3])) : maximum(abs.(Fx))
            if aR<=etol
                converged = true
                if verbose
                    if verb_flag==-1
                        println("---- Solution (P=$(P)GPa | T=$(T)K | ϕ=$(ϕ) | Rs=$(Rs) | Rf=$(Rf) | Source mix = $p | Total carbon = $Ctot) ----")
                    else
                        println("---- Solution (P=$(P)GPa | T=$(T)K | ϕ=$(ϕ) | TOₑₓ=$TOₑₓ | Source mix = $p) ----")
                    end
                    println("Shared fO₂ = $(sfO2(sOₑₓ)) |  residual = $aR")
                    println("Total Oₑₓ budget = $(round(TOₑₓ, digits=4)) kg ($(round((1e2TOₑₓ/Ms), digits=4))% of solid tracer mass)")
                    println("TOₑₓ partitioning → [$(round(1e2x[1], digits=4))% solid, $(round(1e2x[2], digits=4))% melt] stored as Fe₂O₃, $(round(1e2(1 - x[1] - x[2]), digits=4))% stored as melt CO₂")
                    println("Melt XCO₂ = $(x[3])")
                    println("Converged in $it iterations.")
                end
                break
            end
            # Partial derivatives
            ∂S = ∂Sᵢ(sOₑₓ)
            ∂C = ∂C∂XCO₂(XCO₂, _ln10)
            ∂M = ∂M∂mOₑₓ(Φₘ, Ys1, Ys2, α, θₘ, _ln10, a, molXB[3])
            ∂3 = ∂3∂XCO₂(Φ, θ, cα)
            # Jacobian inverse
            J⁻¹ = J + inv([ ∂S -∂M 0.0
                            ∂S 0.0 -∂C
                            -1 -1 -∂3])
            x = x - J⁻¹*Fx*damp

            # Clamp
            x = SA[ min(max(x[1], minsOₑₓ), maxsOₑₓ), min(max(x[2], minmOₑₓ), maxmOₑₓ), min(max(x[3], minXCO₂), maxXCO₂)]
            # Store residuals
            plotevo && (Rs .= Fx)
            # Output if not converged
            if verbose && it==niter
                if verb_flag==-1
                    println("---- Solution (P=$(P)GPa | T=$(T)K | ϕ=$(ϕ) | Rs=$(Rs) | Rf=$(Rf) | Source mix = $p | Total carbon = $Ctot) ----")
                else
                    println("---- Solution (P=$(P)GPa | T=$(T)K | ϕ=$(ϕ) | TOₑₓ=$TOₑₓ | Source mix = $p) ----")
                end
                println("Shared fO₂ = $(sfO2(sOₑₓ)) |  residual = $aR")
                println("Total Oₑₓ budget = $(round(TOₑₓ, digits=4)) kg ($(round((1e2TOₑₓ/Ms), digits=4))% of solid tracer mass)")
                println("TOₑₓ partitioning → [$(round(1e2x[1], digits=4))% solid, $(round(1e2x[2], digits=4))% melt] stored as Fe₂O₃, $(round(1e2(1 - x[1] - x[2]), digits=4))% stored as melt CO₂")
                println("Melt XCO₂ = $(x[3])")
                println("Did not converge in $it iterations.")
            end
        end
        if plotevo
            solidclr = :black
            meltclr  = :red
            co2clr   = :green
            # Plot evolution
            iend = findfirst(mat[:,1,1].==0.0); isnothing(iend) ? (iend=niter) : (iend-=1)
            mOₑₓ = 1 - x[1] - x[2]
            fig = Figure(size=(1800, 800))
            ax = Axis(fig[1,1], ylabel=L"Solution\;residual", xlabel=L"Iterations", xgridvisible=false, ygridvisible=false, ylabelsize=25, xlabelsize=25); 
            scatterlines!(ax, 1:iend, mat[1:iend,1,1], color=meltclr, label="(eq. 1) Solid ↔ melt fO₂ equilibrium (R₁ = $(round(Rs[1], digits=4)))",marker=:rect,strokewidth=1.1)
            scatterlines!(ax, 1:iend, mat[1:iend,2,1], color=co2clr,label="(eq. 2) Solid ↔ EDDOG fO₂ equilibrium (R₂ = $(round(Rs[2], digits=4)))",marker=:rect,strokewidth=1.1)
            scatterlines!(ax, 1:iend, mat[1:iend,3,1], color=:orange,label="(eq. 3) Mass conservation (R₃ = $(round(Rs[3], digits=4)))",marker=:rect,strokewidth=1.1)
            axislegend(ax, position=:rt, "Convergence ϵ = $etol")
            converged && ylims!(ax, -0.5, 1.0)

            ax = Axis(fig[2,1], ylabel=L"log\;fO_2", xlabel=L"Iterations", xgridvisible=false, ygridvisible=false, ylabelsize=25, xlabelsize=25); 
            scatterlines!(ax, 1:iend, mat[1:iend,1,2], color=solidclr, label="Solid (fO₂ = $(round(mat[iend,1,2], digits=4)))",marker=:rect,strokewidth=1.1)
            scatterlines!(ax, 1:iend, mat[1:iend,2,2], color=meltclr,label="Melt (fO₂ = $(round(mat[iend,2,2], digits=4)))",marker=:rect,strokewidth=1.1)
            scatterlines!(ax, 1:iend, mat[1:iend,3,2], color=co2clr,label="EDDOG (fO₂ = $(round(mat[iend,3,2], digits=4)))",marker=:rect,strokewidth=1.1)
            axislegend(ax, position=:rt)
            converged && ylims!(ax, mat[iend,1,2]-1, mat[iend,1,2]+2)

            ax = Axis(fig[1:2,2], xlabel=L"Iterations", ylabel=L"Fraction\;of\;TO_{ex}", rightspinecolor=:green, xgridvisible=false, ygridvisible=false, ylabelsize=25, xlabelsize=25); 
            lines!(ax, 1:iend, mat[1:iend,1,3],label="In solid Fe ($(round(x[1], digits=3)) TOₑₓ)",color=solidclr,linewidth=2.0)
            lines!(ax, 1:iend, mat[1:iend,2,3], color=meltclr,label="In melt Fe ($(round(x[2], digits=3)) TOₑₓ)",linewidth=2.0)
            lines!(ax, 1:iend, mat[1:iend,3,3], color=co2clr,label="In melt CO₂ ($(round(1 - x[1] - x[2], digits=3)) TOₑₓ)",linewidth=2.0)
            scatter!(ax, 1, mat[1,3,3], label="Melt XCO₂ = $(round(x[3], digits=3))", alpha=0.0)
            axislegend(ax, position=:rt)
            # Mark ceilings
                scatterlines!(ax, [1, iend], [maxsOₑₓ, maxsOₑₓ], alpha=0.3, color=solidclr,marker=:rect,strokewidth=1.1); text!(ax, 0.1iend, 1.02maxsOₑₓ, text="Solid Fe cap = $(round(maxsOₑₓ, digits=3))", fontsize=12, font=:italic, color=solidclr)
                scatterlines!(ax, [1, iend], [maxmOₑₓ, maxmOₑₓ], alpha=0.3, color=meltclr,marker=:rect,strokewidth=1.1); text!(ax, 0.4iend, 1.02maxmOₑₓ, text="Melt Fe cap = $(round(maxmOₑₓ, digits=3))", fontsize=12, font=:italic, color=meltclr)
                ax2 = Axis(fig[1:2,2], ylabel=L"XCO_2", yaxisposition=:right, ylabelcolor=co2clr, ytickcolor=co2clr, yticklabelcolor=co2clr, xgridvisible=false, ygridvisible=false); hidespines!(ax2, :l, :t, :b, :r); hidexdecorations!(ax2)
                scatterlines!(ax2, [1, iend], [maxXCO₂, maxXCO₂], alpha=0.3, color=co2clr,marker=:rect,strokewidth=1.1); text!(ax, 0.7iend, 1.02maxXCO₂, text="XCO₂ cap = $(round(maxXCO₂, digits=3))", fontsize=12, font=:italic, color=co2clr)
            # Limits
                ul = max(maxsOₑₓ, maxmOₑₓ, maxXCO₂)
                ylims!(ax, -0.05, 1.3ul); 
                ylims!(ax2, -0.05, 1.3ul);
            display(fig)
        end
        return x[1], x[2], x[3]

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

function Hirsch(T, mOₑₓ, T₀, ΔCₚ, a, b, c, y1, y3, y4, y5, y8, y9, IDV, SymXox, dummy, idxO, _ln10, _T, s, Φₘ, molXB)
    # Checks
    @assert (:O∈SymXox && :FeO∈SymXox) "This function requires FeO + O format!"
    # Reset dummy
    dummy .= molXB
    # Extract mass fraction of Oₑₓ in melt
    mfOₑₓ = Φₘ*mOₑₓ
    # Oxidize
    dummy[idxO] = mfOₑₓ
    # Normalize
    dummy./=(s+mfOₑₓ)
    # Construct bulk and assess hardlimit
    Xl = Cbulk((; zip(SymXox, dummy)...))
    (Xl.O>=0.5Xl.FeO) && (return 14) # Too much oxygen!! Above hard-limit.
    # Compute mfO2
    mfO2 = (log10(Xl.O/(Xl.FeO-2Xl.O)) - b - c*_T + (ΔCₚ/R*_ln10 * (1 - T₀*_T - log(T/T₀))) + IDV/(1e-3R)*_T*_ln10 
                        - _T*(y1*Xl.SiO2 + y3*Xl.MgO + y4*Xl.CaO + y5*Xl.Na2O + y8*Xl.SiO2*Xl.Al2O3 + y9*Xl.SiO2*Xl.MgO))/a
    return mfO2
end

# Stagno and Frost XCO2 equilibrium (2 Oex per 1 CO₂) | cOₑₓ in mass fraction of TOₑₓ
function XCO₂_to_fO2(XCO₂, Pin, Tin)
    # Checks
    @assert XCO₂>=0.0 "XCO₂ cannot be zero."
    # Limit to rexplored ranges
    P = Pin >= 11. ? 11. : Pin
    T = Tin >= c2k(1600.) ? c2k(1600.) : Tin
    # Compute logfO₂
    return 5.44 - 21380/T + 0.078(1e5P-1)/T + log10(XCO₂) - 12
end

function Rx(sfO2, P, T, sOₑₓ, mOₑₓ, XCO₂, T₀, ΔCₚ, a, b, c, y1, y3, y4, y5, y8, y9, IDV, SymXox, dummy, idxO, _ln10, _T, s, Φ, sw, Φₘ, cα, molXB)
    Rₛ = sfO2(sOₑₓ)
    return SA[Rₛ - Hirsch(T, mOₑₓ, T₀, ΔCₚ, a, b, c, y1, y3, y4, y5, y8, y9, IDV, SymXox, dummy, idxO, _ln10, _T, s, Φₘ, molXB)
              Rₛ - XCO₂_to_fO2(XCO₂, P, T)
              1 - sOₑₓ - mOₑₓ - Φ*cα/(sw+cα)]
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

# Variables carbon
@inline evΦ(TOₑₓ, Mf)       = Mf/TOₑₓ               # Conversion factor for XCO₂
@inline evcα(XCO₂)          = 2XCO₂*mm.O            # Molar XCO₂ in melt → non-normalized mass of Oₑₓ in melt stored in carbon
@inline evθ(cα, sw)         = sw + cα               # Molar normalization factor 
# Variables melt
@inline evΦₘ(TOₑₓ, Mf)      = TOₑₓ/Mf/mm.O          # Conversion factor for melt
@inline evα(mOₑₓ, Φₘ)       = Φₘ*mOₑₓ               # Mass fraction of TOₑₓ in melt → non-normalized mass of Oₑₓ in melt
@inline evθₘ(α, s)          = s - α                 # Molar normalization factor 

# Partial derivatives
@inline ∂S∂sOₑₓ(sfO2, sOlist) = @views diff(sfO2)./diff(sOlist)
@inline ∂C∂XCO₂(XCO₂, _ln10) = _ln10 / XCO₂
@inline ∂M∂mOₑₓ(Φₘ, Ys1, Ys2, α, θₘ, _ln10, a, XFeO) = -Φₘ/a * ( θₘ^(-2)*(Ys1 + 2Ys2/θₘ) - _ln10*(1/α + 2/(XFeO - 2α)))
@inline ∂3∂XCO₂(Φ, θ, cα) = 2Φ*mm.O*(1/θ - cα/θ^2)
