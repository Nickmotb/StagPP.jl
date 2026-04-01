# Model:
# Assess partitioning of TOₑₓ into solid and molten tracer + oxidation reduction of parcel carbon through equilibrium XCO₂
# Starts as a root-find solve in 3 dimensions. If XCO₂ touches a boundary, the problems loses an independent variable, and a 2D root-finding solve insues.
#
# Inputs:
# - Pressure (P)
# - Temperature (T)
# - Total oxygen budget either as wt% of solid (TOₑₓ) or solid and melt Fe³⁺/Feᵀ ratios (Rs, Rf)
# - Final melt fraction (ϕ)
# - Available reduced carbon in the parcel as wt% of solid tracer (TC)
#
# =================================
# ========= 3D root-find ==========
# =================================
#
# Independent variables:
# - Equilibrium mass fraction of TOₑₓ in the solid (sOₑₓ)
# - Equilibrium mass fraction of TOₑₓ in the melt  (mOₑₓ)
# - Equilibrium molar XCO₂ in the melt (XCO₂)
#
# Equations:
# (1) Equilibrium constraint    :  solid fO₂(P,T,sOₑₓ) = melt fO₂(P,T,mOₑₓ)
# (2) Equilibrium constraint    :  solid fO₂(P,T,sOₑₓ) = EDDOG fO₂(P,T,XCO₂) <-- Maybe exclude residual if XCO₂<0.0 or XCO₂>1.0
# (3a) Mass conservation         :  1 = sOₑₓ + mfOₑₓ + Oₑₓ_in_CO₂(XCO₂)  →  1 = sOₑₓ + mfOₑₓ + Φ*cα/(sw+cα)
#
# Jacobian : [∂(1)∂sOₑₓ ∂(1)∂mOₑₓ ∂(1)∂XCO₂         [∂S∂sOₑₓ     -∂M∂mOₑₓ        0
#             ∂(2)∂sOₑₓ ∂(2)∂mOₑₓ ∂(2)∂XCO₂    =     ∂S∂sOₑₓ        0       -∂C∂XCO₂
#             ∂(3)∂sOₑₓ ∂(3)∂mOₑₓ ∂(3)∂XCO₂]           -1           -1      -∂(Φ*XCO₂/(1-XCO₂))∂XCO₂]
#
# Variable extentions:
# cα         = 2XCO₂*mm.O  # non-normalized mass fraction Oₑₓ in melt
# Φ          = 2*mm.O*molMf/TOₑₓ
# sw, s      = ∑(oxᵢ*mmᵢ), ∑(oxᵢ/mmᵢ)
# Oₑₓ_in_CO₂ = Φ*cα/(sw+cα)
#
# Tools:
# - Hirschmann 2022 melt mapping from XFe₂O₃ (Oₑₓ) ↔ fO₂
# - Stixrude and Bertelloni 2024 (MAGEMin) solid mapping from XFe₂O₃ (Oₑₓ) ↔ fO₂
# - Stagno and Frost 2010 parameterization of melt EDDOG2 buffer fO₂ ↔ XCO₂
function partition_Oₑₓ(P::K, T::K, p::K, ϕ::K, TOex::K, TC::K; Rs::K=-1.0, Rf::K=-1.0, nr=50, niter=100, 
                        verbose=false, data=nothing, respace=(false, 20), plotevo=false) where {K <: Real}

    # Endmember bulks
    XH      = @SVector [0.4347, 0.4597, 0.0835, 0.0090, 0.0100, 0.0001, 0.0030, 0.0] # mass fraction
    XB      = @SVector [0.5097, 0.0988, 0.0718, 0.1268, 0.1698, 0.0225, 0.0007, 0.0] # mass fraction
    Xox     = @SVector ["SiO2", "MgO", "FeO", "CaO", "Al2O3", "Na2O", "Cr2O3", "O"]
    SymXox  = @SVector [:SiO2, :MgO, :FeO, :CaO, :Al2O3, :Na2O, :Cr2O3, :O]
    mmXox   = @SVector [mm.SiO2, mm.MgO, mm.FeO, mm.CaO, mm.Al2O3, mm.Na2O, mm.Cr2O3, mm.O]
    X       = p*XB + (1-p)*XH
    molXB   = XB./mmXox; molXB = molXB ./ sum(molXB)
    molX    = X./mmXox;  molX  = molX ./sum(molX)

    # Hirschmann 2022 parameters
    a=0.1917; b=-1.961; c=4158.1; ΔCₚ=33.25; T₀=1673.15; y1=-520.46; y2=-185.37; y3=494.39; y4=1838.34; y5=2888.48; y8=-1245.09; y9=-1156.86

    # Solid / Molten tracer mass [kg]
    Mt = 1e3 # 1 kg
    Mf=ϕ*Mt; Ms=Mt*(1-ϕ); Mc=TC*Mt
    molMf = sum(Mf.*molXB)   # Mf in mols (unoxidized)

    # Carbon constraints
    molCav      = Mc/mm.C                # Mols of available carbon
    maxXCO₂_raw = molCav/(molMf+molCav)     # Maximal XCO₂ allowed by carbon

    # IDV = ∫ΔVdP for melts
    IDV = solve_∫ΔVdP([P-0.05P, P, P+0.05P],[T-0.05T, T, T+0.05T])[2,2,1]

    # Flag whether to manage MAGEMin initialisation and finalization
    flag = isnothing(data)

    # Index of oxygen component
    idxO = findfirst(Xox.=="O")
    
    # === Allocate iterative memory and create bulk structures
        sol        = @SVector zeros(3)                                                  # Solution Vector (Raw)
        dummy      = zeros(length(XB))                                                  # for Hirschmann calls
        Xdummy     = (X=Vector{Float64}(X), Xox=Vector{String}(Xox), mm=get_Xoxmm(Xox)) # For oxidizing calls
        if Rs>=0.0 || Rf>=0.0
            (Rs<0.0) && (Rs=0.0); (Rf<0.0) && (Rf=0.0)                                                      # If an Fe³⁺ ratio is passed, set the inputs up
            Xs = oxidize_bulk(X, Rs, Xdummy, FeFormat="FeO_O", wt_out=true, frac=true, SymXox=SymXox);      # Oxidized solid bulk comp
            Xmo = oxidize_bulk(XB, Rf, Xdummy, FeFormat="FeO_O", wt_out=true, frac=true, SymXox=SymXox);    # Oxidized melt bulk comp
        end
        Xm = oxidize_bulk(XB, 0.0, Xdummy, FeFormat="FeO_O", wt_out=true, frac=true);                       # Unoxidized melt bulk comp

    # Sum TOₑₓ contributions from solid and molten tracer Oₑₓ (TOex -> input total mass O budget, TOₑₓ -> system total mass Oₑₓ budget)
    if Rs>=0.0 || Rf>=0.0
        verb_flag = -1; TOₑₓ = Xs.O*Ms + Xmo.O*Mf # If ratios passed, this considers both compositions
    else
        verb_flag = 1; TOₑₓ = TOex*Mt # If TOex passed, a fraction of the total mass becomes the Oₑₓ budget
    end
    if TOₑₓ==0.0 && verbose
        println("---- Solution (P=$(P)GPa | T=$(T)K | ϕ=$(ϕ) | Rs=$(Rs) | Rf=$(Rf) | Source mix = $p) ----")
        println("Total Oₑₓ budget = 0.0 kg")
        println("TOₑₓ partitioning → [0.0% solid, 0.0% melt]")
        return
    end

    # Compute partial derivative variables
    s           = sum(molXB)                                                    # ∑oxᵢ in mols
    Φ, Φₘ       = evΦ(TOₑₓ, molMf), evΦₘ(TOₑₓ, Mf)                              # Conversion factors       
    _ln10, _T   = 1/log(10), 1/T                                                # Pre computed constants
    Ys1         = (y1*molXB[1] + y3*molXB[2] + y4*molXB[4] + y5*molXB[6])*_T    # Sum of linear parameterized molar components
    Ys2         = molXB[1]*(y8*molXB[5] + y9*molXB[2])*_T                       # Sum of non-linear parameterized molar components

    # === Generate solid fO₂ space
        Rlist = LinRange(0.00001, 0.20, nr)
        Xlist = Vector{Vector{Float64}}(undef, nr); 
        sOₑₓlist = zeros(nr)
        for i in 1:nr
            Xl = oxidize_bulk(X, Rlist[i], Xdummy, wt_out=true, frac=true, FeFormat="FeO_O", SymXox=SymXox)
            Xlist[i] = [getfield(Xl, f) for f in SymXox]
            sOₑₓlist[i] = Xl.O
        end
        # -- Define maximal sOₑₓ
            maxsOₑₓ_uncapped  = (0.5molX[3]*mm.O)/(sum(X) - 0.5molX[3]*mm.O)*(Ms/TOₑₓ)
            maxsOₑₓ = min(maxsOₑₓ_uncapped, 1.0)
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
        
        # -- Barrier margin
            margin = 1e-9
        # -- Compute independent boundaries
            maxmOₑₓ_uncapped  = (0.5molXB[3]*mm.O)/(sum(XB) + 0.5molXB[3]*mm.O)*(Mf/TOₑₓ)
            maxmOₑₓ, maxXCO₂  = min(maxmOₑₓ_uncapped, 1.0), max(min(1.0, maxXCO₂_raw), 0.0)
        # -- Initialise (static)solution vector
            y = sol + ([x_to_y(0.3maxsOₑₓ, maxsOₑₓ, margin), x_to_y(0.01maxmOₑₓ, maxmOₑₓ, margin), x_to_y(0.01maxXCO₂, maxXCO₂, margin)])
        # -- Define convergence tolerance (ϵ), correction dampening factor (damp), and maximum remaining residual (aR)
            ϵ, damp = 1e-2, 0.15
        # -- Wrap parameters and call solver
            params = (; verb_flag, P, T, ϕ, Rs, Rf, TOex, TOₑₓ, p, TC, Φ, s, Φₘ,
                            T₀, ΔCₚ, a, b, c, y1, y3, y4, y5, y8, y9, IDV, SymXox, dummy, idxO, _ln10, _T, molXB,
                                Ys1, Ys2, plotevo, verbose, margin, Mt)
            sOₑₓ, mOₑₓ, XCO₂, converged, mat = constrained_smOₑₓ_XCO₂_solver(y, maxsOₑₓ, maxmOₑₓ, maxXCO₂, sfO2, ∂Sᵢ, ϵ, damp, niter; params...)
        # -- Plot evolution if requested
            if plotevo
                x = [sOₑₓ, mOₑₓ, XCO₂]
                solidclr = :black
                meltclr  = :red
                co2clr   = :green
                # Plot evolution
                iend = findfirst(mat[:,1,1].==0.0); isnothing(iend) ? (iend=niter) : (iend-=1)
                fig = Figure(size=(1800, 800))
                ax = Axis(fig[1,1], ylabel=L"Solution\;residual", xlabel=L"Iterations", xgridvisible=false, ygridvisible=false, ylabelsize=25, xlabelsize=25); 
                scatterlines!(ax, 1:iend, mat[1:iend,1,1], color=meltclr, label="(eq. 1) Solid ↔ melt fO₂ equilibrium (R₁ = $(round(mat[iend,1,1], digits=4)))",marker=:rect,strokewidth=1.1)
                scatterlines!(ax, 1:iend, mat[1:iend,2,1], color=co2clr,label="(eq. 2) Solid ↔ EDDOG fO₂ equilibrium (R₂ = $(round(mat[iend,2,1], digits=4)))",marker=:rect,strokewidth=1.1)
                scatterlines!(ax, 1:iend, mat[1:iend,3,1], color=:orange,label="(eq. 3) Mass conservation (R₃ = $(round(mat[iend,3,1], digits=4)))",marker=:rect,strokewidth=1.1)
                lines!(ax, [1, iend], [ϵ, ϵ], linestyle=:dash, color=:gray, alpha=0.5)
                lines!(ax, [1, iend], [-ϵ, -ϵ], linestyle=:dash, color=:gray, alpha=0.5)
                text!(ax, 1, -ϵ-0.1, text="Convergence ϵ = $ϵ", font=:italic, alpha=0.5)
                axislegend(ax, position=:rt)
                converged && ylims!(ax, -0.25, 0.5)

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
                scatter!(ax, 1, mat[1,3,3], label="Melt XCO₂ = $(round(x[3], digits=5))", alpha=0.0)
                axislegend(ax, position=:rt)
                # Mark ceilings
                    scatterlines!(ax, [1, iend], [maxsOₑₓ, maxsOₑₓ], alpha=0.3, color=solidclr,marker=:rect,strokewidth=1.1); text!(ax, 0.1iend, 1.02maxsOₑₓ, text="Solid Fe³⁺ cap = $(round(maxsOₑₓ, digits=3))"*(maxsOₑₓ_uncapped>1.0 ? " ($(round(maxsOₑₓ_uncapped, digits=3)))" : ""), fontsize=12, font=:italic, color=solidclr)
                    scatterlines!(ax, [1, iend], [maxmOₑₓ, maxmOₑₓ], alpha=0.3, color=meltclr,marker=:rect,strokewidth=1.1); text!(ax, 0.4iend, 1.02maxmOₑₓ, text="Melt Fe³⁺ cap = $(round(maxmOₑₓ, digits=3))"*(maxmOₑₓ_uncapped>1.0 ? " ($(round(maxmOₑₓ_uncapped, digits=3)))" : ""), fontsize=12, font=:italic, color=meltclr)
                    ax2 = Axis(fig[1:2,2], ylabel=L"XCO_2", yaxisposition=:right, ylabelcolor=co2clr, ytickcolor=co2clr, yticklabelcolor=co2clr, xgridvisible=false, ygridvisible=false); hidespines!(ax2, :l, :t, :b, :r); hidexdecorations!(ax2)
                    scatterlines!(ax2, [1, iend], [maxXCO₂, maxXCO₂], alpha=0.3, color=co2clr,marker=:rect,strokewidth=1.1); text!(ax, 0.7iend, 1.02maxXCO₂, text="XCO₂ cap = $(round(maxXCO₂, digits=3))", fontsize=12, font=:italic, color=co2clr)
                # Limits
                    ul = max(maxsOₑₓ, maxmOₑₓ, maxXCO₂)
                    ylims!(ax, 0.0, 1.3ul); 
                    ylims!(ax2, 0.0, 1.3ul);
                display(fig)
            end
        # Return partitioning
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

function Hirsch(T, mOₑₓ, T₀, ΔCₚ, a, b, c, y1, y3, y4, y5, y8, y9, IDV, SymXox, dummy, idxO, _ln10, _T, s, Φₘ, molXB)
    # Checks
    @assert (:O∈SymXox && :FeO∈SymXox) "This function requires FeO + O format!"
    # Reset dummy
    dummy .= molXB
    # Compute non-normalized mols of Oₑₓ in melt
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
function XCO₂_to_fO2(XCO₂, P, T)
    # Checks
    @assert XCO₂>=0.0 "XCO₂ cannot be zero."
    # Limit to rexplored ranges
    # P = Pin >= 11. ? 11. : Pin
    # T = Tin >= c2k(1600.) ? c2k(1600.) : Tin
    # Compute logfO₂
    return 5.44 - 21380/T + 0.078(1e5P-1)/T + log10(XCO₂) - 10
end

function Rx(sfO2, P, T, sOₑₓ, mOₑₓ, XCO₂, T₀, ΔCₚ, a, b, c, y1, y3, y4, y5, y8, y9, IDV, SymXox, dummy, idxO, _ln10, _T, s, Φ, Φₘ, molXB, D)
    Rₛ = sfO2(sOₑₓ)
    if D==:D3
        return SA[Rₛ - Hirsch(T, mOₑₓ, T₀, ΔCₚ, a, b, c, y1, y3, y4, y5, y8, y9, IDV, SymXox, dummy, idxO, _ln10, _T, s, Φₘ, molXB)
                Rₛ - XCO₂_to_fO2(XCO₂, P, T)
                1 - sOₑₓ - mOₑₓ - Φ*XCO₂/(1-XCO₂)]
    elseif D==:D2
        return SA[Rₛ - Hirsch(T, mOₑₓ, T₀, ΔCₚ, a, b, c, y1, y3, y4, y5, y8, y9, IDV, SymXox, dummy, idxO, _ln10, _T, s, Φₘ, molXB)
                1 - sOₑₓ - mOₑₓ - Φ*XCO₂/(1-XCO₂)]
    end
end

function constrained_smOₑₓ_XCO₂_solver(y    :: SVector{3,Float64},             # Initial transformed solution vector
                                       slim :: K, mlim :: K, clim :: K,        # Limits for smOₑₓ, mOₑₓ and XCO₂
                                       sfO2, ∂Sᵢ,                              # Solid fO₂ and numerical derivative (MAGEMin)
                                       ϵ    :: K, damp :: K, niter :: Int64;   # Solver tolerance, dampening factor, and maximum iterations
                                       # Verbose parameters
                                       verb_flag::Int64,P::K,T::K,ϕ::K,Rs::K,Rf::K,
                                       TOex::K,TOₑₓ::K,p::K,TC::K,Mt::K,
                                       # Pre-computed parameters
                                       s::K,Φ::K,Φₘ::K,T₀::K,ΔCₚ::K,a::K,b::K,c::K,y1::K,y3::K,margin::K,
                                       y4::K,y5::K,y8::K,y9::K,IDV::K,_ln10::K,_T::K,Ys1::K,Ys2::K,
                                       idxO::Int64, plotevo::Bool,SymXox::SVector{N, Symbol},verbose::Bool,
                                       dummy::Vector{Float64}, molXB::SVector{N, Float64}) where{K<:AbstractFloat, N}
    
    # Residual - fO₂ - Partitioning matrix for plotting
    mat = zeros(niter, 3, 3) # [Residuals, fO₂, Partitioning]

    # y margins
    y_lowc = x_to_y(0.0, clim, margin)
    y_highc = x_to_y(clim, clim, margin)

    # Converged flag
    converged = false

    D       = :D3                 # Start as a 3D solver
    SM3     = @SMatrix zeros(3,3) # Static Matrix 3 × 3
    SM2     = @SMatrix zeros(2,2) # Static Matrix 3 × 3
    for it in 1:niter
        # Evaluate current stage
        y₁, y₂, y₃ = y
        # Transform back to original variables
        sOₑₓ, mOₑₓ, XCO₂ = y_to_x(y₁, slim, margin), y_to_x(y₂, mlim, margin), y_to_x(y₃, clim, margin)
        println("Iteration $it: sOₑₓ = $(round(sOₑₓ, digits=4)) ($(round(slim, digits=4))), mOₑₓ = $(round(mOₑₓ, digits=4)) ($(round(mlim, digits=4))), XCO₂ = $(round(XCO₂, digits=4)) ($(round(clim, digits=4)))")
        # Iteration variables
        α   = evα(mOₑₓ, Φₘ)
        θₘ  = evθₘ(α, s)
        # Compute residual
        Fx = Rx(sfO2, P, T, sOₑₓ, mOₑₓ, XCO₂, T₀, ΔCₚ, a, b, c, y1, y3, y4, y5, y8, y9, IDV, SymXox, dummy, idxO, _ln10, _T, s, Φ, Φₘ, molXB, D)
        D==:D3 && println("             R₁ = $(round(Fx[1], digits=5)),     R₂ = $(round(Fx[2], digits=5)),     R₃ = $(round(Fx[3], digits=5))")
        D==:D2 && println("             R₁ = $(round(Fx[1], digits=5)),     R₃ = $(round(Fx[2], digits=5))")
        aR = maximum(abs.(Fx))
        # Store values for plotting
        if plotevo
            D==:D3 ? (mat[it,:,1] .= Fx) : (mat[it,1,1] = Fx[1]; mat[it,2,1] = 0.0; mat[it,3,1] = Fx[2])
            mat[it,1,2]  = sfO2(sOₑₓ)
            mat[it,2,2]  = Hirsch(T, mOₑₓ, T₀, ΔCₚ, a, b, c, y1, y3, y4, y5, y8, y9, IDV, SymXox, dummy, idxO, _ln10, _T, s, Φₘ, molXB)
            mat[it,3,2]  = XCO₂_to_fO2(XCO₂, P, T)
            mat[it,:,3] .= [sOₑₓ, mOₑₓ, XCO₂]
        end
        # Check convergence
        if aR<=ϵ
            converged = true
            if verbose
                if verb_flag==-1
                    println("---- Solution (P=$(P)GPa | T=$(T)K | ϕ=$(ϕ) | Rs=$(Rs) | Rf=$(Rf) | Mix=$p | TCarbon=$TC) ----")
                else
                    println("---- Solution (P=$(P)GPa | T=$(T)K | ϕ=$(ϕ) | TOₑₓ=$TOex | Mix=$p | TCarbon=$TC) ----")
                end
                println("Shared fO₂ = $(sfO2(sOₑₓ)) |  residual = $aR")
                println("Total Oₑₓ budget = $(round((1e2TOₑₓ/Mt), digits=4))% of total mass")
                println("TOₑₓ partitioning → [$(round(1e2sOₑₓ, digits=4))% solid, $(round(1e2mOₑₓ, digits=4))% melt] stored as Fe₂O₃, $(round(1e2(1 - sOₑₓ - mOₑₓ), digits=4))% stored as melt CO₂")
                println("Melt XCO₂ = $(XCO₂)")
                println("Converged in $it iterations.")
            end
            break
        end
        # Partial derivatives
        ∂S = ∂Sᵢ(sOₑₓ)
        ∂M = ∂M∂mOₑₓ(Φₘ, Ys1, Ys2, α, θₘ, _ln10, a, molXB[3])
        ∂C = ∂C∂XCO₂(XCO₂, _ln10)
        ∂3 = ∂3∂XCO₂(Φ, XCO₂)
        # Jacobian inverse (Chain rule)
        if D==:D3
            ∂x∂y₁, ∂x∂y₂, ∂x∂y₃ = ∂x∂y(y₁, slim, margin), ∂x∂y(y₂, mlim, margin), ∂x∂y(y₃, clim, margin)
            J⁻¹  = SM3 + inv([ ∂S*∂x∂y₁ -∂M*∂x∂y₂ 0.0; ∂S*∂x∂y₁ 0.0 -∂C*∂x∂y₃; -∂x∂y₁ -∂x∂y₂ -∂3*∂x∂y₃]) # J = ∂Rᵢ∂yᵢ =  ∂Rᵢ∂xᵢ * ∂xᵢ∂yᵢ
        elseif D==:D2
            ∂x∂y₁, ∂x∂y₂ = ∂x∂y(y₁, slim, margin), ∂x∂y(y₂, mlim, margin)
            J⁻¹  = SM2 + inv([ ∂S*∂x∂y₁ -∂M*∂x∂y₂; -∂x∂y₁ -∂x∂y₂]) # J = ∂Rᵢ∂yᵢ =  ∂Rᵢ∂xᵢ * ∂xᵢ∂yᵢ
        end
        # Newton step
        dn = J⁻¹*Fx*damp
        # Adaptive step
        if D==:D3
            α = 1.0
            (y[3]-dn[3]<=y_lowc)  && (α = (y[3]-y_lowc)/dn[3]; D=:D2; println("Switched to 2D solver at it=$it"); continue)
            (y[3]-dn[3]>=y_highc) && (α = (y[3]-y_highc)/dn[3]; D=:D2; println("Switched to 2D solver at it=$it"); continue)
        end
        # Take step
        D==:D3 ? (y = y - dn) : (y = y - [dn[1], dn[2], 0.0])

        # Output if not converged
        if verbose && it==niter
            if verb_flag==-1
                println("---- Solution (P=$(P)GPa | T=$(T)K | ϕ=$(ϕ) | Rs=$(Rs) | Rf=$(Rf) | Mix=$p | TCarbon=$TC) ) ----")
            else
                println("---- Solution (P=$(P)GPa | T=$(T)K | ϕ=$(ϕ) | TOₑₓ=$TOex | Mix=$p | TCarbon=$TC) ----")
            end
            println("Shared fO₂ = $(sfO2(sOₑₓ)) |  residual = $aR")
            println("Total Oₑₓ budget = $(round((1e2TOₑₓ/Mt), digits=4))% of total mass")
            println("TOₑₓ partitioning → [$(round(1e2sOₑₓ, digits=4))% solid, $(round(1e2mOₑₓ, digits=4))% melt] stored as Fe₂O₃, $(round(1e2(1 - sOₑₓ - mOₑₓ), digits=4))% stored as melt CO₂")
            println("Melt XCO₂ = $(XCO₂)")
            println("Converged in $it iterations.")
        end
    end

    sOₑₓ, mOₑₓ, XCO₂ = y_to_x(y[1], slim, margin), y_to_x(y[2], mlim, margin), y_to_x(y[3], clim, margin)
    return sOₑₓ, mOₑₓ, XCO₂, converged, mat
end

# Variable transformations
@inline y_to_x(y, maxV, margin)       = margin + (maxV - 2margin)/(1 + exp(-y))
function x_to_y(x, maxV, margin)
    safe = clamp(x, margin, maxV-margin)
    return -log((maxV-2margin)/(safe-margin) - 1.0)
end
@inline ∂x∂y(y, maxV, margin)    = (maxV-2margin)*exp(-y)/(1.0+exp(-y))^2
# Variables carbon
@inline evΦ(TOₑₓ, molMf) = 2*mm.O*molMf/TOₑₓ     # Conversion factor for XCO₂
# Variables melt
@inline evΦₘ(TOₑₓ, Mf)   = TOₑₓ/Mf/mm.O          # Conversion factor for melt
@inline evα(mOₑₓ, Φₘ)    = Φₘ*mOₑₓ               # Mass fraction of TOₑₓ in melt → non-normalized mass of Oₑₓ in melt
@inline evθₘ(α, s)       = s - α                 # Molar normalization factor
# Partial derivatives
@inline ∂S∂sOₑₓ(sfO2, sOlist)                        = @views diff(sfO2)./diff(sOlist)
@inline ∂C∂XCO₂(XCO₂, _ln10)                         = _ln10 / XCO₂
@inline ∂M∂mOₑₓ(Φₘ, Ys1, Ys2, α, θₘ, _ln10, a, XFeO) = -Φₘ/a * ( θₘ^(-2)*(Ys1 + 2Ys2/θₘ) - _ln10*(1/α + 2/(XFeO - 2α)))
@inline ∂3∂XCO₂(Φ, XCO₂)                             = Φ*(1/(1-XCO₂)^2)