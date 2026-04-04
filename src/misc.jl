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
# (1) Equilibrium constraint    :  solid fO₂(P,T,sOₑₓ) - melt fO₂(P,T,mOₑₓ) = 0
# (2) Equilibrium constraint    :  aᵪ[melt fO₂(P,T,mOₑₓ) - EDDOG fO₂(P,T,XCO₂)] = 0
# (3) Mass conservation         :  1 - sOₑₓ - mfOₑₓ - Φ*XCO₂/(1-XCO₂) - fₑₓ = 0
# (4) CaCO₃ oxidation regime    :  (clim - XCO₂)*fₑₓ
#
# Jacobian : [∂(1)∂sOₑₓ ∂(1)∂mOₑₓ ∂(1)∂XCO₂         [∂S∂sOₑₓ     -∂M∂mOₑₓ        0
#             ∂(2)∂sOₑₓ ∂(2)∂mOₑₓ ∂(2)∂XCO₂    =        0         ∂M∂mOₑₓ       -∂C∂XCO₂
#             ∂(3)∂sOₑₓ ∂(3)∂mOₑₓ ∂(3)∂XCO₂]           -1           -1      -∂(Φ*XCO₂/(1-XCO₂))∂XCO₂]
#
# Variable extentions:
# Φ          = 2*mm.O*molMf/TOₑₓ
# Φcac       = mmCaCO3*TOₑₓ/(3mmO)
# sw, s      = ∑(oxᵢ*mmᵢ), ∑(oxᵢ/mmᵢ)
# Mcaco3     = Φcac*fₑₓ
# aᵪ         = (1 - fₑₓ/(fₑₓ+eps))
#
# Tools:
# - Hirschmann 2022 melt mapping from XFe₂O₃ (Oₑₓ) ↔ fO₂
# - Stixrude and Bertelloni 2024 (MAGEMin) solid mapping from XFe₂O₃ (Oₑₓ) ↔ fO₂
# - Stagno and Frost 2010 parameterization of melt EDDOG2 buffer fO₂ ↔ XCO₂
function partition_Oₑₓ(P::K, T::K, p::K, ϕ::K, TOex::K, TC::K; Rs::K=-1.0, Rf::K=-1.0, nr=25, niter=100, 
                        verbose=false, data=nothing, Rspace=false, plotevo=false, damp=0.25, debugging=false) where {K <: Real}

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
    molCav      = Mc/mm.C                   # Mols of available carbon
    maxXCO₂_raw = molCav/(molMf+molCav)     # Maximal XCO₂ allowed by carbon

    # IDV = ∫ΔVdP for melts
    IDV = solve_∫ΔVdP([P-0.05P, P, P+0.05P],[T-0.05T, T, T+0.05T])[2,2,1]

    # Flag whether to manage MAGEMin initialisation and finalization
    flag = isnothing(data)

    # Index of oxygen component
    idxO = findfirst(Xox.=="O")
    
    # === Allocate iterative memory and create bulk structures
        sol        = @SVector zeros(4)                                                  # Solution Vector (Raw)
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

    # -- Solver boundary margin and XCO₂ sharpness parameter
    lowclip   = 1e-12
    sharpness = 0.0

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
        sfO2   = extrapolate(interpolate((sOₑₓlist.*Ms./TOₑₓ,), [out[i].fO2 for i in eachindex(out)], Gridded(Linear())), Flat())
        sfO2⁻¹ = 0.0 #extrapolate(interpolate(([out[i].fO2 for i in eachindex(out)],), sOₑₓlist.*Ms./TOₑₓ, Gridded(Linear())), Flat())
        sample_sOlist = LinRange(0.8lowclip, 1.2maxsOₑₓ, 100)
        sample_sOlist05 = 0.5(sample_sOlist[1:end-1] + sample_sOlist[2:end])
        sampled_sfO2 = sfO2(sample_sOlist)
    # -- Solid partial derivative
        ∂Sᵢ   = extrapolate(interpolate((sample_sOlist05,), ∂S∂sOₑₓ(sampled_sfO2, sample_sOlist), Gridded(Linear())), Line())
        
    # -- Compute independent boundaries
    maxmOₑₓ_uncapped  = (0.5molXB[3]*mm.O)/(sum(XB) + 0.5molXB[3]*mm.O)*(Mf/TOₑₓ)
    maxmOₑₓ, maxXCO₂  = min(maxmOₑₓ_uncapped, 1.0), max(min(1.0, maxXCO₂_raw), 0.0)
    maxfₑₓ            = 1.0 
    # -- Initialise (static)solution vector
    y = sol + [x_to_y(0.5maxsOₑₓ, lowclip, maxsOₑₓ), x_to_y(0.5maxmOₑₓ, lowclip, maxmOₑₓ), x_to_y(0.5maxXCO₂, lowclip, maxXCO₂), x_to_y(lowclip, lowclip, maxfₑₓ)]
    # -- Define convergence tolerance (ϵ)
    ϵ = 1e-2          
    # -- Wrap parameters and call solver
    params = (; verb_flag, P, T, ϕ, Rs, Rf, TOex, TOₑₓ, p, TC, Φ, s, Φₘ,
                    T₀, ΔCₚ, a, b, c, y1, y3, y4, y5, y8, y9, IDV, SymXox, dummy, idxO, _ln10, _T, molXB,
                        Ys1, Ys2, plotevo, verbose, lowclip, Mt, sharpness, debugging)
    converged, mat, itout = constrained_smOₑₓ_XCO₂_solver(y, maxsOₑₓ, maxmOₑₓ, maxXCO₂, maxfₑₓ, sfO2, sfO2⁻¹, ∂Sᵢ, ϵ, damp, niter; params...)
    # -- Plot evolution if requested
    if plotevo
        x = mat[itout, :, 3]
        solidclr = :black
        meltclr  = :red
        co2clr   = :green
        fₑₓclr   = :purple
        # Plot evolution
        fig = Figure(size=(1800, 800))
        ax = Axis(fig[1,1], ylabel=L"Solution\;residual", xlabel=L"Iterations", xgridvisible=false, ygridvisible=false, ylabelsize=25, xlabelsize=25); 
        scatterlines!(ax, 1:itout, mat[1:itout,1,1], color=meltclr, label="(eq. 1) Solid ↔ melt fO₂ equilibrium (R₁ = $(round(mat[itout,1,1], digits=4)))",marker=:rect,strokewidth=1.1)
        scatterlines!(ax, 1:itout, mat[1:itout,2,1], color=co2clr,label="(eq. 2) Solid ↔ EDDOG fO₂ equilibrium (R₂ = $(round(mat[itout,2,1], digits=4)))",marker=:rect,strokewidth=1.1)
        scatterlines!(ax, 1:itout, mat[1:itout,3,1], color=:orange,label="(eq. 3) Mass conservation (R₃ = $(round(mat[itout,3,1], digits=4)))",marker=:rect,strokewidth=1.1)
        lines!(ax, [1, itout], [ϵ, ϵ], linestyle=:dash, color=:gray, alpha=0.5)
        lines!(ax, [1, itout], [-ϵ, -ϵ], linestyle=:dash, color=:gray, alpha=0.5)
        text!(ax, 1, -1.9ϵ, text="Convergence ϵ = $(1e2ϵ)%", font=:italic, alpha=0.5)
        axislegend(ax, position=:rt)
        converged && ylims!(ax, -5ϵ, 7ϵ)

        ax = Axis(fig[2,1], ylabel=L"log\;fO_2", xlabel=L"Iterations", xgridvisible=false, ygridvisible=false, ylabelsize=25, xlabelsize=25); 
        scatterlines!(ax, 1:itout, mat[1:itout,1,2], color=solidclr, label="Solid (fO₂ = $(round(mat[itout,1,2], digits=4)))",marker=:rect,strokewidth=1.1)
        scatterlines!(ax, 1:itout, mat[1:itout,2,2], color=meltclr,label="Melt (fO₂ = $(round(mat[itout,2,2], digits=4)))",marker=:rect,strokewidth=1.1)
        scatterlines!(ax, 1:itout, mat[1:itout,3,2], color=co2clr,label="EDDOG (fO₂ = $(round(mat[itout,3,2], digits=4)))",marker=:rect,strokewidth=1.1)
        axislegend(ax, position=:rt)
        converged && ylims!(ax, mat[itout,3,2]-0.5, mat[itout,3,2]+1)

        ax = Axis(fig[1:2,2], xlabel=L"Iterations", ylabel=L"Fraction\;of\;TO_{ex}", rightspinecolor=:green, xgridvisible=false, ygridvisible=false, ylabelsize=25, xlabelsize=25); 
        lines!(ax, 1:itout, mat[1:itout,1,3],label="In solid Fe ($(round(x[1], digits=3)) TOₑₓ)",color=solidclr,linewidth=2.0)
        lines!(ax, 1:itout, mat[1:itout,2,3], color=meltclr,label="In melt Fe ($(round(x[2], digits=3)) TOₑₓ)",linewidth=2.0)
        lines!(ax, 1:itout, mat[1:itout,4,3], color=fₑₓclr,label="Precipitated ($(round(x[4], digits=3)) TOₑₓ)",linewidth=2.0)
        scatter!(ax, 1, mat[1,3,3], label="Melt XCO₂ = $(round(x[3], digits=5))", alpha=0.0)
        axislegend(ax, position=:rt, framevisible=true)
        # Mark ceilings
        scatterlines!(ax, [1, itout], [maxsOₑₓ, maxsOₑₓ], alpha=0.3, color=solidclr,marker=:rect,strokewidth=1.1); text!(ax, 0.1itout, 1.02maxsOₑₓ, text="Solid Fe³⁺ cap = $(round(maxsOₑₓ, digits=3))"*(maxsOₑₓ_uncapped>1.0 ? " ($(round(maxsOₑₓ_uncapped, digits=3)))" : ""), fontsize=12, font=:italic, color=solidclr)
        scatterlines!(ax, [1, itout], [maxmOₑₓ, maxmOₑₓ], alpha=0.3, color=meltclr,marker=:rect,strokewidth=1.1); text!(ax, 0.4itout, 1.02maxmOₑₓ, text="Melt Fe³⁺ cap = $(round(maxmOₑₓ, digits=3))"*(maxmOₑₓ_uncapped>1.0 ? " ($(round(maxmOₑₓ_uncapped, digits=3)))" : ""), fontsize=12, font=:italic, color=meltclr)
        ax2 = Axis(fig[1:2,2], ylabel=L"XCO_2", yaxisposition=:right, ylabelcolor=co2clr, ytickcolor=co2clr, yticklabelcolor=co2clr, xgridvisible=false, ygridvisible=false); hidespines!(ax2, :l, :t, :b, :r); hidexdecorations!(ax2)
            lines!(ax2, 1:itout, mat[1:itout,3,3], color=co2clr,label="In melt CO₂ ($(round(1 - x[1] - x[2], digits=3)) TOₑₓ)",linewidth=2.0)
            scatterlines!(ax2, [1, itout], [maxXCO₂, maxXCO₂], alpha=0.3, color=co2clr,marker=:rect,strokewidth=1.1); text!(ax2, 0.7itout, 1.1maxXCO₂, text="Carbon limited XCO₂ = $(round(maxXCO₂, digits=3))", fontsize=12, font=:italic, color=co2clr)
        # Limits
            ul = max(maxsOₑₓ, maxmOₑₓ, maxXCO₂)
            ylims!(ax, 0.0, 1.07ul); 
            ylims!(ax2, 0.0, 1.07ul);
        display(fig)
    end

    if Rspace
        # Return partitioning
        return x[1], x[2], x[3]
    else
        return converged, mat
    end

end

function P_T_ϕ_TOₑₓ_TC_topoplogy(data, ivars; p::Float64=0.2, ns=10)
    s = split(ivars, "_")
    @assert length(s)==2 "The function can only accept two independent variables!"
    y, x = s[1], s[2]
    # Resolutions
    defP, defT, defϕ, defTO, defTC = 3.0, 1600., 0.02, 0.01, 0.05
    # Construct map
    P  = contains(ivars, "P")  ? LinRange(1.0, 4.0, ns)      : defP
    T  = contains(ivars, "T")  ? LinRange(1100., 1800., ns)  : defT
    TO = contains(ivars, "TO") ? LinRange(0.05, 0.15, ns)    : defTO
    TC = contains(ivars, "TC") ? LinRange(0.01, 0.1, ns)     : defTC
    ϕ  = contains(ivars, "ϕ")  ? LinRange(0.05, 0.1, ns)     : defϕ
    # Set vecs
    yvec = y=="P" ? P : y=="T" ? T : y=="ϕ" ? ϕ : y=="TO" ? TO : TC
    xvec = x=="P" ? P : x=="T" ? T : x=="ϕ" ? ϕ : x=="TO" ? TO : TC
    # it
    mapYX = zeros(ns,ns,5)
    it, nmax = 0, ns*ns
    for n1 in 1:ns
        for n2 in 1:ns
            it+=1
            println("n = $it / $nmax")
            converged, res =  partition_Oₑₓ( isa(P,Float64)   ? P  :  P[y=="P" ? n1 : n2],
                                            isa(T,Float64)   ? T  :  T[y=="T" ? n1 : n2],
                                            p,
                                            isa(ϕ,Float64)   ? ϕ  :  ϕ[y=="ϕ" ? n1 : n2],
                                            isa(TO,Float64)  ? TO : TO[y=="TO" ? n1 : n2],
                                            isa(TC,Float64)  ? TC : TC[y=="TC" ? n1 : n2], data=data, verbose=false, plotevo=false)
            mapYX[n1,n2,1:4] .= res
            mapYX[n1,n2,5]    = converged ? 1 : 0
        end
    end
    # Set asbolutes
    mapYX.=abs.(mapYX)

    fig = Figure(size=(1200, 1000))
    ax = Axis3(fig[1,1], xlabel=x, ylabel=y, zlabel="R₁ (fO₂ₛₘ)")
        surface!(ax, xvec, yvec, mapYX[:,:,1]', colormap=:Purples, alpha=1.0)
    ax = Axis3(fig[1,2], xlabel=x, ylabel=y, zlabel="R₂ (fO₂ₘᵪ)")
        surface!(ax, xvec, yvec, mapYX[:,:,2]', colormap=:Purples, alpha=1.0)
    ax = Axis3(fig[2,1], xlabel=x, ylabel=y, zlabel="R₃ (mass)")
        surface!(ax, xvec, yvec, mapYX[:,:,3]', colormap=:Purples, alpha=1.0)
    ax = Axis3(fig[2,2], xlabel=x, ylabel=y, zlabel="R₄ (C → CO₃)")
        surface!(ax, xvec, yvec, mapYX[:,:,4]', colormap=:Purples, alpha=1.0)
        # wireframe!(ax, yvec, xvec, mapYX, color=:black)
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

function Hirsch⁻¹(T, eqfO2, T₀, ΔCₚ, a, b, c, y1, y3, y4, y5, y8, y9, IDV, SymXox, dummy, idxO, _ln10, _T, s, Φₘ, molXB, Ys1, Ys2, mlim)
    iguess = 0.05mlim; ϵ = 1e-2
    for it in 1:20
        R      = eqfO2 - Hirsch(T, iguess, T₀, ΔCₚ, a, b, c, y1, y3, y4, y5, y8, y9, IDV, SymXox, dummy, idxO, _ln10, _T, s, Φₘ, molXB)
        (abs(R)<=ϵ) && break
        α      = evα(iguess, Φₘ)
        θₘ     = evθₘ(α, s)
        ∂M     = -∂M∂mOₑₓ(Φₘ, Ys1, Ys2, α, θₘ, _ln10, a, molXB[3])
        iguess = iguess - R/∂M
    end
    return iguess
end

# Modified from Stagno and Frost XCO2 equilibrium (2 Oex per 1 CO₂) | cOₑₓ in mass fraction of TOₑₓ
function XCO₂_to_fO2(XCO₂, Pin, Tin, sharpness, clim)
    # Checks
    @assert XCO₂>=0.0 "XCO₂ cannot be zero."
    # Limit to rexplored ranges
    P = Pin > 11. ? 11. : Pin < 2.5 ? 2.5 : Pin
    T = Tin >= c2k(1600.) ? c2k(1600.) : Tin < c2k(1100.) ? c2k(1100.) : Tin
    # Compute logfO₂
    if sharpness==0.0
        return 5.44 - 21380/T + 0.078(1e4P-1)/T + log10(XCO₂)
    else
        return 5.44 - 21380/T + 0.078(1e4P-1)/T + log10(XCO₂) - sharpness*log10(clim-XCO₂)
    end
end

# From Ariskin and Barmina 1999
function magnetite_buffer_fO2(T, logfO2)
    # Ariskin and Barmina 1999 parameters
    # b0 = -0.6232; b1=0.6826; b2=4.0438; b3=0.8251; c=0.0159; d1=8.3626; d2=35.9674; d3=10.7347
    a=20964.7; b=0.6159; c=-8.8187; d1=-10.1838; d2=32.7525; d3=-14.2506
    # Parameterized melt composition
    SiO2=0.5109/mm.SiO2; TiO2=0.0207/mm.TiO2; Al2O3=0.1534/mm.Al2O3/2
    FeO=0.1330/mm.FeO; MgO=0.0617/mm.MgO; CaO=0.0896/mm.CaO; Na2O= 0.0199/mm.Na2O/2
    K2O=0.0109/mm.K2O/2
    # Normalize
    s = SiO2 + TiO2 + Al2O3 + FeO + MgO + CaO + Na2O + K2O
    SiO2/=s; TiO2/=s; Al2O3/=s; FeO/=s; MgO/=s; CaO/=s; Na2O/=s; K2O/=s
    # Compute logfO₂
    return exp(a/T + b*logfO2 + c +d1*Na2O + d2*K2O)
end

function Rx(sfO2, P, T, sOₑₓ, mOₑₓ, XCO₂, fₑₓ, T₀, ΔCₚ, a, b, c, y1, y3, y4, y5, y8, y9, IDV, SymXox, dummy, idxO, _ln10, _T, s, Φ, Φₘ, molXB, sharpness, clim, Dsat)
    fₛ = sfO2(sOₑₓ)
    fₗ = Hirsch(T, mOₑₓ, T₀, ΔCₚ, a, b, c, y1, y3, y4, y5, y8, y9, IDV, SymXox, dummy, idxO, _ln10, _T, s, Φₘ, molXB)
    if Dsat
        fᵪ = XCO₂_to_fO2(clim, P, T, sharpness, clim)
        return SA[  fₛ - fᵪ
                    fₗ - fᵪ
                    1 - sOₑₓ - mOₑₓ - Φ*clim/(1-clim) - fₑₓ], fₛ, fₗ, fᵪ
    else
        fᵪ = XCO₂_to_fO2(XCO₂, P, T, sharpness, clim)
        return SA[  fₛ - fₗ
                    fₗ - fᵪ
                    1 - sOₑₓ - mOₑₓ - Φ*XCO₂/(1-XCO₂)
                    fₑₓ], fₛ, fₗ, fᵪ
    end
end

function constrained_smOₑₓ_XCO₂_solver(y    :: SVector{4,Float64},                     # Initial transformed solution vector
                                       slim :: K, mlim :: K, clim :: K, fₑₓlim :: K, # Limits for smOₑₓ, mOₑₓ and XCO₂
                                       sfO2, sfO2⁻¹, ∂Sᵢ,                                      # Solid fO₂ and numerical derivative (MAGEMin)
                                       ϵ    :: K, damp :: K, niter :: Int64;           # Solver tolerance, dampening factor, and maximum iterations
                                       # Verbose parameters
                                       verb_flag::Int64,P::K,T::K,ϕ::K,Rs::K,Rf::K,
                                       TOex::K,TOₑₓ::K,p::K,TC::K,Mt::K,
                                       # Pre-computed parameters
                                       s::K,Φ::K,Φₘ::K,T₀::K,ΔCₚ::K,a::K,b::K,c::K,y1::K,y3::K,lowclip::K,
                                       y4::K,y5::K,y8::K,y9::K,IDV::K,_ln10::K,_T::K,Ys1::K,Ys2::K,sharpness::K,
                                       idxO::Int64, plotevo::Bool,SymXox::SVector{N, Symbol},verbose::Bool,debugging::Bool,
                                       dummy::Vector{Float64}, molXB::SVector{N, Float64}) where{K<:AbstractFloat, N}
    
    # Residual - fO₂ - Partitioning matrix for plotting
    mat = zeros(niter, 4, 3) # [Residuals, fO₂, Partitioning]

    # y margins
    y_lows = x_to_y(lowclip, lowclip, slim)
    y_highs = x_to_y(slim-lowclip, lowclip, slim)
    y_lowm = x_to_y(lowclip, lowclip, mlim)
    y_highm = x_to_y(mlim-lowclip, lowclip, mlim)
    y_lowc = x_to_y(lowclip, lowclip, clim)
    y_highc = x_to_y(clim-lowclip, lowclip, clim)
    y_lowfₑₓ = x_to_y(lowclip, lowclip, fₑₓlim)
    y_highfₑₓ = x_to_y(fₑₓlim-lowclip, lowclip, fₑₓlim)

    # Converged flag
    converged = false

    itout    = 0
    Dsat     = false
    switch   = false
    SM3      = @SMatrix zeros(3,3) # Static Matrix 3 × 3
    SM4      = @SMatrix zeros(4,4) # Static Matrix 4 × 4
    for it in 1:niter
        # Evaluate current stage
        y₁, y₂, y₃, y₄ = y
        # Force boundary just after regime transition
        Dsat && (y₃=y_highc)
        !Dsat && (y₄=y_lowfₑₓ)
        (switch&&Dsat)  && (y₄=x_to_y(0.001, lowclip, fₑₓlim))
        (switch&&!Dsat) && (y₃=0.999y_highc)
        switch = false
        # Transform back to original variables
        sOₑₓ, mOₑₓ, XCO₂, fₑₓ = y_to_x(y₁, lowclip, slim), y_to_x(y₂, lowclip, mlim), y_to_x(y₃, lowclip, clim), y_to_x(y₄, lowclip, fₑₓlim)
        # Iteration variables
        α   = evα(mOₑₓ, Φₘ)
        θₘ  = evθₘ(α, s)
        # Compute residual
        Fx, fₛ, fₗ, fᵪ = Rx(sfO2, P, T, sOₑₓ, mOₑₓ, XCO₂, fₑₓ, T₀, ΔCₚ, a, b, c, y1, y3, y4, y5, y8, y9, IDV, SymXox, dummy, idxO, _ln10, _T, s, Φ, Φₘ, molXB, sharpness, clim, Dsat)
        # Store values
        Dsat ? (mat[it,:,1] .= [Fx[1], Fx[2], Fx[3], 0.0]) : (mat[it,:,1] .= Fx)
        mat[it,1,2]  = sfO2(sOₑₓ)
        mat[it,2,2]  = Hirsch(T, mOₑₓ, T₀, ΔCₚ, a, b, c, y1, y3, y4, y5, y8, y9, IDV, SymXox, dummy, idxO, _ln10, _T, s, Φₘ, molXB)
        mat[it,3,2]  = XCO₂_to_fO2(XCO₂, P, T, sharpness, clim)
        mat[it,:,3] .= [sOₑₓ, mOₑₓ, XCO₂, fₑₓ]
        # Check convergence
        aR = maximum(abs.(Fx))
        if aR<=ϵ
            converged = true
            itout=it
            if verbose
                println("")
                if verb_flag==-1
                    println("---- Solution (P=$(P)GPa | T=$(T)K | ϕ=$(ϕ) | Rs=$(Rs) | Rf=$(Rf) | Mix=$p | TCarbon=$TC) ----")
                else
                    println("---- Solution (P=$(P)GPa | T=$(T)K | ϕ=$(ϕ) | TOₑₓ=$TOex | Mix=$p | TC=$TC) ----")
                end
                println("Shared fO₂ = $(sfO2(sOₑₓ)) |  residual = $aR")
                println("Total Oₑₓ budget = $(round((1e2TOₑₓ/Mt), digits=4))% of total mass")
                println("TOₑₓ partitioning → [$(round(1e2sOₑₓ, digits=4))% solid, $(round(1e2mOₑₓ, digits=4))% melt] stored as Fe₂O₃, [$(round(1e2fₑₓ, digits=4))% solid, $(round(1e2(1 - sOₑₓ - mOₑₓ - fₑₓ), digits=4))% melt] stored as CaCO₃ and CO₂ respectively")
                println("Melt XCO₂ = $(XCO₂)")
                println("Converged in $it iterations.")
            end
            return converged, mat, itout
        end
        # Partial derivatives
        ∂S  = ∂Sᵢ(sOₑₓ)
        ∂M  = ∂M∂mOₑₓ(Φₘ, Ys1, Ys2, α, θₘ, _ln10, a, molXB[3])
        ∂C  = ∂C∂XCO₂(XCO₂, _ln10, sharpness, clim)
        ∂3  = ∂3∂XCO₂(Φ, XCO₂)
        if debugging
            @printf "Iteration %d: sOₑₓ = %.4f (%.4f), mOₑₓ = %.4f (%.4f), XCO₂ = %.4f (%.4f), fₑₓ = %.4f (%.4f)\n" it sOₑₓ slim mOₑₓ mlim XCO₂ clim fₑₓ fₑₓlim
            @printf "\t(R₁=%.4f, R₂=%.4f, R₃=%.4f)" Fx[1] Fx[2] Fx[3]
            @printf "  (sfO₂=%.4f, mfO₂=%f, cfO₂=%.4f)\n" fₛ fₗ fᵪ
            @printf "\t(∂S=%.4f, ∂M=%f, ∂C=%.4f, ∂3=%.4f)\n\n" ∂S ∂M ∂C ∂3
        end
        # Jacobian inverse (Chain rule)
        ∂x∂y₁, ∂x∂y₂, ∂x∂y₃, ∂x∂y₄ = ∂x∂y(y₁, lowclip, slim), ∂x∂y(y₂, lowclip, mlim), ∂x∂y(y₃, lowclip, clim), ∂x∂y(y₄, lowclip, fₑₓlim)
        # @printf "  (∂x∂y₁=%.4f, ∂x∂y₂=%f, ∂x∂y₃=%.4f, ∂x∂y₄=%.4f)\n" ∂x∂y₁ ∂x∂y₂ ∂x∂y₃ ∂x∂y₄
        # println("")
        if Dsat
            J⁻¹  = SM3 + inv([   ∂S*∂x∂y₁       0.0         0.0  
                                   0.0        ∂M*∂x∂y₂      0.0   
                                 -∂x∂y₁        -∂x∂y₂      -∂x∂y₄ ]) # J = ∂Rᵢ∂yᵢ =  ∂Rᵢ∂xᵢ * ∂xᵢ∂yᵢ
        else
            J⁻¹  = SM4 + inv([   ∂S*∂x∂y₁    -∂M*∂x∂y₂      0.0        0.0
                                   0.0        ∂M*∂x∂y₂   -∂C*∂x∂y₃     0.0
                                 -∂x∂y₁       -∂x∂y₂     -∂3*∂x∂y₃     0.0
                                   0.0          0.0         0.0        ∂x∂y₄]) # J = ∂Rᵢ∂yᵢ =  ∂Rᵢ∂xᵢ * ∂xᵢ∂yᵢ
        end
        # Newton step
        dn = J⁻¹*Fx*damp; α = 1.0
        # Adaptive stepping
        (y[1]-dn[1]<y_lows)                && (α = min(α, 0.5*(y[1]-y_lows)/dn[1]))
        (y[1]-dn[1]>y_highs)               && (α = min(α, 0.5*(y[1]-y_highs)/dn[1]))
        (y[2]-dn[2]<y_lowm)                && (α = min(α, 0.5*(y[2]-y_lowm)/dn[2]))
        (y[2]-dn[2]>y_highm)               && (α = min(α, 0.5*(y[2]-y_highm)/dn[2]))
        if Dsat
            (y[4]-dn[3]<y_lowfₑₓ)          && (switch=true; Dsat=false)
            # (verbose && !Dsat)             && println("\t • De-Saturating. Releasing XCO₂ constraint.")
            (y[4]-dn[3]>y_highfₑₓ)         && (α = min(α, 0.5*(y[4]-y_highfₑₓ)/dn[3]))
            # Take step
            y = y - α*[dn[1], dn[2], 0.0, dn[3]]
        else
            (y[3]-dn[3]<y_lowc)            && (α = min(α, 0.5*(y[3]-y_lowc)/dn[3]))
            (!Dsat && y[3]-dn[3]>y_highc)  && (switch=true; Dsat=true)
            # (verbose && Dsat)              && println("\t • Saturating. Locking XCO₂ to $clim .")
            # Take step
            y = y - α*dn
        end
        # Output if not converged
        if verbose && (it==niter)
            if verb_flag==-1
                println("---- Solution (P=$(P)GPa | T=$(T)K | ϕ=$(ϕ) | Rs=$(Rs) | Rf=$(Rf) | Mix=$p | TCarbon=$TC) ----")
            else
                println("---- Solution (P=$(P)GPa | T=$(T)K | ϕ=$(ϕ) | TOₑₓ=$TOex | Mix=$p | TC=$TC) ----")
            end
            println("Shared fO₂ = $(sfO2(sOₑₓ)) |  residual = $aR")
            println("Total Oₑₓ budget = $(round((1e2TOₑₓ/Mt), digits=4))% of total mass")
            println("TOₑₓ partitioning → [$(round(1e2sOₑₓ, digits=4))% solid, $(round(1e2mOₑₓ, digits=4))% melt] stored as Fe₂O₃, [$(round(1e2fₑₓ, digits=4))% solid, $(round(1e2(1 - sOₑₓ - mOₑₓ - fₑₓ), digits=4))% melt] stored as CaCO₃ and CO₂ respectively")
            println("Melt XCO₂ = $(XCO₂)")
            println("Did not converged in $it iterations.")
        end
    end
    return false, mat, niter
    
end

# Variable transformations
@inline y_to_x(y, minV, maxV) = minV + (maxV - minV) / (1.0 + exp(-y))
@inline x_to_y(x, minV, maxV) = -log((maxV - x)/(x + minV))
@inline ∂x∂y(y, minV, maxV)   = (maxV - minV)*exp(-y)/(1.0+exp(-y))^2
# Variables carbon
@inline evΦ(TOₑₓ, molMf)   = 2*mm.O*molMf/TOₑₓ     # Conversion factor for XCO₂
# Variables melt
@inline evΦₘ(TOₑₓ, Mf)   = TOₑₓ/Mf/mm.O          # Conversion factor for melt
@inline evα(mOₑₓ, Φₘ)    = Φₘ*mOₑₓ               # Mass fraction of TOₑₓ in melt → non-normalized mass of Oₑₓ in melt
@inline evθₘ(α, s)       = s - α                 # Molar normalization factor
# Partial derivatives
@inline ∂S∂sOₑₓ(sfO2, sOlist)                        = @views diff(sfO2)./diff(sOlist)
@inline ∂C∂XCO₂(XCO₂, _ln10, sharpness, maxV)        = _ln10*(1/XCO₂ + sharpness/(maxV-XCO₂))
@inline ∂M∂mOₑₓ(Φₘ, Ys1, Ys2, α, θₘ, _ln10, a, XFeO) = -Φₘ/a * ( θₘ^(-2)*(Ys1 + 2Ys2/θₘ) - _ln10*(1/α + 2/(XFeO - 2α)))
@inline ∂3∂XCO₂(Φ, XCO₂)                             = Φ*(1/(1-XCO₂)^2)