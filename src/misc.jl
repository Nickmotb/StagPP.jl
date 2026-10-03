# ============================================================
# ================ Excess oxygen partitioning ================
# ============================================================

        # Model:
        # Assess partitioning of TOₑₓ into solid, molten and fluid tracer + oxidation reduction of parcel carbon through equilibrium XCO₂.
        #
        # Inputs:
        # - Pressure (P)
        # - Temperature (T)
        # - Total oxygen budget either as wt% of solid (TOₑₓ) or solid and melt Fe³⁺/Feᵀ ratios (Rs, Rf)
        # - Equilibrium melt fraction (ϕ)
        # - Available reduced carbon in the parcel as wt% of solid tracer (TC)
        #
    # ================================
    # ========== Core idea  ==========
    # ================================
        #
        #   To fully couple the evolution fo the redox state in geodynamical simulations with thermodynamics we quantify the partitioning of a
        #   virtual excess oxygen budget representative of the oxidation potential of a parcel, into solid, molten, and fluid material.
        #   We couple the mass partitioning through the system's fO₂ which acts as a master variable, according to the thermodynamic equilibrium
        #   constraint → μO_solid = μO_melt = μO_COH
        #
        #   Quantifying solid and melt Oₑₓ contributions is straight forward once the governing fO₂ is known, as it fO₂ can be mapped to XFe₂O₃ in
        #   both phases through either Gibbs Free Energy minimization or separate oxidation models for melt. The remaining sink of Oₑₓ required to
        #   properly solve the system at very high fO₂ lies in the COH fluid which has to be treated separately, through an second inner iterative
        #   solve.
        #
    # ==============================================
    # ========== C-O-H fluid description  ==========
    # ==============================================
        #
        #   A mapping must be found between the oxidation state of the C-O-H fluid and the actual contribution to the Oₑₓ mass sink in order
        #   to balance the conservation constraint. Because of this, we define C-O-H fluids as a mixture of 6 different components (more can be easily added):
        #       
        #       • H₂ - C - H₂O - CH₄ - CO₂ - CO
        #
        #   We know the oxidation reactions that govern these compounds
        #
        #       • H₂  +  1/2 O₂  → H₂O  ==>  K₁ = fH₂O / (√fO₂ × fH₂)
        #       • C   +  2 H₂    → CH₄  ==>  K₂ = fH₄  / (aC × fH₂²)        ;   aC being 1 when below the buffer (g/d present)
        #       • C   +  O₂      → CO₂  ==>  K₃ = fCO₂ / (aC × fO₂)
        #       • C   +  1/2 O₂  → CO   ==>  K₄ = fCO  / (aC × √fO₂)
        #
        #   Knowing therefore the molar proportions of all these components (Xᵢ) yields a direct link to the mass of Oₑₓ retained.
        #   In order to find Xᵢ, we must define the problem in a closed form. We can do this by using component partial pressures (Pᵢ).
        #   Partial pressures can be rewritten in terms of component fugacities (fᵢ) and ideality coefficients (γᵢ):
        #
        #       • Pₜₒₜ = ∑ Pᵢ = ∑ fᵢ / γᵢ
        #
        #   One can then find a way to obtain values for the component fugacities, for example by redefining them in terms of known variables like fO₂ and aC.
        #   If the equilibrium constants [Kᵢ(P,T)] are known, this leaves only one unknown fH₂ that must be solved for in order to satisfy the pressure constraint.
        #   Retrieving fugacities and indeality coefficients fully speciates the COH fluid, leading to fOₑₓ being constrained, that can be passed one level above
        #   to assess the convergence of the global mass balance condition (equation 3.).
        #
    # ===========================================
    # ========= Mathematical framework ==========
    # ===========================================
        #                                                                       
        # Independent variables:                                    
        # - Equilibrium mass fraction of TOₑₓ in the solid (sOₑₓ)   
        # - Equilibrium mass fraction of TOₑₓ in the melt  (mOₑₓ)   
        # - Equilibrium molar XCO₂ in the silicate melt (XCO₂)              
        # - Di-hydrogen fugacity (fH₂) → leads to fOₑₓ (Equilibrium mass fraction of TOₑₓ in the fluid)
        #                                                                       
        #
        # Mathematical framework:
        #
        #   • Enstatite - Magnesite - Olivine - Graphite buffer (EMOG) → Stagno et al. (2010) ↔ fO₂(P,T)  => graphite/diamond present only when below EMOG.
        #
        #   • Buffered regime (fO₂ < fO₂_EMOG):
        #
        #       (1) sfO₂(P,T,sOₑₓ) - mfO₂(P,T,mOₑₓ) = 0                                                                      <--- Chemical potential constraint
        #       (2) mfO₂(P,T,mOₑₓ) - fO₂(P,T,XCO₂_solubility_law) = 0                                                        <--- Buffered XCO₂ solubility law 
        #       (3) 1 - sOₑₓ - mOₑₓ - fOₑₓ(fH₂,sfO₂) - Φ*XCO₂/(1-XCO₂) = 0                                                   <--- Global mass balance
        #       (4) Pₜₒₜ = K₁(√sfO₂ × fH₂)/γ₁ + K₂(aC × fH₂²)/γ₂ + K₃(aC × sfO₂)/γ₃ + K₄(aC × √sfO₂)/γ₄ + fO₂/γ₅ + fH₂/γ₆    <--- C-O-H contribution
        #
        #   • Jacobian : [∂(1)∂sOₑₓ ∂(1)∂mOₑₓ ∂(1)∂XCO₂ ∂(1)∂fH₂        [             ∂sfO₂∂sOₑₓ                -∂Hirsch∂mOₑₓ        0             0
        #                 ∂(2)∂sOₑₓ ∂(2)∂mOₑₓ ∂(2)∂XCO₂ ∂(2)∂fH₂                          0                      ∂Hirsch∂mOₑₓ   -∂EMODG∂XCO₂       0
        #                 ∂(3)∂sOₑₓ ∂(3)∂mOₑₓ ∂(3)∂XCO₂ ∂(3)∂fH₂    =    -∂sfO₂∂sOₑₓ - ∂fOₑₓ/∂sfO₂*∂sfO₂∂sOₑₓ        -1           -∂3∂XCO₂      -∂fOₑₓ∂fH₂
        #                 ∂(4)∂sOₑₓ ∂(4)∂mOₑₓ ∂(4)∂XCO₂ ∂(4)∂fH₂]               ∂COH_P/∂sfO₂ * ∂sfO₂∂sOₑₓ             0              0          ∂COH_P∂fH₂
        #
        #
        #   • Un-buffered regime (fO₂ > fO₂_EMOG):
        #
        #       (1) sfO₂(P,T,sOₑₓ) - mfO₂(P,T,mOₑₓ) = 0
        #       (2) 1 - sOₑₓ - mOₑₓ - fOₑₓ = 0
        #       (3) Pₜₒₜ = K₁(√fO₂ × fH₂) + K₂(aC × fH₂²) + K₃(aC × fO₂) + K₄(aC × √fO₂)
        #
        #
        # Variable extentions:
        # Φ          = 2*mm.O*molMf/TOₑₓ
        # sw, s      = ∑(oxᵢ*mmᵢ), ∑(oxᵢ/mmᵢ)
        #
        # Tools:
        # - Hirschmann 2022 melt mapping from XFe₂O₃ (Oₑₓ) ↔ fO₂
        # - Stixrude and Bertelloni 2024 (MAGEMin) solid mapping from XFe₂O₃ (Oₑₓ) ↔ fO₂
        # - Stagno and Frost 2010 parameterization for silicate melt XCO₂ solubility as a function of P-T-fO₂
    
    # ===

    # function partition_Oₑₓ(P::K, T::K, p::K, ϕ::K, TOex::K, TC::K, H::K; Rs::K=-1.0, Rf::K=-1.0, nr=25, niter=1000, 
    #                         verbose=false, data=nothing, Rspace=false, plotevo=false, damp=0.25, debugging=false, saveplotevery=false, savein="") where {K <: Real}

    #     # Hirschmann
    #     a=0.1917; y1=-520.46; y2=-185.37; y3=494.39; y4=1838.34; y5=2888.48; y8=-1245.09; y9=-1156.86

    #     # Endmember bulks
    #     XH      = @SVector [0.4347, 0.4597, 0.0835, 0.0090, 0.0100, 0.0001, 0.0030, 0.0] # mass fraction
    #     XB      = @SVector [0.5097, 0.0988, 0.0718, 0.1268, 0.1698, 0.0225, 0.0007, 0.0] # mass fraction
    #     Xox     = @SVector ["SiO2", "MgO", "FeO", "CaO", "Al2O3", "Na2O", "Cr2O3", "O"]
    #     SymXox  = @SVector [:SiO2, :MgO, :FeO, :CaO, :Al2O3, :Na2O, :Cr2O3, :O]
    #     mmXox   = @SVector [mm.SiO2, mm.MgO, mm.FeO, mm.CaO, mm.Al2O3, mm.Na2O, mm.Cr2O3, mm.O]
    #     X       = p*XB + (1-p)*XH
    #     molXB   = XB./mmXox; molXB = molXB ./ sum(molXB)
    #     molX    = X./mmXox;  molX  = molX ./sum(molX)

    #     # Cell mass [kg] → solid | molten tracer + COH mass
    #     Mt = 1e3 # 1 kg
    #     Mf=ϕ*Mt; Ms=Mt*(1-ϕ); Mc=TC*Mt; Mh=H*Mt
    #     molMf = sum(Mf.*molXB)   # Mf in mols (unoxidized)

    #     # Carbon constraints
    #     molCav      = Mc/mm.C                   # Mols of available carbon
    #     maxXCO₂_raw = molCav/(molMf+molCav)     # Maximal XCO₂ allowed by carbon

    #     # Water constraints
    #     molH2O      = Mh/mm.H2O                 # Mols of available water

    #     # IDV = ∫ΔVdP for melts
    #     IDV = solve_∫ΔVdP([P-0.05P, P, P+0.05P],[T-0.05T, T, T+0.05T])[2,2,1]

    #     # Flag whether to manage MAGEMin initialisation and finalization
    #     flag = isnothing(data)

    #     # Index of oxygen component
    #     idxO = findfirst(Xox.=="O")
        
    #     # === Allocate iterative memory and create bulk structures
    #         sol        = @SVector zeros(4)                                                  # Solution Vector (Raw)
    #         dummy      = zeros(length(XB))                                                  # for Hirschmann calls
    #         Xdummy     = (X=Vector{Float64}(X), Xox=Vector{String}(Xox), mm=get_Xoxmm(Xox)) # For oxidizing calls
    #         if Rs>=0.0 || Rf>=0.0
    #             (Rs<0.0) && (Rs=0.0); (Rf<0.0) && (Rf=0.0)                                                      # If an Fe³⁺ ratio is passed, set the inputs up
    #             Xs = oxidize_bulk(X, Rs, Xdummy, FeFormat="FeO_O", wt_out=true, frac=true, SymXox=SymXox);      # Oxidized solid bulk comp
    #             Xmo = oxidize_bulk(XB, Rf, Xdummy, FeFormat="FeO_O", wt_out=true, frac=true, SymXox=SymXox);    # Oxidized melt bulk comp
    #         end
    #         Xm = oxidize_bulk(XB, 0.0, Xdummy, FeFormat="FeO_O", wt_out=true, frac=true);                       # Unoxidized melt bulk comp

    #     # Sum TOₑₓ contributions from solid and molten tracer Oₑₓ (TOex -> input total mass O budget, TOₑₓ -> system total mass Oₑₓ budget)
    #     if Rs>=0.0 || Rf>=0.0
    #         verb_flag = -1; TOₑₓ = Xs.O*Ms + Xmo.O*Mf # If ratios passed, this considers both compositions
    #     else
    #         verb_flag = 1; TOₑₓ = TOex*Mt # If TOex passed, a fraction of the total mass becomes the Oₑₓ budget
    #     end
    #     if TOₑₓ==0.0 && verbose
    #         println("---- Solution (P=$(P)GPa | T=$(T)K | ϕ=$(ϕ) | Rs=$(Rs) | Rf=$(Rf) | H₂O=($H) | Source mix = $p) ----")
    #         println("Total Oₑₓ budget = 0.0 kg")
    #         println("TOₑₓ partitioning → [0.0% solid, 0.0% melt]")
    #         return
    #     end

    #     # Compute partial derivative variables
    #     s           = sum(molXB)                                                    # ∑oxᵢ in mols
    #     Φ, Φₘ       = evΦ(TOₑₓ, molMf), evΦₘ(TOₑₓ, Mf)                              # Conversion factors       
    #     _ln10, _T   = 1/log(10), 1/T                                                # Pre computed constants
    #     Ys1         = (y1*molXB[1] + y3*molXB[2] + y4*molXB[4] + y5*molXB[6])*_T    # Sum of linear parameterized molar components
    #     Ys2         = molXB[1]*(y8*molXB[5] + y9*molXB[2])*_T                       # Sum of non-linear parameterized molar components

    #     # -- Solver boundary margin and XCO₂ sharpness parameter
    #     lowclip   = 1e-18
    #     sharpness = 0.0

    #     # === Generate solid fO₂ space
    #         Rlist = LinRange(0.00001, 0.20, nr)
    #         Xlist = Vector{Vector{Float64}}(undef, nr); 
    #         sOₑₓlist = zeros(nr)
    #         for i in 1:nr
    #             Xl = oxidize_bulk(X, Rlist[i], Xdummy, wt_out=true, frac=true, FeFormat="FeO_O", SymXox=SymXox)
    #             Xlist[i] = [getfield(Xl, f) for f in SymXox]
    #             sOₑₓlist[i] = Xl.O
    #         end
    #     # -- Define maximal sOₑₓ
    #         maxsOₑₓ_uncapped  = (0.5molX[3]*mm.O)/(sum(X) - 0.5molX[3]*mm.O)*(Ms/TOₑₓ)
    #         maxsOₑₓ = min(maxsOₑₓ_uncapped, 1.0)
    #     # -- Minimizer call
    #         flag && (data = Initialize_MAGEMin("sb24", verbose=false));
    #         out = multi_point_minimization(10P*ones(nr), k2c(T)*ones(nr), data, X=Xlist, Xoxides=Vector{String}(Xox), name_solvus=true, sys_in="wt", progressbar=false)
    #         flag && Finalize_MAGEMin(data);
    #     # -- Create interpolation object
    #         sfO2   = extrapolate(Interpolations.interpolate((sOₑₓlist.*Ms./TOₑₓ,), [out[i].fO2 for i in eachindex(out)], Gridded(Linear())), Flat())
    #         sample_sOlist = LinRange(0.8lowclip, 1.3maxsOₑₓ, 100)
    #         sample_sOlist05 = 0.5(sample_sOlist[1:end-1] + sample_sOlist[2:end])
    #         sampled_sfO2 = sfO2(sample_sOlist)
    #     # -- Solid partial derivative
    #         ∂Sᵢ   = extrapolate(Interpolations.interpolate((sample_sOlist05,), ∂S∂sOₑₓ(sampled_sfO2, sample_sOlist), Gridded(Linear())), Line())
            
    #     # -- Compute independent boundaries
    #     maxmOₑₓ_uncapped  = (0.5molXB[3]*mm.O)/(sum(XB) + 0.5molXB[3]*mm.O)*(Mf/TOₑₓ)
    #     maxmOₑₓ, maxXCO₂  = min(maxmOₑₓ_uncapped, 1.0), max(min(1.0, maxXCO₂_raw), 0.0)
    #     minfH₂, maxfH₂ = -9, 2
    #     # -- Initialise (static)solution vector → y = [sOₑₓ, mOₑₓ, XCO₂, fH₂]
    #     y = sol + [x_to_y(5e-4maxsOₑₓ, lowclip, maxsOₑₓ), x_to_y(0.5maxmOₑₓ, lowclip, maxmOₑₓ), x_to_y(0.01maxXCO₂, lowclip, maxXCO₂), x_to_y(0.0, minfH₂, maxfH₂)]
    #     # -- Define convergence tolerance (ϵ)
    #     ϵ = 2e-15          
    #     # -- Wrap parameters and call solver
    #     params = (; verb_flag, P, T, ϕ, Rs, Rf, TOex, TOₑₓ, p, TC, H, Φ, s, Φₘ,
    #                     IDV, Mh, SymXox, dummy, idxO, _ln10, _T, molXB, a,
    #                         Ys1, Ys2, plotevo, verbose, lowclip, Mt, sharpness, debugging)
    #     converged, mat, itout = constrained_smOₑₓ_XCO₂_solver(y, maxsOₑₓ, maxmOₑₓ, maxXCO₂, minfH₂, maxfH₂, sfO2, ∂Sᵢ, ϵ, damp, niter; params...)
    #     # -- Plot evolution if requested
    #     x = mat[itout, :, 3]
    #     if plotevo
    #         solidclr = :black
    #         meltclr  = :red
    #         co2clr   = :green
    #         fₑₓclr   = :purple
    #         itstart  = converged ? Int(floor(0.1itout)) : 1
    #         if saveplotevery
    #             for it in itstart+1:itout
    #                 x .= mat[it, :, 3]
    #                 # Plot evolution
    #                 fig = Figure(size=(1800, 800))
    #                 ax = Axis(fig[1,1], ylabel=L"Solution\;residual", xlabel=L"Iterations", xgridvisible=false, ygridvisible=false, ylabelsize=25, xlabelsize=25, yscale=log10); 
    #                 scatterlines!(ax, itstart:it, abs.(mat[itstart:it,1,1]), color=meltclr, label="(eq. 1) Solid ↔ melt fO₂ equilibrium (R₁ = $(Float16(mat[it,1,1])))",marker=:rect,strokewidth=1.1)
    #                 scatterlines!(ax, itstart:it, abs.(mat[itstart:it,2,1]), color=co2clr,label="(eq. 2) Solid ↔ EMODG fO₂ equilibrium (R₂ = $(Float16(mat[it,2,1])))",marker=:rect,strokewidth=1.1)
    #                 scatterlines!(ax, itstart:it, abs.(mat[itstart:it,3,1]), color=:orange,label="(eq. 3) Mass conservation (R₃ = $(Float16(mat[it,3,1])))",marker=:rect,strokewidth=1.1)
    #                 lines!(ax, [1, it], [ϵ, ϵ], linestyle=:dash, color=:gray, alpha=0.5)
    #                 text!(ax, 1, 1.9ϵ, text="Convergence ϵ = $ϵ", font=:italic, alpha=0.5)
    #                 axislegend(ax, position=:rt)

    #                 ax = Axis(fig[2,1], ylabel=L"log\;fO_2", xlabel=L"Iterations", xgridvisible=false, ygridvisible=false, ylabelsize=25, xlabelsize=25); 
    #                 scatterlines!(ax, itstart:it, mat[itstart:it,1,2], color=solidclr, label="Solid (fO₂ = $(round(mat[it,1,2], digits=4)))",marker=:rect,strokewidth=1.1)
    #                 scatterlines!(ax, itstart:it, mat[itstart:it,2,2], color=meltclr,label="Melt (fO₂ = $(round(mat[it,2,2], digits=4)))",marker=:rect,strokewidth=1.1)
    #                 scatterlines!(ax, itstart:it, mat[itstart:it,3,2], color=co2clr,label="EMODG (fO₂ = $(round(mat[it,3,2], digits=4)))",marker=:rect,strokewidth=1.1)
    #                 axislegend(ax, position=:rt)
    #                 converged && ylims!(ax, mat[itout,3,2]-0.5, mat[itout,3,2]+1)

    #                 ax = Axis(fig[1:2,2], xlabel=L"Iterations", ylabel=L"Fraction\;of\;TO_{ex}", rightspinecolor=:green, xgridvisible=false, ygridvisible=false, ylabelsize=25, xlabelsize=25); 
    #                 lines!(ax, itstart:it, mat[itstart:it,1,3],label="In solid Fe ($(round(x[1], digits=3)) TOₑₓ)",color=solidclr,linewidth=2.0)
    #                 lines!(ax, itstart:it, mat[itstart:it,2,3], color=meltclr,label="In melt Fe ($(round(x[2], digits=3)) TOₑₓ)",linewidth=2.0)
    #                 lines!(ax, itstart:it, mat[itstart:it,4,3], color=fₑₓclr,label="Precipitated ($(round(x[4], digits=3)) TOₑₓ)",linewidth=2.0)
    #                 scatter!(ax, 1, mat[1,3,3], label="Melt XCO₂ = $(round(x[3], digits=5))", alpha=0.0)
    #                 axislegend(ax, position=:rt, framevisible=true)
    #                 # Mark ceilings
    #                 scatterlines!(ax, [itstart, it], [maxsOₑₓ, maxsOₑₓ], alpha=0.3, color=solidclr,marker=:rect,strokewidth=1.1); text!(ax, 0.1it, 1.02maxsOₑₓ, text="Solid Fe³⁺ cap = $(round(maxsOₑₓ, digits=3))"*(maxsOₑₓ_uncapped>1.0 ? " ($(round(maxsOₑₓ_uncapped, digits=3)))" : ""), fontsize=12, font=:italic, color=solidclr)
    #                 scatterlines!(ax, [itstart, it], [maxmOₑₓ, maxmOₑₓ], alpha=0.3, color=meltclr,marker=:rect,strokewidth=1.1); text!(ax, 0.4it, 1.02maxmOₑₓ, text="Melt Fe³⁺ cap = $(round(maxmOₑₓ, digits=3))"*(maxmOₑₓ_uncapped>1.0 ? " ($(round(maxmOₑₓ_uncapped, digits=3)))" : ""), fontsize=12, font=:italic, color=meltclr)
    #                 ax2 = Axis(fig[1:2,2], ylabel=L"XCO_2", yaxisposition=:right, ylabelcolor=co2clr, ytickcolor=co2clr, yticklabelcolor=co2clr, xgridvisible=false, ygridvisible=false); hidespines!(ax2, :l, :t, :b, :r); hidexdecorations!(ax2)
    #                     lines!(ax2, itstart:it, mat[itstart:it,3,3], color=co2clr,label="In melt CO₂ ($(round(1 - x[1] - x[2], digits=3)) TOₑₓ)",linewidth=2.0)
    #                     scatterlines!(ax2, [itstart, it], [maxXCO₂, maxXCO₂], alpha=0.3, color=co2clr,marker=:rect,strokewidth=1.1); text!(ax2, 0.7it, 1.02maxXCO₂, text="Carbon limited XCO₂ = $(round(maxXCO₂, digits=3))", fontsize=12, font=:italic, color=co2clr)
    #                 # Limits
    #                     ul = max(maxmOₑₓ, maxXCO₂)
    #                     ylims!(ax, 0.0, 1.3ul); 
    #                     ylims!(ax2, 0.0, 1.3ul);
    #                 GLMakie.save("./+img/+toex_seq/it_$it.png", fig)
    #             end
    #         else
    #             # Plot evolution
    #             fig = Figure(size=(1800, 800))
    #             ax = Axis(fig[1,1], ylabel=L"Solution\;residual", xlabel=L"Iterations", xgridvisible=false, ygridvisible=false, ylabelsize=25, xlabelsize=25, yscale=log10); 
    #             scatterlines!(ax, itstart:itout, abs.(mat[itstart:itout,1,1]), color=meltclr, label="(eq. 1) Solid ↔ melt fO₂ equilibrium (R₁ = $(Float16(mat[itout,1,1])))",marker=:rect,strokewidth=1.1)
    #             scatterlines!(ax, itstart:itout, abs.(mat[itstart:itout,2,1]), color=co2clr,label="(eq. 2) Solid ↔ EDDOG fO₂ equilibrium (R₂ = $(Float16(mat[itout,2,1])))",marker=:rect,strokewidth=1.1)
    #             scatterlines!(ax, itstart:itout, abs.(mat[itstart:itout,3,1]), color=:orange,label="(eq. 3) Mass conservation (R₃ = $(Float16(mat[itout,3,1])))",marker=:rect,strokewidth=1.1)
    #             lines!(ax, [1, itout], [ϵ, ϵ], linestyle=:dash, color=:gray, alpha=0.5)
    #             text!(ax, 1, 1.9ϵ, text="Convergence ϵ = $ϵ", font=:italic, alpha=0.5)
    #             axislegend(ax, position=:rt)
    #             # converged && ylims!(ax, -5ϵ, 7ϵ)

    #             ax = Axis(fig[2,1], ylabel=L"log\;fO_2", xlabel=L"Iterations", xgridvisible=false, ygridvisible=false, ylabelsize=25, xlabelsize=25); 
    #             scatterlines!(ax, itstart:itout, mat[itstart:itout,1,2], color=solidclr, label="Solid (fO₂ = $(round(mat[itout,1,2], digits=4)))",marker=:rect,strokewidth=1.1)
    #             scatterlines!(ax, itstart:itout, mat[itstart:itout,2,2], color=meltclr,label="Melt (fO₂ = $(round(mat[itout,2,2], digits=4)))",marker=:rect,strokewidth=1.1)
    #             scatterlines!(ax, itstart:itout, mat[itstart:itout,3,2], color=co2clr,label="EDDOG (fO₂ = $(round(mat[itout,3,2], digits=4)))",marker=:rect,strokewidth=1.1)
    #             axislegend(ax, position=:rt)
    #             converged && ylims!(ax, mat[itout,3,2]-0.5, mat[itout,3,2]+1)

    #             ax = Axis(fig[1:2,2], xlabel=L"Iterations", ylabel=L"Fraction\;of\;TO_{ex}", rightspinecolor=:green, xgridvisible=false, ygridvisible=false, ylabelsize=25, xlabelsize=25); 
    #             lines!(ax, itstart:itout, mat[itstart:itout,1,3],label="In solid Fe ($(round(x[1], digits=5)) TOₑₓ)",color=solidclr,linewidth=2.0)
    #             lines!(ax, itstart:itout, mat[itstart:itout,2,3], color=meltclr,label="In melt Fe ($(round(x[2], digits=5)) TOₑₓ)",linewidth=2.0)
    #             lines!(ax, itstart:itout, mat[itstart:itout,4,3], color=fₑₓclr,label="Precipitated ($(round(x[4], digits=5)) TOₑₓ)",linewidth=2.0)
    #             scatter!(ax, 1, mat[1,3,3], label="Melt XCO₂ = $(round(x[3], digits=5))", alpha=0.0)
    #             axislegend(ax, position=:rt, framevisible=true)
    #             # Mark ceilings
    #             scatterlines!(ax, [itstart, itout], [maxsOₑₓ, maxsOₑₓ], alpha=0.3, color=solidclr,marker=:rect,strokewidth=1.1); text!(ax, 0.1itout, 1.02maxsOₑₓ, text="Solid Fe³⁺ cap = $(round(maxsOₑₓ, digits=3))"*(maxsOₑₓ_uncapped>1.0 ? " ($(round(maxsOₑₓ_uncapped, digits=3)))" : ""), fontsize=12, font=:italic, color=solidclr)
    #             scatterlines!(ax, [itstart, itout], [maxmOₑₓ, maxmOₑₓ], alpha=0.3, color=meltclr,marker=:rect,strokewidth=1.1); text!(ax, 0.4itout, 1.02maxmOₑₓ, text="Melt Fe³⁺ cap = $(round(maxmOₑₓ, digits=3))"*(maxmOₑₓ_uncapped>1.0 ? " ($(round(maxmOₑₓ_uncapped, digits=3)))" : ""), fontsize=12, font=:italic, color=meltclr)
    #             ax2 = Axis(fig[1:2,2], ylabel=L"XCO_2", yaxisposition=:right, ylabelcolor=co2clr, ytickcolor=co2clr, yticklabelcolor=co2clr, xgridvisible=false, ygridvisible=false); hidespines!(ax2, :l, :t, :b, :r); hidexdecorations!(ax2)
    #                 lines!(ax2, itstart:itout, mat[itstart:itout,3,3], color=co2clr,label="In melt CO₂ ($(round(1 - x[1] - x[2], digits=3)) TOₑₓ)",linewidth=2.0)
    #                 scatterlines!(ax2, [itstart, itout], [maxXCO₂, maxXCO₂], alpha=0.3, color=co2clr,marker=:rect,strokewidth=1.1); text!(ax2, 0.7itout, 1.02maxXCO₂, text="Carbon limited XCO₂ = $(round(maxXCO₂, digits=3))", fontsize=12, font=:italic, color=co2clr)
    #                 # Limits
    #                 ul = max(maxmOₑₓ, maxXCO₂) # maxsOₑₓ, 
    #                 d = itout-itstart
    #                 xlims!(ax, itstart-0.05d, itout+0.05d)
    #                 xlims!(ax2, itstart-0.05d, itout+0.05d)
    #                 ylims!(ax, 0.0, 1.3ul); 
    #                 ylims!(ax2, 0.0, 1.3ul);
    #             !isempty(savein) && GLMakie.save("./" * savein * ".png", fig)
    #         end
    #         display(fig)
    #     end

    #     if Rspace
    #         # Return partitioning
    #         return converged, x
    #     else
    #         return converged, mat, itout
    #     end

    # end

    # function Hirsch(T, mOₑₓ, IDV, SymXox, dummy, idxO, _ln10, _T, s, Φₘ, molXB)
    #     # Hirschmann 2022 parameters
    #     a=0.1917; b=-1.961; c=4158.1; ΔCₚ=33.25; T₀=1673.15; y1=-520.46; y2=-185.37; y3=494.39; y4=1838.34; y5=2888.48; y8=-1245.09; y9=-1156.86
    #     # Checks
    #     @assert (:O∈SymXox && :FeO∈SymXox) "This function requires FeO + O format!"
    #     # Reset dummy
    #     dummy .= molXB
    #     # Compute non-normalized mols of Oₑₓ in melt
    #     mfOₑₓ = Φₘ*mOₑₓ
    #     # Oxidize
    #     dummy[idxO] = mfOₑₓ
    #     # Normalize
    #     dummy./=(s+mfOₑₓ)
    #     # Construct bulk and assess hardlimit
    #     Xl = Cbulk((; zip(SymXox, dummy)...))
    #     (Xl.O>=0.5Xl.FeO) && (return 14) # Too much oxygen!! Above hard-limit.
    #     # Compute mfO2
    #     mfO2 = (log10(Xl.O/(Xl.FeO-2Xl.O)) - b - c*_T + (ΔCₚ/R*_ln10 * (1 - T₀*_T - log(T/T₀))) + IDV/(1e-3R)*_T*_ln10 
    #                         - _T*(y1*Xl.SiO2 + y3*Xl.MgO + y4*Xl.CaO + y5*Xl.Na2O + y8*Xl.SiO2*Xl.Al2O3 + y9*Xl.SiO2*Xl.MgO))/a
    #     return mfO2
    # end

    # function Hirsch⁻¹(T, eqfO2, IDV, SymXox, dummy, idxO, _ln10, _T, s, Φₘ, molXB, Ys1, Ys2, mlim)
    #     iguess = 0.05mlim; ϵ = 1e-2
    #     for it in 1:20
    #         R      = eqfO2 - Hirsch(T, iguess, IDV, SymXox, dummy, idxO, _ln10, _T, s, Φₘ, molXB)
    #         (abs(R)<=ϵ) && break
    #         α      = evα(iguess, Φₘ)
    #         θₘ     = evθₘ(α, s)
    #         ∂M     = -∂M∂mOₑₓ(Φₘ, Ys1, Ys2, α, θₘ, _ln10, a, molXB[3])
    #         iguess = iguess - R/∂M
    #     end
    #     return iguess
    # end

    # # Modified from Stagno and Frost XCO2 equilibrium (2 Oex per 1 CO₂) | cOₑₓ in mass fraction of TOₑₓ
    # function XCO₂_to_fO2(XCO₂, P, T, sharpness, clim)
    #     # Checks
    #     @assert XCO₂>=0.0 "XCO₂ cannot be zero."
    #     # Limit to rexplored ranges
    #     # P = Pin > 11. ? 11. : Pin < 2.5 ? 2.5 : Pin
    #     # T = Tin >= c2k(1600.) ? c2k(1600.) : Tin < c2k(1100.) ? c2k(1100.) : Tin
    #     # Compute logfO₂
    #     if sharpness==0.0
    #         return 5.44 - 21380/T + 0.078(1e4P-1)/T + log10(XCO₂)
    #     else
    #         return 5.44 - 21380/T + 0.078(1e4P-1)/T + log10(XCO₂) - sharpness*log10(clim-XCO₂)
    #     end
    # end

    # function Rx(y₁, y₂, y₃, y₄, dummy, params)
    #     (;sfO2, P, T, IDV, Mh, SymXox, idxO, _ln10, _T, s, Φ, Φₘ, molXB, sharpness, Dsat, lowclip, slim, mlim, clim, minfH₂, maxfH₂, TOₑₓ) = params
    #     sOₑₓ, mOₑₓ, XCO₂, fH₂ = y_to_x(y₁, lowclip, slim), y_to_x(y₂, lowclip, mlim), y_to_x(y₃, lowclip, clim), y_to_x(y₄, minfH₂, maxfH₂)
    #     fₛ = sfO2(sOₑₓ)
    #     fₗ = Hirsch(T, mOₑₓ, IDV, SymXox, dummy, idxO, _ln10, _T, s, Φₘ, molXB)
    #     fᵪ = XCO₂_to_fO2(XCO₂, P, T, sharpness, clim)
    #     # if Dsat
    #     #     fᵪ = XCO₂_to_fO2(clim, P, T, sharpness, clim)
    #     #     return SA[  fₛ - fᵪ
    #     #                 fₗ - fᵪ
    #     #                 1 - sOₑₓ - mOₑₓ - Φ*clim/(1-clim) - fₑₓ], fₛ, fₗ, fᵪ, sOₑₓ, mOₑₓ, XCO₂, fₑₓ
    #     # else
    #     #     fᵪ = XCO₂_to_fO2(XCO₂, P, T, sharpness, clim)
    #     #     return SA[  fₛ - fₗ
    #     #                 fₗ - fᵪ
    #     #                 1 - sOₑₓ - mOₑₓ - Φ*XCO₂/(1-XCO₂)], fₛ, fₗ, fᵪ, sOₑₓ, mOₑₓ, XCO₂, fₑₓ
    #     # end
    #     return SA[  fₛ - fₗ                                                         # Equation (1)
    #                 fₗ - fᵪ                                                         # Equation (2)
    #                 1 - sOₑₓ - mOₑₓ - fOₑᵪ(fH₂, fₛ, Mh, P, TOₑₓ) - Φ*XCO₂/(1-XCO₂)        # Equation (3)
    #                 P - COH_P(y₄, fₛ)],                                             # Equation (4) 
    #                 fₛ, fₗ, fᵪ, sOₑₓ, mOₑₓ, XCO₂, fH₂                               # Extras
    # end

    # function constrained_smOₑₓ_XCO₂_solver(y    :: SVector{4,Float64},                              # Initial transformed solution vector
    #                                     slim :: K, mlim :: K, clim :: K, minfH₂ :: K, maxfH₂ :: K,  # Limits for smOₑₓ, mOₑₓ, XCO₂ and fH₂
    #                                     sfO2, ∂Sᵢ,                                                  # Solid fO₂ and numerical derivative (MAGEMin)
    #                                     ϵ    :: K, damp :: K, niter :: Int64;                       # Solver tolerance, dampening factor, and maximum iterations
    #                                     # Verbose parameters
    #                                     verb_flag::Int64,P::K,T::K,ϕ::K,Rs::K,Rf::K,
    #                                     TOex::K,TOₑₓ::K,p::K,TC::K,Mt::K,H::K,Mh::K,
    #                                     # Pre-computed parameters
    #                                     s::K,Φ::K,Φₘ::K,lowclip::K,IDV::K,_ln10::K,_T::K,Ys1::K,Ys2::K,sharpness::K,a::K,
    #                                     idxO::Int64, plotevo::Bool,SymXox::SVector{N, Symbol},verbose::Bool,debugging::Bool,
    #                                     dummy::Vector{Float64}, molXB::SVector{N, Float64}) where{K<:AbstractFloat, N}
        
    #     # Residual - fO₂ - Partitioning matrix for plotting
    #     mat = zeros(niter, 4, 3) # [Residuals, fO₂, Partitioning]

    #     # y margins
    #     y_lows      = x_to_y(lowclip, lowclip, slim)
    #     y_highs     = x_to_y(slim-lowclip, lowclip, slim)
    #     y_lowm      = x_to_y(lowclip, lowclip, mlim)
    #     y_highm     = x_to_y(mlim-lowclip, lowclip, mlim)
    #     y_lowc      = x_to_y(lowclip, lowclip, clim)
    #     y_highc     = x_to_y(clim-lowclip, lowclip, clim)
    #     y_lowH2     = x_to_y(lowclip, minfH₂, maxfH₂)
    #     y_highH2    = x_to_y(maxfH₂-lowclip, minfH₂, maxfH₂)

    #     # Converged flag
    #     converged   = false

    #     itout       = 0
    #     Dsat        = false
    #     switch      = false
    #     # SM3       = @SMatrix zeros(3,3) # Static Matrix 3 × 3
    #     SM4         = @SMatrix zeros(4,4) # Static Matrix 4 × 4
    #     ac          = 1.0
    #     αᵢ          = 1.0

    #     for it in 1:niter
    #         # Evaluate current stage
    #         y₁, y₂, y₃, y₄ = y
    #         # Force boundary just after regime transition
    #         # Dsat && (y₃=Inf)
    #         # !Dsat && (y₄=-Inf)
    #         (switch&&Dsat)  && (y₃=Inf; y₄=x_to_y(1e-9, lowclip, fₑₓlim))
    #         (switch&&!Dsat) && (y₃=x_to_y(0.97clim, lowclip, clim); y₄=-Inf)
    #         switch = false
    #         # Iteration variables
    #         α   = evα(y_to_x(y₂, lowclip, mlim), Φₘ)
    #         θₘ  = evθₘ(α, s)
    #         # Compute residual
    #         params = (;sfO2, P, T, IDV, Mh, SymXox, idxO, _ln10, _T, s, Φ, Φₘ, molXB, sharpness, Dsat, lowclip, slim, mlim, clim, minfH₂, maxfH₂, TOₑₓ)
    #         Fx, fₛ, fₗ, fᵪ, sOₑₓ, mOₑₓ, XCO₂, fH₂ = Rx(y₁, y₂, y₃, y₄, dummy, params)
    #         # Store values
    #         mat[it,:,1] .= [Fx[1], Fx[2], Fx[3], Fx[4]]
    #         mat[it,1,2]  = fₛ
    #         mat[it,2,2]  = fₗ
    #         mat[it,3,2]  = fᵪ
    #         mat[it,:,3] .= [sOₑₓ, mOₑₓ, XCO₂, fH₂]
    #         # Check convergence
    #         aR = maximum(abs.(Fx))
    #         if aR<=ϵ
    #             converged = true
    #             itout=it
    #             if verbose
    #                 println("")
    #                 if verb_flag==-1
    #                     println("---- Solution (P=$(P)GPa | T=$(T)K | ϕ=$(ϕ) | Rs=$(Rs) | Rf=$(Rf) | Mix=$p | TCarbon=$TC | H₂O=($H)) ----")
    #                 else
    #                     println("---- Solution (P=$(P)GPa | T=$(T)K | ϕ=$(ϕ) | TOₑₓ=$TOex | Mix=$p | TC=$TC | H₂O=($H)) ----")
    #                 end
    #                 println("Shared fO₂ = $(sfO2(sOₑₓ)) |  residual = $aR")
    #                 println("Total Oₑₓ budget = $(round((1e2TOₑₓ/Mt), digits=4))% of total mass")
    #                 println("TOₑₓ partitioning → [$(round(1e2sOₑₓ, digits=4))% solid, $(round(1e2mOₑₓ, digits=4))% melt] stored as Fe₂O₃, [$(round(1e2fₑₓ, digits=4))% solid, $(round(1e2(1 - sOₑₓ - mOₑₓ - fₑₓ), digits=4))% melt] stored as CaCO₃ and CO₂ respectively")
    #                 println("Melt XCO₂ = $(XCO₂)")
    #                 println("Converged in $it iterations.")
    #             end
    #             return converged, mat, itout
    #         end
    #         # Jacobian inverse (Chain rule)
    #         ∂x∂y₁, ∂x∂y₂, ∂x∂y₃, ∂x∂y₄ = ∂x∂y(y₁, lowclip, slim), ∂x∂y(y₂, lowclip, mlim), ∂x∂y(y₃, lowclip, clim), ∂x∂y(y₄, minfH₂, maxfH₂)
    #         # Partial derivatives
    #         #   ∂(1)
    #             ∂1∂sOₑₓ = ∂Sᵢ(sOₑₓ)*∂x∂y₁
    #             ∂1∂mOₑₓ = -∂M∂mOₑₓ(Φₘ, Ys1, Ys2, α, θₘ, _ln10, a, molXB[3])*∂x∂y₁
    #             ∂1∂XCO₂ = 0.0; ∂1∂fH₂  = 0.0
    #         #   ∂(2)
    #             ∂2∂mOₑₓ = ∂1∂mOₑₓ*∂x∂y₂
    #             ∂2∂XCO₂ = sharpness==0.0 ? -∂C∂XCO₂_noedge(XCO₂, _ln10)*∂x∂y₂ : -∂C∂XCO₂(XCO₂, _ln10, sharpness, clim)*∂x∂y₂
    #             ∂2∂sOₑₓ = 0.0; ∂2∂fH₂  = 0.0
    #         #   ∂(3)
    #             ∂3∂sOₑₓ = (-∂1∂sOₑₓ - ∂fOₑᵪ∂fO₂(fH₂, fO₂, Mf, Ptot, TOₑₓ)*∂1∂sOₑₓ)*∂x∂y₃
    #             ∂3∂fH₂  = -∂fOₑᵪ∂fH₂(fH₂, fO₂, Mf, Ptot, TOₑₓ)*∂x∂y₃
    #             ∂3∂XCO₂ = -∂Stagno∂XCO₂(Φ, XCO₂)
    #             ∂3∂mOₑₓ = -1.0
    #         #   ∂(4)
    #             ∂4∂sOₑₓ = ∂COH_P∂fO₂(fH₂, fₛ)*∂1∂sOₑₓ*∂x∂y₄
    #             ∂4∂fH₂  = ∂COH_P∂fH₂(fH₂, fₛ)*∂x∂y₄
    #             ∂4∂mOₑₓ = 0.0; ∂4∂XCO₂ = 0.0
    #         if debugging
    #             @printf "Iteration %d: sOₑₓ = %.4f (%.4f), mOₑₓ = %.4f (%.4f), XCO₂ = %.4f (%.4f), fH₂ = %.4f (%.4f)\n" it sOₑₓ slim mOₑₓ mlim XCO₂ clim fH₂ maxfH₂
    #             @printf "\t(R₁=%.16f, R₂=%.16f, R₃=%.16f, R₄=%.16f)" Fx[1] Fx[2] Fx[3] Fx[4]
    #             @printf "  (sfO₂=%.4f, mfO₂=%f, cfO₂=%.4f)\n" fₛ fₗ fᵪ
    #             @printf "Jacobian:\n"
    #             @printf "[%.16f %.16f %.16f %.16f\n" ∂1∂sOₑₓ ∂1∂mOₑₓ ∂1∂XCO₂ ∂1∂fH₂
    #             @printf " %.16f %.16f %.16f %.16f\n" ∂2∂sOₑₓ ∂2∂mOₑₓ ∂2∂XCO₂ ∂2∂fH₂
    #             @printf " %.16f %.16f %.16f %.16f\n" ∂3∂sOₑₓ ∂3∂mOₑₓ ∂3∂XCO₂ ∂3∂fH₂
    #             @printf " %.16f %.16f %.16f %.16f]\n" ∂4∂sOₑₓ ∂4∂mOₑₓ ∂4∂XCO₂ ∂4∂fH₂
    #             @printf "  (∂x∂y₁=%.4f, ∂x∂y₂=%f, ∂x∂y₃=%.4f, ∂x∂y₄=%.4f)\n" ∂x∂y₁ ∂x∂y₂ ∂x∂y₃ ∂x∂y₄
    #         end
    #         # Jacobian
    #         J⁻¹  = SM4 + inv([  ∂1∂sOₑₓ ∂1∂mOₑₓ ∂1∂XCO₂ ∂1∂fH₂
    #                             ∂2∂sOₑₓ ∂2∂mOₑₓ ∂2∂XCO₂ ∂2∂fH₂
    #                             ∂3∂sOₑₓ ∂3∂mOₑₓ ∂3∂XCO₂ ∂3∂fH₂
    #                             ∂4∂sOₑₓ ∂4∂mOₑₓ ∂4∂XCO₂ ∂4∂fH₂])
    #         # Newton step
    #         dn = J⁻¹*Fx
            
    #         # Adaptive stepping
    #         (y[1]-dn[1]<y_lows)                && (α = min(α, 0.5*(y[1]-y_lows)/dn[1]))
    #         (y[1]-dn[1]>y_highs)               && (α = min(α, 0.5*(y[1]-y_highs)/dn[1]))
    #         (y[2]-dn[2]<y_lowm)                && (α = min(α, 0.5*(y[2]-y_lowm)/dn[2]))
    #         (y[2]-dn[2]>y_highm)               && (α = min(α, 0.5*(y[2]-y_highm)/dn[2]))
    #         # if Dsat
    #         #     (y[4]-dn[3]<y_lowfₑₓ)          && (switch=true; Dsat=false; ac=1)
    #         #     (debugging && !Dsat)             && println("\t • De-Saturating. Releasing XCO₂ constraint.")
    #         #     (y[4]-dn[3]>y_highfₑₓ)         && (α = min(α, 0.5*(y[4]-y_highfₑₓ)/dn[3]))
    #         #     # Take step
    #         #     α  = 0.7αᵢ*(ac/150)
    #         #     y = y - α*[dn[1], dn[2], 0.0, dn[3]]
    #         # else
    #             (y[3]-dn[3]<y_lowc)            && (α = min(α, 0.5*(y[3]-y_lowc)/dn[3]))
    #             (y[3]-dn[3]>y_highc)           && (α = min(α, 0.5*(y[3]-y_lowc)/dn[3]))
    #             (y[4]-dn[4]<y_lowc)            && (α = min(α, 0.5*(y[4]-y_lowH2)/dn[4]))
    #             (y[4]-dn[4]>y_highc)           && (α = min(α, 0.5*(y[4]-y_lowH2)/dn[4]))
    #             # (!Dsat && y[3]-dn[3]>y_highc)  && (switch=true; Dsat=true; ac=1)
    #             # (debugging && Dsat)            && println("\t • Saturating. Locking XCO₂ to $clim .")
    #             # Take step
    #             α  = αᵢ*(ac/150)
    #             y = y - α*[dn[1], dn[2], dn[3], dn[4]]
    #         # end
    #         debugging && @printf "\t(α = %.3f)\n\n" α
    #         ac += 1

    #         # Output if not converged
    #         if verbose && (it==niter)
    #             if verb_flag==-1
    #                 println("---- Solution (P=$(P)GPa | T=$(T)K | ϕ=$(ϕ) | Rs=$(Rs) | Rf=$(Rf) | Mix=$p | TCarbon=$TC | H₂O=($H)) ----")
    #             else
    #                 println("---- Solution (P=$(P)GPa | T=$(T)K | ϕ=$(ϕ) | TOₑₓ=$TOex | Mix=$p | TC=$TC | H₂O=($H)) ----")
    #             end
    #             println("Shared fO₂ = $(sfO2(sOₑₓ)) |  residual = $aR")
    #             println("Total Oₑₓ budget = $(round((1e2TOₑₓ/Mt), digits=4))% of total mass")
    #             println("TOₑₓ partitioning → [$(round(1e2sOₑₓ, digits=4))% solid, $(round(1e2mOₑₓ, digits=4))% melt] stored as Fe₂O₃, [$(round(1e2fₑₓ, digits=4))% solid, $(round(1e2(1 - sOₑₓ - mOₑₓ - fₑₓ), digits=4))% melt] stored as CaCO₃ and CO₂ respectively")
    #             println("Melt XCO₂ = $(XCO₂)")
    #             println("Did not converged in $it iterations.")
    #         end
    #     end
    #     return false, mat, niter
        
    # end

    # # Variable transformations
    # @inline y_to_x(y, minV, maxV) = minV + (maxV - minV) / (1.0 + exp(-y))
    # @inline x_to_y(x, minV, maxV) = -log((maxV - x)/(x + minV))
    # @inline ∂x∂y(y, minV, maxV)   = (maxV - minV)*exp(-y)/(1.0+exp(-y))^2
    # # Variables carbon
    # @inline evΦ(TOₑₓ, molMf)   = 2*mm.O*molMf/TOₑₓ   # Conversion factor for XCO₂
    # # Variables melt
    # @inline evΦₘ(TOₑₓ, Mf)   = TOₑₓ/Mf/mm.O          # Conversion factor for melt
    # @inline evα(mOₑₓ, Φₘ)    = Φₘ*mOₑₓ               # Mass fraction of TOₑₓ in melt → non-normalized mass of Oₑₓ in melt
    # @inline evθₘ(α, s)       = s - α                 # Molar normalization factor
    # # Partial derivatives
    # @inline ∂S∂sOₑₓ(sfO2, sOlist)                        = @views diff(sfO2)./diff(sOlist)
    # @inline ∂C∂XCO₂(XCO₂, _ln10, sharpness, maxV)        = _ln10*(1/XCO₂ + sharpness/(maxV-XCO₂))
    # @inline ∂C∂XCO₂_noedge(XCO₂, _ln10)                  = _ln10*(1/XCO₂)
    # @inline ∂M∂mOₑₓ(Φₘ, Ys1, Ys2, α, θₘ, _ln10, a, XFeO) = -Φₘ/a * ( θₘ^(-2)*(Ys1 + 2Ys2/θₘ) - _ln10*(1/α + 2/(XFeO - 2α)))
    # @inline ∂Stagno∂XCO₂(Φ, XCO₂)                        = Φ*(1/(1-XCO₂)^2)
    # function COH_P(fH₂, fO₂)
    #     # Equilibrium constants
    #         K₁, K₂, K₃, K₄ = 1.0, 1.0, 1.0, 1.0
    #     # Ideality coefficients
    #         γ₁, γ₂, γ₃, γ₄, γ₅, γ₆ = 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
    #     # Carbon activity
    #         aC = 1.0
    #     # Unpack fO₂ and fH₂
    #         ufO₂ = 10^fO₂
    #         ufH₂ = 10^fH₂
    #     # Pₜₒₜ = ∑ fᵢ/γᵢ
    #         return K₁*(√ufO₂*ufH₂)/γ₁ + K₂*(aC*ufH₂^2)/γ₂ + K₃*(aC*ufO₂)/γ₃ + K₄*(aC*√ufO₂)/γ₄ + ufO₂/γ₅ + ufH₂/γ₆
    #     #               |                   |                 |                 |              |        |
    #     #             fH₂O                 fCH₄              fCO₂              fCO            fO₂      fH₂   
    # end
    # function fOₑᵪ(fH₂, fO₂, Mh, Ptot, TOₑₓ)
    #     # Equilibrium constants
    #         K₁, K₃, K₄ = 1.0, 1.0, 1.0
    #     # Ideality coefficients
    #         γ₁, γ₃, γ₄, γ₅ = 1.0, 1.0, 1.0, 1.0
    #     # Carbon activity
    #         aC = 1.0
    #     # Unpack fO₂
    #         ufO₂ = 10^fO₂
    #         ufH₂ = 10^fH₂
    #     # Compute MMf
    #         X₁ = K₁*ufH₂*√ufO₂/γ₁ # H₂O
    #         X₃ = K₃*aC*ufO₂/γ₃   # CO₂
    #         X₄ = K₄*aC*√ufO₂/γ₄  # CO
    #         X₅ = ufO₂/γ₅         # O₂
    #         MMf = X₁*mm.H2O + X₃*mm.CO2 + X₄*(mm.C+mm.O) + X₅*2mm.O
    #     # Compute COH_Oₑₓ_mass / TOₑₓ
    #         return (1e-3*Mh/MMf/Ptot*(X₁ + 2X₃ + X₄ + X₅)*mm.O) / TOₑₓ
    # end
    # ∂COH_P∂fH₂(fH₂, fO₂) = ForwardDiff.derivative(fH₂ -> COH_P(fH₂, fO₂), fH₂)
    # ∂COH_P∂fO₂(fH₂, fO₂) = ForwardDiff.derivative(fO₂ -> COH_P(fH₂, fO₂), fO₂)
    # ∂fOₑᵪ∂fH₂(fH₂, fO₂, Mf, Ptot, TOₑₓ) = ForwardDiff.derivative(fH₂ -> fOₑᵪ(fH₂, fO₂, Mf, Ptot, TOₑₓ), fH₂)
    # ∂fOₑᵪ∂fO₂(fH₂, fO₂, Mf, Ptot, TOₑₓ) = ForwardDiff.derivative(fO₂ -> fOₑᵪ(fH₂, fO₂, Mf, Ptot, TOₑₓ), fO₂)

    # Dampening
    # function bt_line_search(Δx, y₁, y₂, y₃, y₄, r, dummy, params; α = 1.0, ρ = 0.5, lstol = 0.9, α_min = 1.0e-8)

    #     x = SA[y₁, y₂, y₃, y₄]
    #     perturbed_x = @. x + α * Δx
    #     r, = Rx(x.data..., dummy, params)
    #     rnorm = mynorm(r, x)

    #     # Iterate unless step length becomes too small
    #     while α > α_min
    #         # Apply scaled update
    #         perturbed_x = @. x + α * Δx

    #         # Get updated residual
    #         perturbed_r, = Rx(perturbed_x.data..., dummy, params)
    #         perturbed_rnorm = mynorm(perturbed_r, x)

    #         # Check whether residual is sufficiently reduced
    #         if perturbed_rnorm ≤ lstol * rnorm
    #             break
    #         end

    #         # Bisect step length
    #         α *= ρ
    #     end

    #     return α
    # end
    # @generated function mynorm(x::SVector{N, T}, y::SVector{N}) where {N, T}
    #     return quote
    #         @inline
    #         v = zero(T)
    #         Base.@nexprs $N i -> begin
    #             xi = @inbounds x[i]
    #             yi = @inbounds y[i]
    #             v += !iszero(yi) * abs(xi / yi)
    #         end
    #         return v
    #     end
    # end

#

# ==============================================================
# ===== Silicate melt density (Jing and Karato 2011, 2012) =====
# ==============================================================
#  Audited 2026-09-10 against sources: TOex_math.txt PART 15 has the equations, the numbers and
#  what was wrong. The hard-sphere law is kept as calibrated; the anhydrous SOLID below now comes
#  from MAGEMin's rho with the full sb24 oxide list, and the tables are written basalt-first.
    function ev_P(σᵢ₀, T, Tᵣ, ηᵢ, ξᵢ, α₀, V, V₀, X)
        σᵢ      = @. σᵢ₀ * (T/Tᵣ)^(ηᵢ) * exp( ξᵢ/3 * α₀ * (T - Tᵣ) ) * (V/V₀)^(ξᵢ/3) # cm (Hard-sphere diameter)
        Vₘᵢ     = Nₐ * π / 6 * σᵢ.^3 .* X    # cm³ (Component monatomic volume)
        Vₘᵢ₀    = Nₐ * π / 6 * σᵢ₀.^3 .* X   # cm³ (Component reference monatomic volume)
        Vₘ      = sum(Vₘᵢ)                   # cm³ (Total monatomic volume)
        Vₘ₀     = sum(Vₘᵢ₀)                  # cm³ (Total reference monatomic volume)
        ξ       = sum(ξᵢ.*Vₘᵢ)/Vₘ            # Deformability of spheres
        f       = Vₘ/V                       # Packing fraction of spheres
        f₀      = Vₘ₀/V₀                     # Reference packing fraction of spheres
        Φ       = (1 + f + f^2)/(1-f)^3      # Deviation of ideality due to packing
        Φ₀      = (1 + f₀ + f₀^2)/(1-f₀)^3   # Deviation of ideality due to packing at reference
        P_calc  = 1e-3R*T/V * ( (1-ξ)*Φ - Φ₀*(V₀/V)^(4/3 - 1) + ξ*Φ₀*(Vₘ₀/Vₘ)^(8/3 - 1) ) # GPa
        return P_calc
    end
    ev_∂P∂V(σᵢ₀, T, Tᵣ, ηᵢ, ξᵢ, α₀, V, V₀, X) = ForwardDiff.derivative(V -> ev_P(σᵢ₀, T, Tᵣ, ηᵢ, ξᵢ, α₀, V, V₀, X), V)

    # α₀ is the law's OWN 1-bar expansivity parameter: 2.5e-4 gives α(1 bar) = 9.2e-5/K against Lange &
    #  Carmichael's 10.2e-5/K for a peridotitic melt; handing it the solid's α (~3e-5) turns α NEGATIVE.
    function Jing_Karato_silicate_melt_density(P, T, Xin, Xoxin; reg=1, niter=50, α₀ = 2.5e-4, verbose=false)

        # Force Xox list (Na2O, K2O, TiO2 have no hard-sphere parameters here and are renormalised away)
        Xox = ["SiO2", "Al2O3", "FeO", "MgO", "CaO", "H2O"]
        X   = zeros(length(Xox))
        MM  = get_Xoxmm(Xox)
        X   .= filter_ox_list(Xin, Xoxin, Xox)

        # Parameters
        mma = X[1]*mm.SiO2 + X[2]*mm.Al2O3 + X[3]*mm.FeO + X[4]*mm.MgO + X[5]*mm.CaO + X[6]*mm.H2O # g/mol (Molar mass of the mixture)
        Tᵣ  = 1673.          # K (Reference temperature)
        V₀ᵢ = SA[26.86, 37.11, 13.65, 11.69, 16.57, 22.90] # cm³/mol (Reference liquid component volume)
        V₀  = sum(V₀ᵢ.*X)    # cm³/mol (Reference liquid volume)
        V   = 10.0           # cm³/mol (Initial V guess)

        # Regression parameters [SiO2, Al2O3, FeO, MgO, CaO, H2O]
        if reg==1
            σᵢ₀ = 1e-7SA[0.369, 0.326, 0.286, 0.272, 0.315, 0.172] # cm (Reference hard-sphere diameters)
            ηᵢ  = SA[-0.04, -0.03,  0.01,  0.03, -0.07,  0.06]
            ξᵢ  = SA[ 0.70,  0.58,  0.17,  0.12,  0.09, -0.35]
        elseif reg==2
            σᵢ₀ = 1e-7SA[0.369, 0.326, 0.285, 0.272, 0.315, 0.180]
            ηᵢ  = SA[-0.04, -0.03,  0.01,  0.04, -0.06,  0.13]
            ξᵢ  = SA[ 0.69,  0.58,  0.15,  0.14,  0.07,   0.0]
        end

        # P(V) falls monotonically with V and diverges at the packing wall, so a root always exists:
        #  bracket it walking down from 1.5 V₀ in 1 % steps, bisect, then polish with Newton. A bare
        #  Newton from V = 10 cm³/mol failed at cold, water-rich, high-P corners (50 GPa, 500 K, 8 wt%).
        Pf(V) = ev_P(σᵢ₀, T, Tᵣ, ηᵢ, ξᵢ, α₀, V, V₀, X)
        hi = 1.5V₀; lo = hi; plast = -Inf
        while true
            lo = 0.99lo
            p  = Pf(lo)
            (isfinite(p) && p >= P) && break
            # P(V) has a tension minimum above the 1-bar volume; only once positive must it keep rising
            (!isfinite(p) || (plast > 0 && p < plast) || lo < 1e-2V₀) && error("Jing_Karato_silicate_melt_density: no bracket at P=$P GPa, T=$T K")
            plast = p
        end
        hi = lo/0.99
        for _ in 1:80
            m = 0.5(lo+hi); (Pf(m) > P) ? (lo = m) : (hi = m)
            (hi - lo) < 1e-9hi && break
        end
        V = 0.5(lo+hi)
        for it in 1:niter
            dP = Pf(V) - P
            verbose && (@printf "Iteration %d: V = %.4f cm³/mol, ρ = %.4f g/cm³, P_calc = %.4f GPa, P_target = %.4f GPa\n" it V mma/V P+dP P)
            (abs(dP) < 1e-4) && break
            Vn = V - dP/ev_∂P∂V(σᵢ₀, T, Tᵣ, ηᵢ, ξᵢ, α₀, V, V₀, X)
            (isfinite(Vn) && lo < Vn < hi) || break   # keep the bracketed bisection answer
            V = Vn
        end
        return mma/V
    end

    # Inverse 3rd-order Birch-Murnaghan: the V at which P(V) = P, by bisection.
    function bm3_volume(P, V₀, K₀, K′; P₀=0.0)
        f(V) = 1.5K₀*((V₀/V)^(7/3) - (V₀/V)^(5/3))*(1 + 0.75(K′-4)*((V₀/V)^(2/3) - 1)) - (P - P₀)
        lo, hi = 0.05V₀, 1.5V₀
        for _ in 1:200; m = 0.5(lo+hi); (f(m)*f(lo) < 0) ? (hi = m) : (lo = m); end
        return 0.5(lo+hi)
    end
    # Partial molar volume of H₂O in the melt for melt_H2O_model=:ideal (math 15.4): the FPMD compression of
    #  Karki, Ghosh & Karato (2021; V₀ 25.0 cm³/mol, K₀ 1.91 GPa, K′ 5.80 at 2000 K) scaled to the law's own
    #  1-bar value V₁, so the two melt laws agree where the hard-sphere fit is calibrated and part above it.
    melt_H2O_partial_volume(P, V₁) = V₁*bm3_volume(P, 25.0, 1.91, 5.80)/25.0

# ======================================================
# ===== Solid H₂O density correction (Gerya, 2004) =====
# ======================================================
    xlogx(x) = x > 0 ? x*log(x) : zero(x)   # 0·log 0 = 0, so Xliq = 1 is a finite point of Gₛ
    function Gₛ(P, T, Xliq)
        # Gₛ Parameters (Gerya, Podlesskii, Perchuk & Maresch 2004, H₂O fluid; P in bar, T in K)
            H298             = -286831.56               # J
            S298             = 65.188                   # J K⁻¹
            Vₛ               = 1.71382                  # J bar⁻¹
            ϕ                = 6209                     # bar
            c₁, c₂, c₃       = 7.23576, 0.31482, 0.0
            ΔHₛ₁, ΔHₛ₂, ΔHₛ₃ = 4586.46, 0.0, 0.0        # J
            ΔVₛ₁, ΔVₛ₂, ΔVₛ₃ = 0.04310884, 0.0, 0.0     # J bar⁻¹
            ΔH⁰ₒᵣ            = -44838.80                # J
            ΔS⁰ₒᵣ            = -122.397                 # J K⁻¹
            ΔCₚ⁰ₒᵣ           = 21.486                   # J K⁻¹
            ΔV⁰ₒᵣ            = 0.0                      # J bar⁻¹
            Wh₁              = -28793.19                # J⁻¹
            Ws₁              = -11.704                  # J K⁻¹
            Wcp₁             = 5.086                    # J K⁻¹
            Wv₁              = 0.0                      # J bar⁻¹
            ΔHₛλ⁰            = 0.0                      # J 
            ΔVₛλ⁰            = 0.0                      # J bar⁻¹
        # Reference parameters
            T₀ = 298.15 # K
            P₀ = 1      # bar
            n  = 2
        # Auxilliaries
            Ψ   = (5/4)*(P₀ + ϕ)^(1/5)*((P + ϕ)^(4/5) - (1 + ϕ)^(4/5))
            e₁  = exp(-(ΔHₛ₁ + ΔVₛ₁*Ψ)/R/T)
            e₂  = exp(-ΔHₛ₁/R/T)
            e₀  = exp(-ΔHₛ₁/R/T₀)
        # Gibbs free energy (J/mol)
            G   = H298 - T*S298 + Vₛ*Ψ + R*T*(c₁*log(1-e₁)+c₂*log(1-e₂)) -
                    (c₁ + c₂)*(ΔHₛ₁*(1-T/T₀)*e₀/(1-e₀) + R*T*log(1-e₀)) +
                    R*T*(xlogx(1-Xliq)+xlogx(Xliq)) +
                    (1-Xliq)*R*T*log(ϕ*(Xliq^2)+P) -
                    (1-Xliq)*(ΔH⁰ₒᵣ - T*ΔS⁰ₒᵣ + ΔCₚ⁰ₒᵣ*(T - T₀ - T*log(T/T₀))) +
                    (Wh₁ - T*Ws₁ + Wcp₁*(T - T₀ - T*log(T/T₀)))*Xliq*(1-Xliq)
            return G
    end

    ∂Gₛ∂P(P,T,Xliq)     = ForwardDiff.derivative(P -> Gₛ(P, T, Xliq), P)
    ∂Gₛ∂Xliq(P,T,Xliq)  = ForwardDiff.derivative(Xliq -> Gₛ(P, T, Xliq), Xliq)

    # Ordering parameter that minimises Gₛ: coarse scan, then bisection on ∂G/∂Xliq.
    function Gₛ_minXliq(P, T)
        xs = LinRange(1e-4, 1-1e-4, 200); i = argmin(Gₛ.(P, T, xs))
        lo, hi = xs[max(i-1,1)], xs[min(i+1,length(xs))]
        for _ in 1:50; m = 0.5(lo+hi); (∂Gₛ∂Xliq(P, T, m) > 0) ? (hi = m) : (lo = m); end
        return 0.5(lo+hi)
    end

    # Molar volume of H₂O, cm³/mol, from V = ∂G/∂P. :condensed is the liquid-like branch (Xliq = 1) used for
    #  water BOUND in minerals - it cannot boil off at low P; :equilibrium minimises G and reproduces Zhang &
    #  Duan (2009) to 0.6 % mean over 0.2-10 GPa (math 15.3), but that is the volume of a FREE fluid.
    function Gerya_H2O_volume(P_GPa, T; branch=:condensed)
        Pbar = 1e4P_GPa
        Xl   = branch == :condensed ? 1.0 : Gₛ_minXliq(Pbar, T)
        return 10∂Gₛ∂P(Pbar, T, Xl)   # 10 J bar⁻¹ mol⁻¹ = cm³ mol⁻¹
    end

    # Xin   = mole fractions on Xox, LAST entry H2O (0 for the anhydrous solid)
    # Xox   = ["SiO2", "MgO", "FeO", "CaO", "Al2O3", "Na2O", "Cr2O3", "O", "H2O"]
    function Gerya_solid_H2O_density_correction_interpolator(Xin,Xox; nP=100, nT=100, out=nothing, verbose=false,
                                                             Prange=(0.1, 75.), Trange=(500., 3000.), H2O_branch=:condensed)
        # 1. Anhydrous DENSITY of the solid from MAGEMin sb24 (out.rho - out.V is on another mole basis).
        # 2. Molar volume of H₂O from the Gibbs energy of Gerya et al. (2004): V_H2O = ∂G/∂P.
        # 3. Hydrated volume by ideal mixing, V = (1-x)·M_anh/ρ_anh + x·V_H2O  (math 15.2).
            @assert Xox[end]=="H2O" "Last oxide entry must be H2O"
            @assert ("Na2O" in Xox) && ("Cr2O3" in Xox) "sb24 needs Na2O and Cr2O3 in the list (absent, MAGEMin returns 99% akimotoite at 5.2 g/cm³)"
            X = sum(Xin)>1.0 ? 1e-2Xin : Xin
            X = X./sum(X)
            noWX    = X[1:end-1]
            @assert noWX[findfirst(Xox.=="Na2O")] > 0 && noWX[findfirst(Xox.=="Cr2O3")] > 0 "sb24 needs Na2O and Cr2O3 > 0 (a trace is enough)"
        # Vectors
            nPnT    = nP*nT
            P       = LinRange(Prange[1], Prange[2], nP)
            T       = LinRange(Trange[1], Trange[2], nT)
        # Vectorize arguments
            Xv = Vector{Vector{Float64}}(undef, nPnT)
            for i in 1:nPnT
                Xv[i] = noWX
            end
            vP, vT = repeat(P, outer=nT), repeat(T, inner=nP)
        # Map
            Vh2o, ρanh, α = zeros(nP, nT), zeros(nPnT), zeros(nPnT)
        # Compute H₂O molar volume
            for (ip, p) in enumerate(P)
                for (it, t) in enumerate(T)
                    Vh2o[ip, it] = Gerya_H2O_volume(p, t; branch=H2O_branch)
                end 
            end
        # Retrieve anhydrous solid density
            if isnothing(out)
                data   = Initialize_MAGEMin("sb24", verbose=false);
                out    = multi_point_minimization(10vP, vT.-273.15, data, X=Xv, Xoxides=Xox[1:end-1], sys_in="mol", name_solvus=true, progressbar=verbose) # kbar and K
                Finalize_MAGEMin(data)
            end
            # Extract information
            nbad = 0
            for i in 1:nPnT
                ρanh[i] = out[i].rho/1e3    # g/cm³
                α[i]    = out[i].alpha[1]   # 1/K
                (ρanh[i] < 2.4 || ρanh[i] > 6.0) && (nbad += 1)
            end
            nbad > 0 && @warn "Gerya_solid_H2O_density_correction_interpolator: $nbad of $nPnT MAGEMin points have rho outside 2.4-6.0 g/cm³ - check the oxide list"
        # Reshape
            ρanh = reshape(ρanh, nP, nT)
            α = reshape(α, nP, nT)

        # Construct anhydrous ρ | H₂O V | α interpolator objects
            VVρ = [ Interpolations.extrapolate(Interpolations.interpolate((P,T),ρanh,Gridded(Linear())), Flat()), 
                    Interpolations.extrapolate(Interpolations.interpolate((P,T),Vh2o,Gridded(Linear())), Flat()), 
                    Interpolations.extrapolate(Interpolations.interpolate((P,T),α,Gridded(Linear())), Flat())]
            return VVρ

    end

    function Gerya_solid_H2O_density_correction(P::Float64,T::Float64,Xin::Vector{Float64},Xox::Vector{String}; nP=50, nT=30, VVρ=nothing, MM=nothing, verbose=false)
        # Construct interpolator if not passed in
            isnothing(VVρ) && (VVρ = Gerya_solid_H2O_density_correction_interpolator(Xin,Xox; nP=nP, nT=nT, verbose=verbose))
            isnothing(MM)  && (MM = get_Xoxmm(Xox))
            ρA = VVρ[1](P,T)     # g/cm³, anhydrous solid
            VH = VVρ[2](P,T)     # cm³/mol, H₂O
        # Mole fractions with and without water
            X       = sum(Xin)>1.0 ? 1e-2Xin : Xin
            XH      = X./sum(X)
            XA      = @views X[1:end-1]./sum(X[1:end-1])
            massA   = sum(XA.*(@view MM[1:end-1]))
            massH   = sum(XH.*MM)
        # Ideal mixing of the anhydrous molar volume and the H₂O molar volume
            VA      = massA/ρA
            rV      = (1-XH[end])*VA + XH[end]*VH
        # ρ decrease in %
            return 1e2(1 - (massH/rV)/ρA)
    end

# =============================
# ===== Table constructor =====
# =============================
    # Xin   = mole fractions on Xox, last entry H2O = 0;  Xox = [..., "H2O"]
    # Returns SΔρ, MΔρ [nP × nT × nH × slot] as % density DECREASE relative to dry, P (GPa), T (K), Hwt (wt% H₂O).
    #  slot 1 = basalt, slot 2 = harzburgite (default=true) - the order thermochem_int weights by (basalt, harz).
    #  Hwt is UNIFORM in wt%: StagYY rebuilds the axis from (min, max, n) alone.
    #  melt_H2O_model = :JK (hard-sphere H₂O sphere, floored at 0 when melt_floor) | :ideal (math 15.4).
    function PT_H2O_ρ(Pi,Pf,Ti,Tf,Xin,Xox; nH=50, nP=50, nT=50, default=false, verbose=true,
                      Hmin_wt=0.0, Hmax_wt=8.0, melt_H2O_model=:JK, melt_floor=true, reg=1, H2O_branch=:condensed)

        # Default --> the sf_H2O compositions (Stixrude & Lithgow-Bertelloni 2024 basalt | harzburgite), wt%
            if default
                Xox  = ["SiO2", "MgO", "FeO", "CaO", "Al2O3", "Na2O", "Cr2O3", "O", "H2O"]
                XB   = [50.42, 9.77, 7.1, 12.54, 16.8, 2.23, 0.07, 0.0, 0.0]   # basalt
                XH   = [43.43, 45.93, 8.34, 0.9, 1.0, 0.01, 0.3, 0.0, 0.0]     # harzburgite
                Xin  = filter_ox_list(XB, Xox, Xox, wt_in=true, wt_out=false)  # mole fractions
                XinH = filter_ox_list(XH, Xox, Xox, wt_in=true, wt_out=false)
            else
                # Start Anhydrous
                @assert Xox[end]=="H2O" "Last oxide entry must be H2O"
                @assert Xin[end]==0.0 "Last oxide entry must be 0.0 (H₂O)"
            end
            melt_H2O_model in (:JK, :ideal) || error("PT_H2O_ρ: melt_H2O_model must be :JK or :ideal")
        # Grid
            P   = LinRange(Pi, Pf, nP)
            T   = LinRange(Ti, Tf, nT)
            Hwt = collect(LinRange(Hmin_wt, Hmax_wt, nH))
        # Result matrices
            comps = default ? (Xin, XinH) : (Xin,)
            nslot = length(comps)
            SΔρ  = zeros(nP, nT, nH, nslot)
            MΔρ  = zeros(nP, nT, nH, nslot)
            anhM = zeros(nP, nT, nslot)
        # Ordered MM vector
            MM = get_Xoxmm(Xox)
        # Initiate constructor
            for (is, Xc) in enumerate(comps)
                Xl = sum(Xc)>1.0 ? 1e-2Xc : copy(Xc)
                Xl = Xl./sum(Xl)
                MA = sum(Xl.*MM)                         # g per mole of oxides, anhydrous
                verbose && (@printf "Slot %d: anhydrous solid (MAGEMin sb24) and H₂O volumes...\n" is)
                VVρ = Gerya_solid_H2O_density_correction_interpolator(Xl, Xox; nP=nP, nT=nT, verbose=verbose,
                                                                      Prange=(max(Pi,1e-3), Pf), Trange=(Ti, Tf), H2O_branch=H2O_branch)
                # Anhydrous melt
                for it in 1:nT, ip in 1:nP
                    anhM[ip, it, is] = Jing_Karato_silicate_melt_density(P[ip], T[it], Xl, Xox; reg=reg)
                end
                for iH in 1:nH
                    verbose && (@printf "Slot %d: iH = %d/%d...\n" is iH nH)
                    w  = Hwt[iH]
                    x  = (w/MM[end])/((w/MM[end]) + (100-w)/MA)   # mole fraction H₂O
                    X  = (1-x).*Xl;  X[end] = x
                    MH = (1-x)*MA + x*MM[end]
                    # 1-bar partial molar volume of H₂O from the law itself, per T (used by :ideal)
                    V₁ = zeros(nT)
                    if melt_H2O_model == :ideal && x > 0.0
                        for it in 1:nT
                            ρd1 = Jing_Karato_silicate_melt_density(1e-4, T[it], Xl, Xox; reg=reg)
                            ρh1 = Jing_Karato_silicate_melt_density(1e-4, T[it], X,  Xox; reg=reg)
                            V₁[it] = (MH/ρh1 - (1-x)*MA/ρd1)/x
                        end
                    end
                    x == 0.0 && continue      # dry node: both corrections are exactly zero
                    for it in 1:nT, ip in 1:nP
                        # Solid H₂O density correction (% decrease)
                        SΔρ[ip, it, iH, is] = Gerya_solid_H2O_density_correction(P[ip], T[it], X, Xox; VVρ=VVρ, MM=MM)
                        # Hydrous melt density
                        ρd = anhM[ip, it, is]
                        if melt_H2O_model == :JK
                            ρh = Jing_Karato_silicate_melt_density(P[ip], T[it], X, Xox; reg=reg)
                        else
                            ρh = MH/((1-x)*MA/ρd + x*melt_H2O_partial_volume(P[ip], V₁[it]))
                        end
                        d = 1e2*(1 - ρh/ρd)
                        MΔρ[ip, it, iH, is] = melt_floor ? max(d, 0.0) : d
                    end
                end
            end
            return SΔρ, MΔρ, collect(P), collect(T), Hwt
    end
