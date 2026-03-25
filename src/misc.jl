# Stagno and Frost XCO2 equilibrium (2 Oex per 1 CO₂) (Melt C⁰ -> C⁴⁺ ==> fO₂ to XCO₂ mapping)
function eq_XCO2(logfO2, Pin, Tin)
    P = Pin >= 11. ? 11. : Pin
    T = Tin >= c2k(1600.) ? c2k(1600.) : Tin
    return max(min(10^(logfO2 - 5.44 + 21380/T - 0.078(1e-4P-1)/T), 1.0), 0.0)
end

Hirsc(P, T, X, T₀, ΔCₚ, a, b, c, y1, y2, y3, y4, y5, y6, y7, y8, y9, IDV) = (log10(X["Fe2O3"]/(X["FeO"])) - b - c/T + (ΔCₚ/R/log(10) * (1 - T₀/T - log(T/T₀))) + IDV/(1e-3R)/T/log(10) 
                        - (1/T)*(y1*X["SiO2"] + y2*0.0 + y3*X["MgO"] + y4*X["CaO"] + y5*X["Na2O"] + y6*0.0 + y7*0.0 + y8*X["SiO2"]*X["Al2O3"] + y9*X["SiO2"]*X["MgO"]))/a
function reverseHirsch!(logfO2, P, T, Xc, IDV)
    a=0.1917; b=-1.961; c=4158.1; ΔCₚ=33.25; T₀=1673.15
    y1=-520.46; y2=-185.37; y3=494.39; y4=1838.34; y5=2888.48; y6=3473.68; y7=-4473.6; y8=-1245.09; y9=-1156.86
    if Xc["FeO"] > 0.0
        Femax = 0.5Xc["FeO"]; Femin = 1e-8; Feguess = 0.5(Femax+Femin); residual=0.0
        for n in 1:10
            # Initial guess
            Xcc = copy(Xc)
            Xcc["FeO"] -= 2Feguess; Xcc["Fe2O3"] = Feguess
            # Model
            guess = Hirsc(P, T, Xcc, T₀, ΔCₚ, a, b, c, y1, y2, y3, y4, y5, y6, y7, y8, y9, IDV)
            # Convergence
            residual = abs(guess - logfO2)
            (residual < 1e-6) && break
            (guess > logfO2) ? (Femax = Feguess) : (Femin = Feguess)
            Feguess = 0.5(Femax + Femin)
        end
        Xc["FeO"] -= 2Feguess; Xc["Fe2O3"] = Feguess
        # Construct copy
        Xc = Dict(k => v/sum(values(Xc)) for (k,v) in Xc)
        return residual
    else
        # nothing
    end
end

# Δ function
function ΔR(logfO2, P, T, X, IDV)
    ns = length(logfO2)
    if ns > 1
        Δ = zeros(ns, 3)
        for i in 1:ns
            # Copy composition
            Xc = copy(X)
            # Residual
            Δ[i, 3] = reverseHirsch!(logfO2[i], P, T, Xc, IDV) # Moves XFeO -> XFe₂O₃ to accomodate the required log(XFe₂O₃/XFeO) to explain the selected fO₂
            # XFe₂O₃
            Δ[i, 1] = Xc["FeO"] 
            Δ[i, 2] = Xc["Fe2O3"] 
        end; return Δ
    else
        # Copy composition
        Xc = copy(X)
        d = reverseHirsch!(logfO2, P, T, Xc, IDV) # Moves XFeO -> XFe₂O₃ to accomodate the required log(XFe₂O₃/XFeO) to explain the selected fO₂
        return [Xc["FeO"], Xc["Fe2O3"], d]
    end
end

function Hirschmann_fO2_to_R(; RF=0.02, fO2range=(-25., -2), P=3.5, T=1400., FMQ=false, ns=60, fsize=(1400, 800), savein="", p=0.2)

    # IDV = ∫ΔVdP
    IDV = solve_∫ΔVdP([P-0.05P, P, P+0.05P],[T-0.05T, T, T+0.05T])[2,2,1]

    # Generate fO2 and XCO₂ arrays
    logfO2 = LinRange(fO2range..., ns)
    logfO2F = LinRange(fO2range..., 3ns)

    # Composition dictionary
    XH = oxidize_bulk(Dict( "SiO2" => 0.4343, "MgO" => 0.4593, "FeO" => 0.0834, "CaO" => 0.0090, "Al2O3" => 0.0100, "Na2O" => 0.0001, "Cr2O3" => 0.0030, "Fe2O3" => 0.0), RF, frac=true, FeFormat="FeO_Fe2O3", dict=true) # Harzburgite
    XBu = oxidize_bulk(Dict( "SiO2" => 0.5042, "MgO" => 0.0977, "FeO" => 0.0710, "CaO" => 0.1254, "Al2O3" => 0.1680, "Na2O" => 0.0223, "Cr2O3" => 0.0007, "Fe2O3" => 0.0), 0.0, frac=true, FeFormat="FeO_Fe2O3", dict=true) # MORB
    XBo = oxidize_bulk(Dict( "SiO2" => 0.5042, "MgO" => 0.0977, "FeO" => 0.0710, "CaO" => 0.1254, "Al2O3" => 0.1680, "Na2O" => 0.0223, "Cr2O3" => 0.0007, "Fe2O3" => 0.0), RF, frac=true, FeFormat="FeO_Fe2O3", dict=true) # MORB
    
    # Compose parent tracer bulk composition
    Xs = Dict{String, Float64}(); [Xs[k] = p*XBo[k] + (1-p)*XH[k] for k in keys(XH)]
    Xm = copy(XBu)
    # Parent fO₂ + FMQ buffer
    data = Initialize_MAGEMin("sb24", verbose=false);
    out = single_point_minimization(10P, T-273.15, data, X=collect(values(Xs)), Xoxides=collect(keys(Xs)), name_solvus=true);
    Finalize_MAGEMin(data); ΔFMQ = out.fO2 - out.dQFM
    logfO2FMQ = logfO2 .- ΔFMQ
    logfO2FFMQ = logfO2F .- ΔFMQ

    # Compute equilibrated XFeO | XFe₂O₃ | Residual
    Δ = ΔR(logfO2, P, T, Xm, IDV)

    # Plot
    fig = Figure(size=fsize)
    nf = 8;

    # ==== Iron =====
        Fe2clr = :firebrick2
        Fe3clr = :dodgerblue2
        XCOclr = :green
        # ==== Solution space ====
            Fe2, Fe3, res = ΔR(out.fO2, P, T, Xm, IDV)
            XCO2spec = eq_XCO2.(logfO2F, P, T)
            XCO2specsol = eq_XCO2(out.fO2, P, T)
            ax = Axis(fig[2,1], xlabel=FMQ ? L"fO_2\;(\Delta FMQ)" : L"fO_2\;\mathrm{[abs]}", ylabel=L"X\;\mathrm{[mol\;fr]}", xlabelsize=20, ylabelsize=20, yticklabelsize=18, xticklabelsize=18, ygridvisible=false, xgridvisible=false)
            ax3 = Axis(fig[2,1], ylabel=L"XCO_2", ylabelsize=20, yticklabelsize=18, yaxisposition=:right, ytickcolor=XCOclr, yticklabelcolor=XCOclr, ylabelcolor=XCOclr, ygridvisible=false, xgridvisible=false); hidespines!(ax3); hidexdecorations!(ax3)
            dx, dy = 0.1(maximum(FMQ ? logfO2FMQ : logfO2)-minimum(FMQ ? logfO2FMQ : logfO2)), 0.1(maximum(Δ[:,1])-minimum(Δ[:,1]))
            scatterlines!(ax, FMQ ? logfO2FMQ : logfO2, Δ[:,1], color=Fe2clr, strokewidth=0.4, marker=:utriangle, markersize=15, label="XFeO   → $(round(Fe2,digits=4))") # FeO
            scatterlines!(ax, FMQ ? logfO2FMQ : logfO2, Δ[:,2], color=Fe3clr, strokewidth=0.4, marker=:utriangle, markersize=15, label="XFe₂O₃ → $(round(Fe3,digits=4))") # Fe2O3
            xlims!(ax, minimum(FMQ ? logfO2FMQ : logfO2)-0.4dx, maximum(FMQ ? logfO2FMQ : logfO2)+0.4dx)
            xlims!(ax3, minimum(FMQ ? logfO2FMQ : logfO2)-0.4dx, maximum(FMQ ? logfO2FMQ : logfO2)+0.4dx)
            ylims!(ax, 0.0, max(maximum(Δ[:,1]), maximum(Δ[:,2]))+dy); ylims!(ax3, 0.0, 1.0);
            lines!(ax3, FMQ ? logfO2FFMQ : logfO2F, XCO2spec, color=XCOclr, alpha=0.5, linewidth=2.0)
            lines!(ax, [0.0, 0.1], [-0.1, -0.2], color=XCOclr, alpha=0.5, linewidth=2.0, label="Eq. XCO₂ = $(round(XCO2specsol, digits=5))")

        # ==== Solution ====
            lines!(ax, FMQ ? [out.dQFM, out.dQFM] : [out.fO2, out.fO2], [0.0, maximum(Δ[:,1])+dy], linestyle=:dash, color=:black, alpha=0.4)
            scatter!(ax, FMQ ? out.dQFM : out.fO2, Fe2, color=:grey, markersize=15, strokewidth=0.7, marker=:rect, label=("Eq. log₁₀ fO₂ = $(round(FMQ ? out.dQFM : out.fO2, digits=4))" * (FMQ ? " (ΔFMQ)" : " (abs)")))
            scatter!(ax, FMQ ? out.dQFM : out.fO2, Fe3, color=:grey, markersize=15, strokewidth=0.7, marker=:rect)
            scatter!(ax3, FMQ ? out.dQFM : out.fO2, XCO2specsol, color=:green, markersize=15, strokewidth=0.7, marker=:rect)
            scatter!(ax, FMQ ? out.dQFM : out.fO2, Fe3, alpha=0.0, label="(Solid)  ∑Fe³⁺/∑Feᵀ → $RF")
            scatter!(ax, FMQ ? out.dQFM : out.fO2, Fe3, alpha=0.0, label="(Melt)   ∑Fe³⁺/∑Feᵀ → $(round(2Fe3/(2Fe3+Fe2), digits=4))")
            Legend(fig[2,2], ax, framevisible=false, "P = $P GPa | T = $T K")
            display(fig)

    # ==== Carbon =====
        ax2 = Axis3(fig[1,1], xlabel= FMQ ? L"log_{10}\;fO_2\;\mathrm{[\Delta FMQ]}" : L"log_{10}\;fO_2\;\mathrm{[abs]}", ylabel=L"Pressure\;\mathrm{[GPa]}", zlabel=L"XCO_2", azimuth=-π*(3.6/3), elevation=π/12); 
        arrf, arrt, arrp = zeros(ns, ns, nf), zeros(ns, ns, nf), zeros(ns, ns, nf)
        Pr, Prc = LinRange(0.0, 10., ns), LinRange(0.0, 10., nf)
        Tr, Trc = LinRange(c2k(1100.), c2k(1600.), ns), LinRange(c2k(1100.), c2k(1600.), nf)
        logfO2, logfO2c = LinRange(-9, -3, ns), LinRange(-9, -3, nf)
        for f in 1:nf
            for i in 1:ns
                for j in 1:ns
                    XCO2 = eq_XCO2(logfO2[j], Pr[i], Trc[f])
                    arrt[i, j, f] = XCO2; # P-T-fO2surf
                end
            end
        end
        cmap = cpalette(:vik, nf)
        for f in 1:nf
            larr = arrt[:,:,f]'; surface!(ax2, FMQ ? logfO2FMQ : logfO2, Pr, larr, color=fill(cmap[f], ns, ns));
            scatter!(ax2, logfO2[1]-10, Pr[1]-5, -1, color=cmap[f], label="$(round(Trc[f], digits=1))")
        end
        # P-logfO2-Tsurf
        xlims!(ax2, first(FMQ ? logfO2FMQ : logfO2), last(FMQ ? logfO2FMQ : logfO2)); ylims!(ax2, first(Pr), last(Pr)); zlims!(ax2, 0.1, 1.01)
        Legend(fig[1,2], ax2, "Temperature (K)", framevisible=false)
    display(fig)

    return Δ

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
@inline move05(x) = 0.5(x[1:end-1] .+ x[2:end])
function partition_Oₑₓ(P::K, T::K; p::K=0.1, ϕ::K=0.01, Rs::K=0.02, Rf::K=0.0, nr=50, niter=100) where {K <: Real}

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
    Xox = ["SiO2", "MgO", "FeO", "CaO", "Al2O3", "Na2O", "Cr2O3", "Fe2O3"]
    X   = Vector(p.*XB .+ (1-p).*XH)
    
    # Compose composition dictionaries
    Xs = oxidize_bulk(X, Xox, Rs, FeFormat="FeO_O", dict=true, wt_out=true, frac=true); 
    Xm = oxidize_bulk(Vector(XB), Xox, Rf, FeFormat="FeO_O", dict=true, wt_out=true, frac=true);

    # Compute total oxygen budget
    TOₑₓ = Xs["O"]*Ms + Xm["O"]*Mf # kg
    sfOₑₓ = [0.99TOₑₓ, 0.01TOₑₓ]

    # Generate solid fO₂ space
    Rlist = LinRange(0.0001, 0.3, nr)
    Xlist = Vector{Vector{Float64}}(); [push!(Xlist, oxidize_bulk(X, Xox, Rlist[i], wt_out=true, frac=true, FeFormat="FeO_O", onlyvals=true)) for i in 1:nr]
    sOₑₓlist = [Xlist[i][end] for i in 1:nr] # mass fraction Oₑₓ in solid
    # -- Minimizer call
        data    = Initialize_MAGEMin("sb24", verbose=false);
        out = multi_point_minimization(10P*ones(nr), k2c(T)*ones(nr), data, X=Xlist, Xoxides=Xox, name_solvus=true, sys_in="wt")
        Finalize_MAGEMin(data);
    # -- Create interpolation object
        sfO2 = interpolate((sOₑₓlist,), [out[i].fO2 for i in eachindex(out)], Gridded(Linear()))

    # Newton solver
    residual = zeros(5, niter); residual[1,:].*=NaN
    for iter = 1:niter
        sOₑₓ = sfOₑₓ[1]/Ms
        residual[4,iter] = sfO2(sOₑₓ)
        residual[5,iter] = Hirsch(P, T, Xm, sOₑₓ, TOₑₓ, Ms, Mf, T₀, ΔCₚ, a, b, c, y1, y2, y3, y4, y5, y6, y7, y8, y9, IDV)
        residual[1,iter] = residual[5,iter] - residual[4,iter]
        sfOₑₓ[1]-=1e-3TOₑₓ; sfOₑₓ[2]+=1e-3TOₑₓ
        residual[2,iter] = sfOₑₓ[1]; residual[3,iter] = sfOₑₓ[2]
    end
    idx = argmin(filter(x -> !isnan(x), abs.(residual[1,:])))
    println("residual = $(round(abs(residual[1,idx]), digits=3)) | fO₂ solid = $(round(residual[4,idx], digits=3)) | fO₂ melt = $(round(residual[5,idx], digits=3)) [Oxygen budget partitioning → solid=$(round(residual[2,idx]/TOₑₓ, digits=3))TOₑₓ, melt=$(round(residual[3,idx]/TOₑₓ, digits=3))TOₑₓ]")

    plt, ax = plot(residual[1,:])
    lines!(ax, [0, niter], [0, 0], color=:red)
    display(plt)
    
end

function Hirsch(P, T, Xm, sOₑₓ, TOₑₓ, Ms, Mf, T₀, ΔCₚ, a, b, c, y1, y2, y3, y4, y5, y6, y7, y8, y9, IDV)
    # Extract melt fOₑₓ
    mfOₑₓ = (TOₑₓ-sOₑₓ*Ms)/Mf
    # Work on dummy copy while passing the mfOₑₓ
    X = copy(Xm); X["O"] = mfOₑₓ; [X[k]/=mm[k] for k in keys(X)]
    # Retrace mfOₑₓ → XOₑₓ
    n=sum(values(X)); [X[k]/=n for k in keys(X)]; XFe₂O₃ = X["O"] # 1 mol of [Fe₂O₃] every 1 mol of [Oₑₓ]
    if XFe₂O₃>=0.5X["FeO"] # Too much oxygen!! Above hard-limit.
        return NaN
    else
        # Compute mfO2
        return (log10(XFe₂O₃/(X["FeO"]-2XFe₂O₃)) - b - c/T + (ΔCₚ/R/log(10) * (1 - T₀/T - log(T/T₀))) + IDV/(1e-3R)/T/log(10) 
                            - (1/T)*(y1*X["SiO2"] + y2*0.0 + y3*X["MgO"] + y4*X["CaO"] + y5*X["Na2O"] + y6*0.0 + y7*0.0 + y8*X["SiO2"]*X["Al2O3"] + y9*X["SiO2"]*X["MgO"]))/a
    end
end