# Stagno and Frost XCO2 equilibrium (2 Oex per 1 CO₂) (Melt C⁰ -> C⁴⁺ ==> fO₂ to XCO₂ mapping)
function eq_XCO2(logfO2, P, T, Xc)
    XCO2min=0.0; XCO2max=1.0; XCO2guess=0.5(XCO2min+XCO2max); residual=0.0
    for n in 1:10
        # Initial guess
        XCO2 = XCO2guess
        # Model
        guess = 5.44 - 21380/T + 0.078(1e4P-1)/T + log10(XCO2)
        # Convergence
        residual = abs(guess - logfO2)
        (residual < 1e-6) && break
        (guess > logfO2) ? (XCO2max = XCO2guess) : (XCO2min = XCO2guess)
        XCO2guess = 0.5(XCO2max + XCO2min)
    end
    # XCO2
    return XCO2guess, residual
end

Hirsc(P, T, X, T₀, ΔCₚ, a, b, c, y1, y2, y3, y4, y5, y6, y7, y8, y9, IDV) = (log10(X["Fe2O3"]/(X["FeO"])) - b - c/T + (ΔCₚ/R/log(10) * (1 - T₀/T - log(T/T₀))) + IDV/(1e-3R)/T/log(10) 
                        - (1/T)*(y1*X["SiO2"] + y2*0.0 + y3*X["MgO"] + y4*X["CaO"] + y5*X["Na2O"] + y6*0.0 + y7*0.0 + y8*X["SiO2"]*X["Al2O3"] + y9*X["SiO2"]*X["MgO"]))/a
function reverseHirsch!(logfO2, P, T, Xc, IDV)
    a=0.1917; b=-1.961; c=4158.1; ΔCₚ=33.25; T₀=1673.15
    y1=-520.46; y2=-185.37; y3=494.39; y4=1838.34; y5=2888.48; y6=3473.68; y7=-4473.6; y8=-1245.09; y9=-1156.86
    if Xc["FeO"] > 0.0
        Femax = 0.5Xc["FeO"]; Femin = 1e-9; Feguess = 0.5(Femax+Femin); residual=0.0
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
function ΔR(logfO2, P, T, X, IDV, slimit, crblim)
    ns = length(logfO2)
    if ns > 1
        Δ = zeros(ns, 3)
        for i in 1:ns
            # Copy composition
            Xc = copy(X)
            # Residual
            # eq_XCO2!(logfO2[i], P, T, Xc, slimit, crblim) # Moves XFeO -> XFe₂O₃ to accomodate the XCO2 resulting from equilibrium at selected fO₂
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

function Hirschmann_fO2_to_R(; RF=0.02, fO2range=(-15., 0.), P=3.5, T=1400., slimit=0.5, crblim=0.8, FMQ=false, ns=60, fsize=(1400, 800), savein="", p=0.2)

    # IDV = ∫ΔVdP
    IDV = solve_∫ΔVdP([P-0.05P, P, P+0.05P],[T-0.05T, T, T+0.05T])[2,2,1]

    # Generate fO2 and XCO₂ arrays
    logfO2 = LinRange(fO2range..., ns)

    # Composition dictionary
    XH = oxidize_bulk(Dict( "SiO2" => 0.4343, "MgO" => 0.4593, "FeO" => 0.0834, "CaO" => 0.0090, "Al2O3" => 0.0100, "Na2O" => 0.0001, "Cr2O3" => 0.0030, "Fe2O3" => 0.0), RF, frac=true, FeO_Fe2O3=true) # Harzburgite
    XB = oxidize_bulk(Dict( "SiO2" => 0.5042, "MgO" => 0.0977, "FeO" => 0.0710, "CaO" => 0.1254, "Al2O3" => 0.1680, "Na2O" => 0.0223, "Cr2O3" => 0.0007, "Fe2O3" => 0.0), 0.0, frac=true, FeO_Fe2O3=true) # MORB

    # Compose parent tracer bulk composition
    Xs = Dict{String, Float64}(); [Xs[k] = p*XB[k] + (1-p)*XH[k] for k in keys(XH)]
    Xm = copy(XB)

    # Parent fO₂ + FMQ buffer
    data = Initialize_MAGEMin(P<=7. ? "um" : "sb24", verbose=false);
    out = single_point_minimization(10P, T-273.15, data, X=collect(values(Xs)), Xoxides=collect(keys(Xs)));
    Finalize_MAGEMin(data); ΔFMQ = out.fO2 - out.dQFM
    logfO2FMQ = logfO2 .- ΔFMQ

    # Compute equilibrated XFeO | XFe₂O₃ | Residual
    Δ = ΔR(logfO2, P, T, Xm, IDV, slimit, crblim)

    fig = Figure(size=fsize)
    Fe2clr = :firebrick2
    Fe3clr = :dodgerblue2
        # ==== Solution space ====
        Fe2, Fe3, res = ΔR(out.fO2, P, T, Xm, IDV, slimit, crblim)
        ax = Axis(fig[1,1], xlabel=FMQ ? L"fO_2\;(\Delta FMQ)" : L"fO_2", ylabel=L"X\;\mathrm{[mol\;fr]}", xlabelsize=20, ylabelsize=20, yticklabelsize=18, xticklabelsize=18, ygridvisible=false, xgridvisible=false)
        # ax2 = Axis(fig[1,1], ylabel=L"XFe_2O_3", ylabelsize=20, yticklabelsize=18, yaxisposition=:right, ytickcolor=Fe3clr, yticklabelcolor=Fe3clr, ylabelcolor=Fe3clr, ygridvisible=false, xgridvisible=false); hidespines!(ax2); hidexdecorations!(ax2)
        dx, dy = 0.1(maximum(FMQ ? logfO2FMQ : logfO2)-minimum(FMQ ? logfO2FMQ : logfO2)), 0.1(maximum(Δ[:,1])-minimum(Δ[:,1]))
        scatterlines!(ax, FMQ ? logfO2FMQ : logfO2, Δ[:,1], color=Fe2clr, strokewidth=0.4, marker=:utriangle, markersize=15, label="XFeO   → $(round(Fe2,digits=4))") # FeO
        scatterlines!(ax, FMQ ? logfO2FMQ : logfO2, Δ[:,2], color=Fe3clr, strokewidth=0.4, marker=:utriangle, markersize=15, label="XFe₂O₃ → $(round(Fe3,digits=4))") # Fe2O3
        xlims!(ax, minimum(FMQ ? logfO2FMQ : logfO2)-0.4dx, maximum(FMQ ? logfO2FMQ : logfO2)+0.4dx)
        ylims!(ax, 0.0, maximum(Δ[:,1])+dy);

        # ==== Solution ====
        lines!(ax, FMQ ? [out.dQFM, out.dQFM] : [out.fO2, out.fO2], [min(Fe3, Fe2)-2dy, max(Fe3,Fe2)+2dy], linestyle=:dash, color=:black, alpha=0.4)
        text!(ax, (FMQ ? out.dQFM : out.fO2)+0.3dx, Fe2+0.1dy, text="Eq. fO₂ = $(round(FMQ ? out.dQFM : out.fO2, digits=8))", rotation=π/4, color=:grey, fontsize=16, font=:bold)
        text!(ax, (FMQ ? out.dQFM : out.fO2)+0.3dx+0.16dx, Fe2+0.1dy-0.28dy, text="XFeO = $(round(Fe2, digits=4))", rotation=π/4, color=Fe2clr, fontsize=16, font=:bold)
        text!(ax, (FMQ ? out.dQFM : out.fO2)+0.3dx+0.32dx, Fe2+0.1dy-0.56dy, text="XFe₂O₃ = $(round(Fe3, digits=4))", rotation=π/4, color=Fe3clr, fontsize=16, font=:bold)
        scatter!(ax, FMQ ? out.dQFM : out.fO2, Fe2, color=:grey, markersize=15, strokewidth=0.7, marker=:rect, label="Eq. fO₂ = $(round(FMQ ? out.dQFM : out.fO2, digits=8))")
        scatter!(ax, FMQ ? out.dQFM : out.fO2, Fe3, color=:grey, markersize=15, strokewidth=0.7, marker=:rect)
        scatter!(ax, FMQ ? out.dQFM : out.fO2, Fe3, alpha=0.0, label="(Solid)  ∑Fe³⁺/∑Feᵀ → $RF")
        scatter!(ax, FMQ ? out.dQFM : out.fO2, Fe3, alpha=0.0, label="(Melt)   ∑Fe³⁺/∑Feᵀ → $(round(2Fe3/(2Fe3+Fe2), digits=4))")
        axislegend(position=:rt, labelsize=20, labelfont=:bold, labelhalign=:center)
        # display(fig)

    fig = Figure(size=fsize)
    ax = Axis(fig[1,1]); ax2 = Axis(fig[1,2])
    arr = zeros(ns, 2)
    logfO2 = LinRange(-4.0, -0.7, ns) .+ ΔFMQ
    for i in 1:ns
        Xc = copy(Xm)
        XCO2, res = eq_XCO2(logfO2[i], P, T, Xc)
        arr[i,1] = XCO2; arr[i,2] = res
    end
    plot!(ax, logfO2, arr[:,1], color=:red)
    plot!(ax2, logfO2, arr[:,2], color=:blue)
    display(fig)

    return Δ

end

