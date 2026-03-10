# --- Solid-melt fO2-XCO₂ solution space
"""
    Computes and presents the solution space of the fO2-Oexcess system solved within stagYY to define melt equilibrium.

    \t Basic usage: \t SM_fO2_Oex_solution_space(ns; kwargs...)

    Optional arguments (kwargs):

    • -- General arguments --

        - P :: \t\t\t-->\t Pressure in GPa [default: 3.]
        - T :: \t\t\t-->\t Temperature in K [default: 1000.]
        - slimit :: \t\t\t-->\t Melt XCO₂ solubility threshold [default: 0.5]
        - Oexrange :: \t\t\t-->\t fO2 range explored by the funtion in log units [default: (-10, 1)]
        - Oexrange :: \t\t\t-->\t Molar proportion range of the solid tracer excess oxygen budget [default: (0.0, 0.5)]
        - Ocut :: \t\t\t-->\t Curve cut at fixed Oex where the program solves for. Nothing -> Automatic [default: nothing]

    • -- Plot arguments --

        - fsize :: \t\t\t-->\t Figure size [default: (1200, 800)]
        - savein :: \t\t\t-->\t Path where image is saved. "" to skip [default: ""]

"""
function SM_fO2_Oex_solution_space(; ns=60, fsize=(1200, 800), savein="", fO2range=(-13,1), P=3.5, T=1821., slimit=0.8, Oexrange=(1e-7,1e-3), Ocut=nothing, crblim=0.4)

    smass, fmass = 3.342e16, 2.120e13
    mmO = 16.00

    # IDV = ∫ΔVdP
    IDV = solve_∫ΔVdP([P-0.05P, P, P+0.05P],[T-0.05T, T, T+0.05T])[2,2,1]

    # Generate fO2 and XCO₂ arrays
    logfO2 = LinRange(fO2range..., ns)
    Oex = LinRange(Oexrange..., ns)

    # Composition dictionary
    # X = Dict( "SiO2" => 0.4343, "MgO" => 0.4593, "FeO" => 0.0834, "CaO" => 0.0090, "Al2O3" => 0.0100, "Na2O" => 0.0001, "Cr2O3" => 0.0030, "Fe2O3" => 0.0)
    X = Dict( "SiO2" => 0.5042, "MgO" => 0.0977, "FeO" => 0.0710, "CaO" => 0.1254, "Al2O3" => 0.1680, "Na2O" => 0.0223, "Cr2O3" => 0.0007, "Fe2O3" => 0.0)
    mm = Dict("SiO2" => 60.08, "Al2O3" => 101.96, "CaO" => 56.08, "MgO" => 40.30, "FeO" => 71.85, "Fe2O3" => 159.69, "K2O" => 94.2, "Na2O" => 61.98, 
                "TiO2" => 79.88, "O" => 16.0, "Cr2O3" => 151.99, "MnO" => 70.937, "H2O" => 18.015, "CO2" => 44.01, "S" => 32.06, "P2O5" => 141.9445, "Fe" => 55.845)

    # BaseFMQ
    data = Initialize_MAGEMin(P<=7. ? "um" : "sb24", verbose=false);
    out = single_point_minimization(10P, T-273.15, data, X=collect(values(X)), Xoxides=collect(keys(X)));
    Finalize_MAGEMin(data); ΔFMQ = out.fO2 - out.dQFM

    # Sun and Yao 2024 (1 Oex per 1 Fe₂O₃)
    a0=-0.0001; a1=0.0002; a2=-0.0003; a3=-0.0004; a4=-0.0005; a5=-0.0006; a6=-0.0007; a7=-0.0008; a8=-0.0009; a9=-0.0010; a10=-0.0011; a11=-0.0012; a12=-0.0013; h=2.1410
    function reverseSY(logfO2)
        Ω = a1 + a2*T^(1.5) + a3*log(T)
        Femax = 0.5X["FeO"]; Femin = 0.0; Feguess = 0.5(Femax+Femin)
        for n in 1:10
            XFe2O3 = Feguess
            ψ = a4*log(X["FeO"] - 2XFe2O3) + a5*((X["FeO"] - 2XFe2O3)^0.5) + a6*X["SiO2"]^3 + a7*X["Al2O3"] + a8*0.0 + a9*X["CaO"] + a10*X["MgO"] + (a11 + a12*(X["FeO"] - 2XFe2O3))*(X["Na2O"] + 0.0)
            guess = (a0*((X["FeO"] - 2XFe2O3)^0.5)*log10(XFe2O3/(X["FeO"] - 2XFe2O3)) + Ω + ψ + h*IDV + 4*log10(XFe2O3/(X["FeO"] - 2XFe2O3))) - ΔFMQ
            (abs(guess) < 1e-6) && break
            (guess > logfO2) ? (Femax = Feguess) : (Femin = Feguess)
            Feguess = 0.5(Femax + Femin)
        end
        # Construct copy
        Xc = copy(X); Xc["FeO"] -= 2Feguess; Xc["Fe2O3"] = Feguess
        Xc = Dict(k => v/sum(values(Xc)) for (k,v) in Xc)
        # Convert to mass fraction
        Xc = Dict(k => v*mm[k] for (k,v) in Xc)
        Xc = Dict(k => v/sum(values(Xc)) for (k,v) in Xc)
        return Xc["Fe2O3"]
    end

    # Stagno and Frost XCO2 equilibrium (2 Oex per 1 CO₂)
    function eq_XCO2(logfO2)
        XCO2 = min(10^(logfO2 - 5.44 + 21380/T - 0.078(1e4P-1)/T + ΔFMQ), min(slimit, crblim))
        XCO2 = 2(XCO2/(sum(values(X))+XCO2)) # Normalized mols of O
        Xc = copy(X)
        Xc = Dict(k => v*mm[k] for (k,v) in Xc)
        XCO2 = (XCO2*mmO)/(sum(values(Xc))+XCO2*mmO) # Normalized mass f of O
    end

    # Δ function
    function ΔR(logfO2, Oex)
        R = zeros(length(logfO2),length(Oex))
        for i in 1:ns
            for j in 1:ns
                R[i,j] = fmass*(reverseSY(logfO2[i]) + eq_XCO2(logfO2[i])) - smass*Oex[j]
            end
        end
        return R
    end
    res = ΔR(logfO2, Oex)

    # Plot
    if isnothing(Ocut)
        ncut = Int(floor(ns/2))
        Ocut= Oex[ncut]
        crv = res[:,ncut]
    else
        itpmap = interpolate((logfO2, Oex), res, Gridded(Linear()))
        crv = itpmap(logfO2, Ocut)
    end
    fig = Figure(size=fsize)
    ax = Axis3(fig[1:5, 1], xlabel=L"Equilibrium\;fO_2", ylabel=L"O_{excess}\;budget\;[prop.\;tracer\;mass]", zlabel=L"Residual");
        surface!(ax, logfO2, Oex, res, colormap=:vik100)
        scatter!(ax, logfO2, Ocut*ones(ns), crv, color=:yellow, markersize=12, strokewidth=0.5, strokecolor=:black)

    ax = Axis(fig[2:4, 2], xlabel=L"Equilibrium\;fO_2");
        scatterlines!(ax, logfO2, crv, color=:red, strokewidth=0.4, strokecolor=:black, label="Residual ($P GPa | $T K)")
        dy, dx = 0.1(crv[end] - crv[1]), 0.1(logfO2[end] - logfO2[1])
        ylims!(ax, minimum(crv)-dy, maximum(crv)+dy)

        # resulting fO2
        itp = interpolate((logfO2,), crv, Gridded(Linear()))
        logfO2fine = LinRange(logfO2[1], logfO2[end], 10ns)
        idx = argmin(abs.(itp(logfO2fine)))
        lines!(ax, [logfO2fine[idx], logfO2fine[idx]], [minimum(crv), maximum(crv)], color=:black, linestyle=:dash)
        text!(ax, logfO2fine[idx]-0.1, maximum(crv)-2.6dy, text="log fO₂ = $(round(logfO2fine[idx], digits=3))", color=:black, rotation=0.5π)

        # XCO2
        scatter!(ax, logfO2, eq_XCO2.(logfO2), color=:green, strokewidth=0.4, strokecolor=:black, marker=:utriangle, label="XCO₂ = $(round(eq_XCO2(logfO2fine[idx]), digits=3))")
        lines!(ax, [logfO2[1], logfO2[end]], [slimit, slimit], color=:green, linestyle=:dash)
        text!(ax, logfO2[1]+0.5dx, slimit+0.2dy, text="XCO₂ saturation = $slimit", color=:green)
        scatter!(ax, logfO2[1], crv[1], alpha=0.0, label="Solid O budget = $(round(log10(smass*Ocut), digits=3))")
        rOex = fmass*(reverseSY(logfO2fine[idx]) + eq_XCO2(logfO2fine[idx]))
        scatter!(ax, logfO2[1], crv[1], alpha=0.0, label="Requested O budget = $(round(log10(rOex), digits=3))")
        axislegend(ax, position=:lb)

    display(fig)
    !isempty(savein) && save(joinpath(savein, "SM_fO2_XCO2_solution_space.png"), fig)

end