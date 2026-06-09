using StagPP, GLMakie

function f()
    nlines = 7
    it = (15, 40)

    # Solver
    a = solve_sH2O_fO2(5, 5, exprt=false, s=false, fO2=false, verbose=false)

    # do stuff
    fig = Figure(size=(1100,550))
    ax = Axis(fig[2,1], 
                xlabelsize=20, ylabelsize=20,
                xlabel=L"Pressure\;\mathrm{[GPa]}", ylabel=L"\Delta \rho - decrease\;\mathrm{[\%]}",
                xgridvisible=false, ygridvisible=false)
    ax2 = Axis(fig[2,2], 
                xlabelsize=20, ylabelsize=20,
                xlabel=L"Pressure\;\mathrm{[GPa]}", ylabel=L"\Delta \rho - decrease\;\mathrm{[\%]}",
                xgridvisible=false, ygridvisible=false)
    ax3 = Axis(fig[3,1], 
                xlabelsize=20, ylabelsize=20,
                xlabel=L"Pressure\;\mathrm{[GPa]}", ylabel=L"\Delta \rho - decrease\;\mathrm{[\%]}",
                xgridvisible=false, ygridvisible=false)
    ax4 = Axis(fig[3,2], 
                xlabelsize=20, ylabelsize=20,
                xlabel=L"Pressure\;\mathrm{[GPa]}", ylabel=L"\Delta \rho - decrease\;\mathrm{[\%]}",
                xgridvisible=false, ygridvisible=false)
    l = length(a.ρH)
    T = "$(Int(floor((a.ρT[it[1]])))) K"
    T2 = "$(Int(floor((a.ρT[it[2]])))) K"
    step = (l-1) / (nlines)
    clr = cpalette(:Blues, l);
    for il in 1:step:l
        i = Int(floor(il))
        v = round(a.ρH[i],digits=3)
        scatterlines!(ax, a.ρP,  a.SΔρ[:,it[1],i,1], label="$v", strokewidth=0.9, color=clr[i], marker=:rect)
        scatterlines!(ax2, a.ρP, a.MΔρ[:,it[1],i,1], strokewidth=0.9, color=clr[i], marker=:rect)
        scatterlines!(ax3, a.ρP, a.SΔρ[:,it[2],i,1], label="$v", strokewidth=0.9, color=clr[i], marker=:rect)
        scatterlines!(ax4, a.ρP, a.MΔρ[:,it[2],i,1], strokewidth=0.9, color=clr[i], marker=:rect)

    end
    text!(ax, (0.7, 0.8), text="Solid ($T)", font=:bold_italic, space=:relative, fontsize=17)
    text!(ax2, (0.58, 0.8), text="Hydrous Melt ($T)", font=:bold_italic, space=:relative, fontsize=17)
    text!(ax3, (0.7, 0.8), text="Solid ($T2)", font=:bold_italic, space=:relative, fontsize=17)
    text!(ax4, (0.58, 0.8), text="Hydrous Melt ($T2)", font=:bold_italic, space=:relative, fontsize=17)
    Legend(fig[1,:], ax, orientation=:horizontal, L"\textbf{H_2O\;Concentration\;(wt\%)}", labelsize=14, titlesize=20)

    idx = argmin(abs.(a.MΔρ[:,it[1],1,1]))
    lines!(ax2, [a.ρP[idx], a.ρP[idx]], [-4., 4.], color=:gray, linestyle=:dash)
    t = round(a.ρP[idx], digits=2)
    text!(ax2, 1.01a.ρP[idx], 4.1, text="$t GPa", color=:gray)
    idx = argmin(abs.(a.MΔρ[:,it[2],1,1]))
    lines!(ax4, [a.ρP[idx], a.ρP[idx]], [-4., 4.], color=:gray, linestyle=:dash)
    t = round(a.ρP[idx], digits=2)
    text!(ax4, 1.01a.ρP[idx], 4.1, text="$t GPa", color=:gray)
    GLMakie.save("./scripts/Paper1/imgs/densities.png", fig)
    display(fig)
end
f()
