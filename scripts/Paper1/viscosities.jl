using StagPP, LaTeXStrings, GLMakie

etafac(cH2O, ∂η, r, c₀) = @.  min(max(1. / (∂η + (1-∂η)*(cH2O/c₀)^r), 10^(-1.5)), 1.0)

function f()
    c₀  = 0.01
    ∂η  = 0.01 # low value -> stronger reduction
    r   = 0.75 #0.5(0.35+0.63)
    H2O = LinRange(5e-3, 3, 10000) # wt%
    η   = etafac(H2O, ∂η, r, c₀)
    xlm = (first(H2O), last(H2O))

    # -- Plot --
        # do stuff
        fig = Figure(size=(1000,900))
        ax = Axis(fig[1,1], yscale=log10, xscale=log10,
                    xlabelsize=35, ylabelsize=35, xticklabelsize=20, yticklabelsize=20,
                    xlabel=L"c^{H_2O}\;\mathrm{[wt\%]}", ylabel=L"\eta_\text{factor}",
                    xgridvisible=false, ygridvisible=false)
        # Lims
            ylm = (1e-2, 1e1)
            ylims!(ax, ylm[1], ylm[2]); xlims!(ax, xlm[1], xlm[2])
        # Identity cross
            yc₀ = etafac(c₀, ∂η, r, c₀)
            c₀_ppm = 1e4c₀
            lines!(ax, [c₀, c₀], [yc₀-0.4, yc₀+0.7], color=:red, linestyle=:dash)
            lines!(ax, [c₀, c₀], [ylm[1], yc₀-0.4], color=:red, linestyle=:dash, alpha=0.15)
            lines!(ax, [c₀, c₀], [yc₀+0.7, 10^(0.5)], color=:red, linestyle=:dash, alpha=0.15)
            text!(ax, 0.6c₀, 10^(0.55), text=L"c_0^{H_2O}\;chosen\;at\;%$c₀\;\mathrm{wt\%}\;(%$c₀_ppm\;ppm)", color=:red, fontsize=16)
        # Lower limit (Druzhbin et al., 2022)
            text!(ax, xlm[1] + 0.5*xlm[1] , 10^(-1.44), text=L"Druzhbin\;et\;al.\;(2022) - \sim\;1.5\;orders\;of\;magnitude", color=:gray, fontsize=18)
            lines!(ax, [xlm[1], xlm[2]], [10^(-1.5), 10^(-1.5)], color=:gray, alpha=0.4)
        # Ringwoodite (Fei & Katsura, 2020)
            lines!(ax, [0.8, 0.8], [1e-2, 2e1], color=:skyblue2, linestyle=:dash)
            lines!(ax, [1.2, 1.2], [1e-2, 2e1], color=:skyblue2, linestyle=:dash)
            lines!(ax, [1.0, 1.0], [1e-2, 2e1], color=:skyblue2, linestyle=:dash)
            band!(ax, [0.8, 1.2], [1e-2, 1e-2], [2e1, 2e1], color=:skyblue2, alpha=0.25)
            text!(ax, 1.2 + 0.35, 10^(-0.7), rotation=π/2, text=L"Ringwoodite\;H_2O\;saturation", fontsize=18, color=:royalblue1)
            text!(ax, 1.2 + 0.65, 10^(-0.6), rotation=π/2, text=L"(Fei\;&\;Katsura,\;2022)", fontsize=18, color=:royalblue1)
        # Inset
            ax2 = Axis(fig, bbox=BBox(460, 770, 600, 790), xticklabelsize=15, yticklabelsize=15, xgridvisible=false, ygridvisible=false,
                        xlabel=L"c^{H_2O}\;\mathrm{[ppm]}", ylabel=L"\eta_\text{factor}",
                        xlabelsize=20, ylabelsize=20)
            idx = findfirst(H2O.>=1.0)
            scatter!(ax2, 1e4H2O[1:idx], η[1:idx], color=:black, markersize=2)
        # Plot
        scatter!(ax, H2O, η, color=:black)
        scatter!(ax, 1.0, 10^(-1.5), color=:royalblue1, markersize=17)
        scatter!(ax, c₀, yc₀, color=:red, markersize=17)
    display(fig)
    # Save
    GLMakie.save("./scripts/Paper1/imgs/viscosities.png", fig)
end
f()