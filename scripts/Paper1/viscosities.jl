using StagPP

etafac(cH2O, ∂η, r, c₀) = @.  max(1. / (∂η + (1-∂η)*(cH2O/c₀))^r, 10^(-1.5))

function f()
    c₀  = 0.01
    ∂η  = 0.3 # low value -> stronger reduction
    r   = 0.81
    H2O = LinRange(1e-7, 10.0, 10000)
    η   = etafac(H2O, ∂η, r, c₀)

    # do stuff
    fig = Figure(size=(900,550))
    ax = Axis(fig[1,1], yscale=log10, xscale=log10,
                xlabelsize=30, ylabelsize=30, xticklabelsize=16, yticklabelsize=16,
                xlabel=L"c^{H_2O} / c_0^{H_2O}", ylabel=L"\eta_\text{factor}",
                xgridvisible=false, ygridvisible=false)
    # Lims
    ylims!(ax, 1e-2, 2e1); xlims!(ax, 1e-3, 1e3)
    # Identity cross
    scatter!(ax, 1.0, 1.0, color=:red, markersize=20)
    lines!(ax, [1.0, 1.0], [0.6, 1.7], color=:red, linestyle=:dash)
    lines!(ax, [1.0, 1.0], [1e-2, 0.6], color=:red, linestyle=:dash, alpha=0.15)
    lines!(ax, [0.6, 1.7], [1.0, 1.0], color=:red, linestyle=:dash)
    lines!(ax, [H2O[1], 0.6], [1.0, 1.0], color=:red, linestyle=:dash, alpha=0.15)
    lines!(ax, [1.0, 1.0], [1.7, 10^(0.5)], color=:red, linestyle=:dash, alpha=0.15)
    text!(ax, 10^(-0.3), 10^(0.55), text=L"c_0^{H_2O}\;chosen\;at\;%$c₀", color=:red, fontsize=16)
    # Lower limit (Druzhbin et al., 2022)
    text!(ax, 10^(-2.7) , 10^(-1.44), text=L"Druzhbin\;et\;al.\;(2022) - \sim\;1.5\;orders\;of\;magnitude", color=:gray, fontsize=18)
    lines!(ax, [1e-3, 1e3], [10^(-1.5), 10^(-1.5)], color=:gray, alpha=0.4)
    # Ringwoodite (Fei & Katsura, 2020)
    lines!(ax, [0.8/c₀, 0.8/c₀], [1e-2, 2e1], color=:skyblue2, linestyle=:dash)
    lines!(ax, [1.2/c₀, 1.2/c₀], [1e-2, 2e1], color=:skyblue2, linestyle=:dash)
    band!(ax, [0.8/c₀, 1.2/c₀], [1e-2, 1e-2], [2e1, 2e1], color=:skyblue2, alpha=0.25)
    text!(ax, 10^(2.25), 10^(-0.7), rotation=π/2, text=L"Ringwoodite\;H_2O\;saturation", fontsize=16, color=:royalblue1)
    text!(ax, 10^(2.42), 10^(-0.6), rotation=π/2, text=L"(Fei\;&\;Katsura,\;2022)", fontsize=16, color=:royalblue1)
    # Plot
    scatterlines!(ax, H2O./c₀, η, color=:black)
    display(fig)
    # Save
    GLMakie.save("./scripts/Paper1/imgs/viscosities.png", fig)
end
f()