using CairoMakie, ColorSchemes

function cpalette(scheme_sym::Symbol, n::Integer)
    cs = getfield(ColorSchemes, scheme_sym)         # e.g. :vik100 → ColorSchemes.vik100
    ts = range(0.0, 1.0; length=n)                  # n evenly spaced points
    return get.(Ref(cs), ts)                        # cs(t) for each t
end

function Katz(p, c)
    χ1 = 12.; χ2 = 1.; λ = 0.6; k = 43.; γ = 0.75
    x_sat = χ1 * p^λ + χ2 * p
    x = min(c / 0.01, x_sat)
    return k * x^γ
end

c = LinRange(1e-4, 1.0, 20)
p = LinRange(0.1, 10.0, 150)

mat = zeros(length(p), length(c))
for i in eachindex(p)
    mat[i, :] = Katz.(p[i], c)
end

clr = cpalette(:Blues, length(c))
fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1], xlabel = "Pressure (MPa)", ylabel = "ΔT (K)", title = "Katz vs Pressure")
for i in eachindex(c)
    lines!(ax, p, mat[:, i], color = clr[i], label = "Katz(p, 0.01)", linewidth=5.0)
end
display(fig)

