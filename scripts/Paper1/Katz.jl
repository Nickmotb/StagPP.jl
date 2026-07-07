using CairoMakie, ColorSchemes, LsqFit, Interpolations

function cpalette(scheme_sym::Symbol, n::Integer)
    cs = getfield(ColorSchemes, scheme_sym)         # e.g. :vik100 → ColorSchemes.vik100
    ts = range(0.0, 1.0; length=n)                  # n evenly spaced points
    return get.(Ref(cs), ts)                        # cs(t) for each t
end

function Katz(p, c)
    # Default
    kp = p<=8.0 ? p : 8.0
    χ1 = 12.; χ2 = 1.; λ = 0.6; k = 43.; γ = 0.75
    x_sat = χ1 * kp^λ + χ2 * kp
    x = min(c / 0.01, x_sat)
    ΔT = k * x^γ
    # g(P) correction - Generalized Hill function
    P_start = 7.0; P_scale = 6.8; n=5.0
    dp = (p - P_start) / P_scale
    ΔT *= p<=P_start ? 1.0 : 0.15 + (1.0 - 0.15) / (1.0 + dp^n)

    return ΔT
end

c = [1e-3, 1e-2, 1e-1, 2e-1, 2.0]    # H₂O content wt%
p = LinRange(0.1, 25.0, 150)    # Pressure GPa
xlims = (0.0, p[end])
ylims = (-50, 950.0)

mat = zeros(length(p), length(c)) # ΔT(GPa,  H₂O wt%)
for i in eachindex(p)
    mat[i, :] = Katz.(p[i], c)
end

# Litasov ΔT brackets (Pyrolite + 2 H₂O wt%) + 22GPa approx ΔT 2300 to 2150 for 50% to 38%, approx 120C ?
    lit = [ 13.5 16.25 19.5 24.5;   # Pressure GPa
            2050  2080 1950 2220;   # Last dry no-melt
            2100  2120 2160 2350;   # First dry melt
            1550  1600 1800 2000;   # Last wet no-melt
            1600  1900 2050 2180]   # First wet melt
    litΔT = zeros(4, size(lit, 2))
    litΔT[1,:] .= lit[1,:]                                                  # Pressure GPa
    litΔT[2,:] .= 0.5(lit[2,:] .+ lit[3,:]) .- 0.5(lit[4,:] .+ lit[5,:])    # ΔT mid-point
    litΔT[3,:] .= lit[3,:] .- lit[4,:]                                      # Δ upper error
    litΔT[4,:] .= lit[2,:] .- lit[5,:]                                      # Δ lower error
    litΔT = cat(litΔT, [22, 100, 150, 50], dims=2)
    litΔT = max.(litΔT, 0.0)

# Iwamori peridotite 0.2 H₂O wt% ΔT bracket
    an_iwa_p = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 
                9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 
                    16.5, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 24.0, 24.5, 25.0]
    an_iwa_t = [1180, 1270, 1351, 1383, 1451, 1533, 1597, 1647, 1687, 1727, 1763, 1796, 1825, 1853, 
                1878, 1900, 1921, 1941, 1957, 1971, 1986, 1996, 2008, 2015, 2025, 2031, 2036, 2042, 2047, 
                    2047, 2065, 2083, 2101, 2118, 2132, 2147, 2162, 2172, 2183, 2193, 2201, 2214, 2217, 2227, 2240, 2258, 2273, 2288, 2303, 2317]
    h_iwa_p = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 
                9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 
                    17.5, 18.0, 18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0]
    h_iwa_t = [1062, 1083, 1102, 1108, 1013, 1049, 1085, 1117, 1147, 1172, 1192, 1208, 1223, 1235, 1242, 
                1249, 1254, 1259, 1264, 1269, 1274, 1277, 1281, 1283, 1288, 1289, 1448, 1633, 1818, 1962, 1980, 
                    1998, 2015, 2038, 2049, 2063, 2078, 2101, 2108, 2122, 2133, 2146, 2157, 2161]
    an_iwa_itp = Interpolations.interpolate((an_iwa_p,), an_iwa_t, Gridded(Linear()))
    h_iwa_itp = Interpolations.interpolate((h_iwa_p,), h_iwa_t, Gridded(Linear()))
    p_iwa = LinRange(14.5, 22.0, 20)
    ΔT_iwa = an_iwa_itp.(p_iwa) .- h_iwa_itp.(p_iwa)

# Decorations
    clr = cpalette(:balance, length(c))
    fig = Figure(size = (1000, 500))
    ax = Axis(fig[1, 1], xgridvisible=false, ygridvisible=false)
    hidedecorations!(ax)
    xlims!(ax, xlims[1], xlims[2])
    ylims!(ax, ylims[1], ylims[2])
    lines!(ax, [8.0, 8.0], [-25, 900.0], color=:gray55, linestyle=:dash)
    lines!(ax, [14, 14], [-25, 750.0], color=:olive, linestyle=:dash)
    text!(ax, 8.3, 840, text="Katz (2003) boundary", color=:gray55, font=:italic)
    text!(ax, 14.3, 680, text="Spinel onset", color=:olive, font=:italic)
    # H2O labels
    text!(ax, 2.0, -45, text=L"1\times 10^{-3}\;wt\%\;H_2O", color=clr[5])
    text!(ax, 2.5, 50, text=L"1\times 10^{-2}\;wt\%\;H_2O", color=clr[4])
    text!(ax, 3.5, 250, text=L"1\times 10^{-1}\;wt\%\;H_2O", color=clr[3])
    text!(ax, 4.5, 415, text=L"0.2\;wt\%\;H_2O", color=clr[2])
    text!(ax, 4.5, 620, text=L"2.0\;wt\%\;H_2O", color=clr[1], rotation=π/5)

# Katz
    for i in eachindex(c)
        lines!(ax, p, mat[:, i], color = clr[length(c)-i+1], label = "$(round(c[i], digits=3)) wt%", linewidth=3.0)
    end
    # Legend(fig[1,2], ax, "H₂O content")
    ax = Axis(fig[1, 1], titlesize=25, xlabel = "Pressure [MPa]", ylabel = "ΔT [K]", xgridvisible=false, ygridvisible=false,
                xlabelsize=25, ylabelsize=25, xticklabelsize=20, yticklabelsize=20)
    xlims!(ax, xlims[1], xlims[2])
    ylims!(ax, ylims[1], ylims[2])
# ΔT brackets
    for i in axes(litΔT, 2)
        # Errors
        scatter!(ax, litΔT[1,i], litΔT[3,i], color=clr[1], marker=:hline, strokewidth=0.0, markersize=15)
        scatter!(ax, litΔT[1,i], litΔT[4,i], color=clr[1], marker=:hline, strokewidth=0.0, markersize=15)
        lines!(ax, [litΔT[1,i], litΔT[1,i]], [litΔT[3,i], litΔT[4,i]], color=clr[1], linewidth=0.5, linestyle=:dash)
        # Points
        scatter!(ax, litΔT[1,i], litΔT[2,i], color=clr[1], marker=:utriangle, strokewidth=0.5, markersize=15, label="Litasov 2002 (Pyrolite + 2 wt% H₂O)")
    end
    scatter!(ax, p_iwa, ΔT_iwa, color=clr[2], marker=:rect, strokewidth=0.5, markersize=15, label="Iwamori 2004 (Peridotite + 0.2 wt% H₂O)")
axislegend(ax, "Study measurements", unique=true, position=:rt, titlefont=:bold_italic)
display(fig)

