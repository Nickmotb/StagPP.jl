using GLMakie

ns = 1000
k0 = 1e-16
totmelt = 0.0

water_percolation_threshold = 0.0
water_relperm_exponent = 3.0
water_perm_refporosity = 1e-3
exp_cutoff = 0.01

A = 1e5
dt = 1e13
dynvisc = 3.5e-5
gP = 8500
w_rho = 1e3
phi_w = LinRange(1e-5, 1.0, 2000)
n = 3 - (min(totmelt, exp_cutoff))/exp_cutoff #  Kozeny-Carman exponent, depends on melt fraction. 2->tube-like percolation | 3->film-like percolation (see Rudge 2018)


phi_weff = @. max(phi_w - water_percolation_threshold, 0.0) / max(1.0 - water_percolation_threshold, 1e-2)
phi_t = @. totmelt + phi_weff

S_w  = @. phi_weff/phi_t
k_rw = @. S_w^water_relperm_exponent
k    = @. k0 * (phi_t/water_perm_refporosity)^n * k_rw

allowed_flow = @.  max((k * A * gP * w_rho * dt) / dynvisc, 0.0)
darcyV = k * gP / dynvisc

fig = Figure(size=(800,800))
ax = Axis(fig[1,1], xlabel="water volume fraction", ylabel="Darcy velocity [cm/yr]", yscale=log10, xscale=log10)
ax2 = Axis(fig[2,1], xlabel="water volume fraction", ylabel="Effective permeability", yscale=log10, xscale=log10)
replace!(darcyV, 0.0 => 1e-15)
replace!(k, 0.0 => 1e-20)
darcy9 = 1e9.*darcyV
plot!(ax, phi_t, darcy9)
plot!(ax2, phi_t, k)

# additions
idx = argmin(abs.(phi_t.-water_perm_refporosity))
text!(ax, 4.0phi_t[1], 1.8darcy9[idx], text="$(round(log10(darcy9[idx]), digits=2)) cm/yr", color=:red, font=:italic)
text!(ax, 4.0phi_t[1], 1e-1darcy9[idx], text="ref = $(round(water_perm_refporosity, digits=6))", color=:red, font=:italic)
lines!(ax, [phi_t[1], phi_t[idx]], [darcy9[idx], darcy9[idx]], color=:red, linestyle=:dash)
lines!(ax, [water_perm_refporosity, water_perm_refporosity], [minimum(darcy9), maximum(darcy9)], color=:red, linestyle=:dash)
text!(ax2, 4.0phi_t[1], 1.4k[idx], text="$(round(log10(k[idx]), digits=2)) m⁻²", color=:red, font=:italic)
text!(ax2, 4.0phi_t[1], 1e-1k[idx], text="ref = $(round(water_perm_refporosity, digits=6))", color=:red, font=:italic)
lines!(ax2, [phi_t[1], phi_t[idx]], [k[idx], k[idx]], color=:red, linestyle=:dash)
lines!(ax2, [water_perm_refporosity, water_perm_refporosity], [minimum(k), maximum(k)], color=:red, linestyle=:dash)

display(fig)


