using GLMakie

n1, n2 = 100, 1000
# Dstress
# R = LinRange(0.1, 0.9, n1)
C_wet = 5.0
c0eta_water = 0.01
R = [0.90, 0.85, 0.80]
dstressY_eta_water = @. log(1/R) / (C_wet - c0eta_water)

# yield stress
# ys = 30e6 # MPa
# H2O = LinRange(1e-3, 4.0, n2)
# mat = zeros(n1, n2)
# for i in 1:n1
#     for j in 1:n2
#         mat[i,j] = ys*exp(-dstressY_eta_water[i] * (H2O[j].-c0eta_water))
#     end
# end

# fig = Figure(size=(1000, 700))
# ax = Axis(fig[1,1], xlabel="cH2O", ylabel="Reduction factor")
# hm = heatmap!(ax, H2O, R, mat', colormap=:vik100)
# Colorbar(fig[1,2], hm)
# display(fig)