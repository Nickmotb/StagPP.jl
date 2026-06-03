using StagPP, GLMakie

# P = 2.0
# T = 1300.
# p = 0.0
# ϕ = 0.02
# TOex = 0.10
# TC = 0.05
# partition_Oₑₓ(P, T, p, ϕ, TOex, TC, plotevo=true)

a = solve_sH2O_fO2(150, 150, fO2=true, s=true)
plot_sf(a)