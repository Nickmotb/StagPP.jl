using StagPP, GLMakie, MAGEMin_C

P = 4.0
T = 1500.
p = 0.0
ϕ = 0.05
TOex = 0.005
TC = 0.05

# data = Initialize_MAGEMin("sb24", verbose=false)
# partition_Oₑₓ(P, T, p, ϕ, TOex, TC, plotevo=true, debugging=true, data=data)

a = solve_sH2O_fO2(150, 150)