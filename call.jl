using StagPP, GLMakie, MAGEMin_C

# solve_sH2O_fO2(100, 100;
#                s          = false,          # solubility maps unchanged - do not regenerate them
#                fO2        = true,
#                densities  = false,
#                nR         = 30,
#                Rrange     = (-1.96, 0.4), # MAGEMin ceiling, see above
#                Rlog       = true,           # log-space the Fe³⁺/ΣFe axis
#                verbose    = true,
#                exprt      = true)

P       = 8.0
T       = 1400.
p       = 0.0
ϕ       = 0.95        # Melt fraction k
TOex    = 1e-4        # Excess oxygen as a wt% of tracer mas
TC      = 1e-1         # Carbon content as wt% of tracer mass
H       = 0.001         # Free H₂O mass in wt%

# data = Initialize_MAGEMin("sb24", verbose=false)
partition_Oₑₓ(P, T, p, ϕ, TOex, TC, H, data=data)