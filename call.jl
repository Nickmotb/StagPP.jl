using StagPP, GLMakie, MAGEMin_C

P       = 25.0       # GPa
T       = 1900      # K
p       = 0.0       # 0 -> depleted solid || 1 -> enriched solid
ϕ       = 0.01      # Melt fraction k
TOex    = 2e-2      # Excess oxygen as a wt% of tracer mass
TC      = 1e-1      # Carbon content as wt% of tracer mass
H       = 2e-4      # Free H₂O mass in wt%

# data = Initialize_MAGEMin("sb24", verbose=false)
partition_Oₑₓ(P, T, p, ϕ, 1e-2TOex, 1e-2TC, 1e-2H, verbose=true, data=data, maxiter=1000, debug=false)
# Finalize_MAGEMin(data)