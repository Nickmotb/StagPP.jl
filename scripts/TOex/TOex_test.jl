using StagPP, GLMakie, StaticArrays, MAGEMin_C, ColorSchemes

# Params
PT      = SA[km2GPa(10.) 1300.
           km2GPa(150.) 1400.
           km2GPa(410.) 1650.
           km2GPa(660.) 1900.]
p       = 0.0
ϕ       = SA[1e-2, 1e-1]        # Melt fraction k
TOex    = SA[1e-1, 1.0]        # Excess oxygen as a wt% of tracer mas
TC      = SA[1e-3, 1.0]         # Carbon content as wt% of tracer mass
H       = SA[1e-2, 1e-1]         # Free H₂O mass in wt%

# Initialise MAGEMin domain
data = Initialize_MAGEMin("sb24", verbose=false)

# Initialise output data array
# sz3 -> convergence (1) + eqfO2/fmq (3) + sFe/sCarb/mFe/mCO2/COH portions (8) + COH speciation (13) + melt XCO2 (14) + carbon fractions (18) + solid iron speciation (23) + dIW (24)
sz1 = size(PT, 1); sz2 = 4; sz3 = 24; sz4 = length(TOex)
mat = zeros(sz1, sz2, sz3, sz4)

# i -> iterator over PT conditions
# j -> iterator over singular increments (1: default | 2: more melt | 3: more carbon | 4: more water)
# k -> iterator over TOex

for i in axes(mat, 1)
    for j in axes(mat, 2)
        for k in axes(mat, 4)
            out = partition_Oₑₓ(PT[i,1], PT[i,2], p, ϕ[j==2 ? 2 : 1], 1e-2TOex[k], 1e-2TC[j==3 ? 2 : 1], 1e-2H[j==4 ? 2 : 1], 
                                    verbose=false, data=data, out=true)
            mat[i,j,1,k] = out.converged ? 1 : 0    # Converged ?
            mat[i,j,2,k] = out.fO2                  # Equilibrium fO₂
            mat[i,j,3,k] = out.FMQ                  # fO₂ relative to the FMQ buffer
            mat[i,j,4,k] = out.sOₑₓ                 # fraction of TOex going to solid iron
            mat[i,j,5,k] = out.carbO                # fraction of TOex precipitating as solid carbonate
            mat[i,j,6,k] = out.mOₑₓ                 # fraction of TOex going to melt iron
            mat[i,j,7,k] = out.carbonO              # fraction of TOex going to melt as CO2
            mat[i,j,8,k] = out.cohsink              # fraction of TOex going to the COH fluid
            mat[i,j,9,k] = out.XH2O                 # Molar fraction of H₂O in the COH fluid
            mat[i,j,10,k] = out.XCO2                # Molar fraction of CO₂ in the COH fluid
            mat[i,j,11,k] = out.XCH4                # Molar fraction of CH₄ in the COH fluid
            mat[i,j,12,k] = out.XCO                 # Molar fraction of CO in the COH fluid
            mat[i,j,13,k] = out.XH2                 # Molar fraction of H₂ in the COH fluid
            mat[i,j,14,k] = out.XCO₂                # Molar fraction of CO₂ in the melt
            mat[i,j,15,k] = out.fC_gph              # fracton of carbon precipitating as graphite/diamond
            mat[i,j,16,k] = out.fC_carb             # fracton of carbon precipitating as carbonate
            mat[i,j,17,k] = out.fC_melt             # fracton of carbon dissolving in the melt
            mat[i,j,18,k] = out.fC_fluid            # fracton of carbon dissolving in the fluid
            # Fe³⁺ and Fe⁰ can now BOTH be non-zero: metal comes from MAGEMin's assemblage, not
            # from sOₑₓ, so disproportionated iron is counted even though it moves no oxygen.
            mat[i,j,19,k] = out.FeR_solid           # solid Fe³⁺/ΣFe
            mat[i,j,20,k] = out.FeM_solid           # solid Fe⁰/ΣFe, total
            mat[i,j,21,k] = out.FeM_diss            # Fe⁰/ΣFe made by dissociation (FeO → Fe⁰ + ½O₂)
            mat[i,j,22,k] = out.FeM_disp            # Fe⁰/ΣFe made by disproportionation (3Fe²⁺ → 2Fe³⁺ + Fe⁰)
            mat[i,j,23,k] = out.nFe_solid           # moles Fe in the solid — 0 marks an ironless cell
            mat[i,j,24,k] = out.IW                  # fO₂ relative to the IW buffer

        end
    end
end
Finalize_MAGEMin(data)
