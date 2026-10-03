# Lunar H2O storage-capacity table StagH2O_moon.dat (StagH2O layout, no water axis), 0-5 GPa: the NAM laws of sfO2_fns.jl
#  on the melt-free assemblage of the lunar mantle-basalt series (StagPP.jl LunarWater.jl)
using MAGEMin_C, Printf, Serialization
include(joinpath(@__DIR__, "sfO2_fns.jl"))

# ================================ INPUTS (production values) ================================
# Run control: the table overwrites Tables/ beside the script (watercarb), else goes to out_moon/; SMOON_* override for tests
outdir   = get(ENV, "SMOON_OUT", joinpath(@__DIR__, isdir(joinpath(@__DIR__, "Tables")) ? "Tables" : "out_moon"))
savemaps = !endswith(rstrip(outdir, '/'), "Tables")                  # .jls maps for figures, kept out of Tables/
fname    = "StagH2O_moon.dat"

# P-T grid: three sectors down to the core-mantle boundary (~4.7 GPa), nP x nT nodes each
nP       = parse(Int, get(ENV, "SMOON_NP", "100"))                     # pressure nodes per sector
nT       = parse(Int, get(ENV, "SMOON_NT", "100"))                     # temperature nodes per sector
Psec     = (1e-4, 0.5, 2.0, 5.0)                                       # GPa: crust and lid | upper mantle | deep mantle
Tsec     = ((250.0, 2200.0), (250.0, 2200.0), (250.0, 2200.0))         # K, from the surface (Tsurf_eqm 260 K) past the solidus
Teq_min  = 800.0                                                       # K, assemblage equilibrated at max(T, Teq_min): frozen lid

# Compositions: the endmembers of LunarWater.jl (mantle at p = 0, basalt at p = 1) and mixes in between, wt%
db       = "ig"                                                        # igneous set: every lunar major oxide but MnO
Xox      = ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "O", "Cr2O3", "H2O"]
mantle   = "Taylor"                                                    # bulk silicate Moon: "Taylor" (1999) | "Warren" (2005)
XH       = Dict("Taylor" => [43.4, 6.0, 4.5, 32.0, 13.0, 0.0, 0.09, 0.3, 0.0, 0.61, 0.0],
                "Warren" => [46.8, 3.87, 3.06, 36.0, 9.2, 0.0, 0.05, 0.18, 0.0, 0.44, 0.0])[mantle]  # MnO dropped; O = 0: no Fe3+
XB       = [50.47, 14.70, 11.39, 7.58, 10.43, 0.0, 2.79, 0.0, 0.0, 0.0, 0.0]  # LunarWater.jl basalt (its 1.68 MnO dropped)
Xmix     = [0.05, 0.1, 0.2, 0.5, 0.65, 0.75, 0.85, 0.95]               # basalt fractions of the mix blocks, file order

# Water in solids
hydrous  = false                                                       # false: dry bulk, so no hydrous phase; true: aH2O = 1
pl_wt    = 0.0                                                         # plagioclase water at saturation (wt%), no calibrated law
# ==============================================================================================

sectors = Tuple((LinRange(Psec[q], Psec[q+1], nP), LinRange(Tsec[q]..., nT)) for q in 1:3)
pX      = vcat(0.0, 1.0, Xmix)
nEM     = length(pX)
mkpath(outdir)
d       = Initialize_MAGEMin(db, verbose=false, buffer=hydrous ? "aH2O" : "")
Xc      = [p .* XB .+ (1 - p) .* XH for p in pX]
hydrous && (Xc = [(x = copy(x); x[end] = 20.0; x) for x in Xc])        # excess water under the aH2O buffer

# ---- Minimisation: one call over every (P, T, sector, composition), melt (and, dry, fluid) removed
t0  = time()
idx = [(i, j, q, c) for q in 1:3 for c in 1:nEM for j in 1:nT for i in 1:nP]
Pg  = [sectors[q][1][i] for (i, j, q, c) in idx]; Tg = [sectors[q][2][j] for (i, j, q, c) in idx]
rml = remove_phases(hydrous ? ["liq"] : ["liq", "fl"], db)
out = multi_point_minimization(10 .* Pg, max.(Tg, Teq_min) .- 273.15, d, X=[Xc[t[4]] for t in idx], Xoxides=Xox, sys_in="wt",
                               B=hydrous ? ones(length(idx)) : nothing, rm_list=rml, name_solvus=true, progressbar=false)
Finalize_MAGEMin(d)
@printf("%d points in %.0f s, %d not converged\n", length(out), time() - t0, count(o -> o.status > 1, out))

# ---- Capacity: hydrous phases (if any) plus NAM laws at the node's own (P, T); plagioclase from pl_wt
L = [merge(l, (pl=pl_wt,)) for l in nam_laws(Pg, Tg, nam_fits(); alpv=false)]
n = length(out); cap = zeros(n); hyd = zeros(n); host = fill("", n); unknown = Set{String}()
assemble!(cap, hyd, host, unknown, out, L; alpv=false)
isempty(unknown) || error("unclassified phases $(collect(unknown)), add them to NAMLAW or DRY in sfO2_fns.jl")
S = fill(NaN, nP, nT, 3, nEM); Shyd = similar(S); Shost = fill("", nP, nT, 3, nEM)
for (k, (i, j, q, c)) in enumerate(idx); S[i, j, q, c] = cap[k]; Shyd[i, j, q, c] = hyd[k]; Shost[i, j, q, c] = host[k]; end

# ---- Export: StagH2O layout, as the Earth table (header nP nT nEM + mixes; per sector basalt, residue, mixes)
savemaps && serialize(joinpath(outdir, "sH2O_moon_maps.jls"), (sectors=sectors, pX=pX, S=S, Shyd=Shyd, Shost=Shost, XB=XB, XH=XH, Xox=Xox,
                                                               Teq_min=Teq_min, hydrous=hydrous, pl_wt=pl_wt))
open(joinpath(outdir, fname), "w") do io
    println(io, nP, " ", nT, " ", nEM, " ", join(Xmix, " "), "\n")
    for (P, T) in sectors; println(io, P[1], " ", P[end], " ", T[1], " ", T[end]); end
    println(io)
    for q in 1:3, c in vcat([2, 1], 3:nEM)
        for i in 1:nP; for j in 1:nT; print(io, fmtv(max(1e-10, S[i, j, q, c]), 4), " "); end; println(io, ""); end
        println(io, "")
    end
end
@printf("Wrote %s: capacity %.2e to %.2e wt%%\n", joinpath(outdir, fname), minimum(S), maximum(S))
