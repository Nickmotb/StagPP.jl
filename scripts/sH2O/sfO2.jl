# Stand-alone generator of the StagYY lookup tables: StagH2O (storage capacity), StagfO2 + StagMetal, StagfO2step,
#  StagSdens + StagMdens. Kernels in sfO2_fns.jl; method in StagPP.jl scripts/sH2O/sH2O_math.md
using MAGEMin_C, Printf, Serialization
include(joinpath(@__DIR__, "sfO2_fns.jl"))

# ================================ INPUTS (production values) ================================
# Run control: tables overwrite Tables/ beside the script (watercarb), else go to out/; SFO2_* variables override for tests
todo        = Symbol.(split(get(ENV, "SFO2_DO", "sH2O_um,sH2O_deep,fO2,step,dens"), ","))  # passes to run
outdir      = get(ENV, "SFO2_OUT", joinpath(@__DIR__, isdir(joinpath(@__DIR__, "Tables")) ? "Tables" : "out"))
savemaps    = !endswith(rstrip(outdir, '/'), "Tables")         # .jls maps for sH2O_plot.jl, kept out of Tables/
chunk       = 200_000                                          # points per MAGEMin call in the fO2 pass

# P-T grid: three sectors (upper mantle | transition zone | lower mantle), nP x nT nodes each
nP          = parse(Int, get(ENV, "SFO2_NP", "100"))           # pressure nodes per sector
nT          = parse(Int, get(ENV, "SFO2_NT", "100"))           # temperature nodes per sector
Psec        = (1e-4, 7.0, 25.0, 135.0)                         # GPa, sector bounds
Tsec        = ((500.0, 2000.0), (700.0, 4000.0), (700.0, 4000.0))  # K, per sector

# Compositions: harzburgite, basalt and mixes p XB + (1 - p) XH (wt%), shared by every composition table
Xox         = ["SiO2", "MgO", "FeO", "CaO", "Al2O3", "Na2O", "Cr2O3", "H2O"]
XH          = [43.43, 45.93, 8.34, 0.9, 1.0, 0.01, 0.3, 80.0]     # harzburgite (Stixrude & Lithgow-Bertelloni 2024)
XB          = [50.42, 9.77, 7.1, 12.54, 16.8, 2.23, 0.07, 800.0]  # basalt; H2O is the excess of the ume/mbe runs
fO2_Xmix    = [0.05, 0.1, 0.2, 0.5, 0.65, 0.75, 0.85, 0.95]    # basalt fractions of the mix blocks, file order
sys_in      = "wt"

# Storage capacity (StagH2O)
Rsol        = 0.01                                             # Fe3+/ΣFe of the capacity bulks
pmbe        = 0.65                                             # basalt fraction from which a composition runs on mbe, not ume
Phand       = (5.0, 6.0)                                       # GPa, ume compositions hand over to sb24 + DHMS
Pmbe        = (9.5, 10.0)                                      # GPa, mbe compositions hand over to sb24 (lawsonite limit)
Al_pv       = true                                             # bridgmanite water from its Al2O3 (Litasov 2003); false: 30 ppm
dhms        = true                                             # DHMS term (Komabayashi & Omori 2006 MSH, Brodholt & Wood 1993 H2O)
PDmax, PHmax = 50.0, 60.0                                      # GPa: no phase D beyond (Nishi 2014), no DHMS beyond (Ohtani 2014)
TDmax       = 1873.0                                           # K, K&O calibration limit of the DHMS term

# fO2 and native iron (StagfO2, StagMetal): signed Fe3+/ΣFe axis, asinh to the knee then logit
nR          = 36                                               # R nodes
Rrange      = (-1.96, 0.95)                                    # R range (< 0: metal branch)
Rknee, nRhi = 0.4, 6                                           # asinh -> logit knee, logit nodes above it
R0          = 2e-3                                             # asinh scale

# Silica-saturation step (StagfO2step)
step_grid   = Tuple(parse.(Int, split(get(ENV, "SFO2_STEP", "40,16"), ",")))  # P x T nodes per sector
step_R      = [-1.0, -0.2, 0.0, 0.01, 0.03, 0.08, 0.15, 0.3, 0.5, 0.8]

# Hydrous density corrections (StagSdens, StagMdens): 50 x 50 x 50 over P to 50 GPa, T to 3000 K, H2O to Hmax_wt
Hmax_wt, melt_H2O_model, melt_floor = 8.0, :ideal, true        # :ideal = shipped StagMdens (:JK: V_H2O <= 0 > 5 GPa)
# ==============================================================================================

# ---- Axes, compositions, oxide lists
sectors  = Tuple((LinRange(Psec[q], Psec[q+1], nP), LinRange(Tsec[q]..., nT)) for q in 1:3)
(Pum, Tum), (Ptz, Ttz), (Plm, Tlm) = sectors
Rdom     = Rhybridspace(Rrange[1], Rrange[2], nR; R0=R0, Rknee=Rknee, nhi=nRhi)
Xf       = vcat([XH, XB], [x.*XB .+ (1.0-x).*XH for x in fO2_Xmix])       # every table: harz, basalt, mixes
nEM      = length(Xf)
pX       = vcat(0.0, 1.0, fO2_Xmix)                                      # their basalt fractions
Xox24    = ["SiO2", "MgO", "FeO", "CaO", "Al2O3", "Na2O", "Cr2O3", "O"]
Xox24f   = oxidized_Xox(Xox24, "Fe_O")                                   # sb24 is handed Fe + O
dummy()  = (X=zeros(length(Xox24)), Xox=copy(Xox24), mm=get_Xoxmm(Xox24))
sb24bulk(X, R) = oxidize_bulk(X, R, dummy(); wt_in=(sys_in=="wt"), wt_out=(sys_in=="wt"), FeFormat="Fe_O", oldXox=Xox)
mkpath(outdir)

S     = fill(NaN, nP, nT, 3, nEM)          # capacity wt% [P, T, sector, composition]
Shyd  = fill(NaN, nP, nT, 3, nEM)          # part held by hydrous phases
Shost = fill("", nP, nT, 3, nEM)           # phase holding most water
F     = fill(NaN, nP, nT, nR, nEM, 3)      # log10 fO2 [P, T, R, composition, sector]
M     = fill(NaN, nP, nT, nR, nEM, 3)      # native-iron weight fraction, same layout
if :sH2O_um in todo || :sH2O_deep in todo
    t0 = time(); Ffit = nam_fits(); Afit, Kfit = msh_setup(Pmax=PHmax); ph = Kfit[end]
    @printf("NAM fits + DHMS setup: %.1f s (phase H: ΔfH %.2f kJ, S %.2f J/K, V0 %.4f J/bar, K' %.2f)\n",
            time() - t0, ph.H/1e3, ph.S, ph.V0, ph.Kp)
    # ume/mbe runs: bulk with excess water, pure-H2O fluid (aH2O buffer), melt and DEW fluid removed; mbe has no Cr2O3
    function um_run(db, X, P, T)
        ox = vcat(Xox24, "H2O")
        Xb = vcat(oxidize_bulk(X, Rsol, dummy(); wt_in=true, wt_out=true, FeFormat="FeO_O", oldXox=Xox), X[end])
        db == "mbe" && (keep = ox .!= "Cr2O3"; ox = ox[keep]; Xb = Xb[keep])
        inf = retrieve_solution_phase_information(db)
        rml = remove_phases([n for n in ("liq", "DEW") if n in inf.ss_name || n in inf.data_pp], db)
        du  = Initialize_MAGEMin(db, verbose=false, buffer="aH2O")
        ou  = multi_point_minimization(10 .* P, T .- 273.15, du, X=[Xb for _ in P], Xoxides=ox, sys_in=sys_in,
                                       B=ones(length(P)), rm_list=rml, name_solvus=true, progressbar=false)
        Finalize_MAGEMin(du)
        return ou
    end
end

# ---- Storage capacity, upper mantle: ume below pmbe basalt fraction, mbe from it
if :sH2O_um in todo
    Pg = vec([Pum[i] for i in 1:nP, j in 1:nT]); Tg = vec([Tum[j] for i in 1:nP, j in 1:nT])
    L  = nam_laws(Pg, Tg, Ffit; alpv=Al_pv)
    for c in 1:nEM
        db = pX[c] >= pmbe ? "mbe" : "ume"; tc = time(); unk = Set{String}()
        ou = um_run(db, Xf[c], Pg, Tg)
        assemble!(view(S, :, :, 1, c), view(Shyd, :, :, 1, c), view(Shost, :, :, 1, c), unk, ou, L; alpv=Al_pv)
        isempty(unk) || error("UM p = $(pX[c]) ($db): unclassified phases $(collect(unk)), add them to NAMLAW or DRY")
        @printf("sH2O UM   p = %.2f %s: %d points in %.1f s, %d not converged\n", pX[c], db, length(ou), time() - tc, count(o -> o.status > 1, ou))
    end
    # Over Phand the ume compositions blend linearly (wt%) into sb24 + DHMS, held beyond: they leave ume before
    #  its brucite field (no phase A there). mbe compositions stay on mbe up to the TZ sector
    rix = [(i, j, c) for c in 1:nEM if pX[c] < pmbe for j in 1:nT for i in 1:nP if Pum[i] > Phand[1]]
    Pr  = [Pum[i] for (i, j, c) in rix]; Tr = [Tum[j] for (i, j, c) in rix]
    Xc  = [sb24bulk(Xf[c], Rsol) for c in 1:nEM]
    Lr  = nam_laws(Pr, Tr, Ffit; alpv=Al_pv)
    d   = Initialize_MAGEMin("sb24", verbose=false)
    out = multi_point_minimization(10 .* Pr, Tr .- 273.15, d, X=[Xc[t[3]] for t in rix], Xoxides=Xox24f, sys_in=sys_in,
                                   name_solvus=true, progressbar=false)
    Finalize_MAGEMin(d)
    n = length(out); cap = zeros(n); hyd = zeros(n); host = fill("", n); unknown = Set{String}(); sub = zeros(4, n)
    assemble!(cap, hyd, host, unknown, out, Lr; sub=sub, alpv=Al_pv)
    isempty(unknown) || error("UM hand-over: unclassified phases $(collect(unknown)), add them to NAMLAW or DRY")
    dhms && dhms!(cap, hyd, host, sub, Pr, Tr, Lr, Afit, Kfit; Pmax=PHmax, Tmax=TDmax, PDmax=PDmax)
    for (k, (i, j, c)) in enumerate(rix)
        w = min((Pum[i] - Phand[1])/(Phand[2] - Phand[1]), 1.0)
        S[i, j, 1, c] = (1 - w)*S[i, j, 1, c] + w*cap[k]; Shyd[i, j, 1, c] = (1 - w)*Shyd[i, j, 1, c] + w*hyd[k]
        w >= 0.5 && (Shost[i, j, 1, c] = host[k])
    end
    @printf("sH2O UM   hand-over: %d sb24%s points above %.2f GPa\n", n, dhms ? " + DHMS" : "", Phand[1])
end

# ---- Storage capacity, transition zone + lower mantle: sb24 NAM sum + DHMS, one call for both sectors and all blocks
if :sH2O_deep in todo
    idx = [(i, j, q, c) for q in 2:3 for c in 1:nEM for j in 1:nT for i in 1:nP]
    Pg  = [sectors[q][1][i] for (i, j, q, c) in idx]; Tg = [sectors[q][2][j] for (i, j, q, c) in idx]
    Xc  = [sb24bulk(Xf[c], Rsol) for c in 1:nEM]
    L   = nam_laws(Pg, Tg, Ffit; alpv=Al_pv)
    d   = Initialize_MAGEMin("sb24", verbose=false); t0 = time()
    out = multi_point_minimization(10 .* Pg, Tg .- 273.15, d, X=[Xc[t[4]] for t in idx], Xoxides=Xox24f, sys_in=sys_in,
                                   name_solvus=true, progressbar=false)
    Finalize_MAGEMin(d)
    n = length(out); cap = zeros(n); hyd = zeros(n); host = fill("", n); unknown = Set{String}(); sub = zeros(4, n)
    assemble!(cap, hyd, host, unknown, out, L; sub=sub, alpv=Al_pv)
    isempty(unknown) || error("TZ/LM: unclassified phases $(collect(unknown)), add them to NAMLAW or DRY")
    dhms && dhms!(cap, hyd, host, sub, Pg, Tg, L, Afit, Kfit; Pmax=PHmax, Tmax=TDmax, PDmax=PDmax)
    for (k, (i, j, q, c)) in enumerate(idx); S[i, j, q, c] = cap[k]; Shyd[i, j, q, c] = hyd[k]; Shost[i, j, q, c] = host[k]; end
    @printf("sH2O TZ+LM: %d points in %.1f s, %d not converged\n", n, time() - t0, count(o -> o.status > 1, out))
    # mbe compositions below Pmbe[2]: mbe as in the UM (continuous 7 GPa seam), blended linearly into sb24 over Pmbe
    ix = [(i, j) for j in 1:nT for i in 1:nP if Ptz[i] < Pmbe[2]]
    Pb = [Ptz[i] for (i, j) in ix]; Tb = [Ttz[j] for (i, j) in ix]; Lb = nam_laws(Pb, Tb, Ffit; alpv=Al_pv)
    cb = zeros(length(ix)); hb = zeros(length(ix)); sb = fill("", length(ix))
    for c in 1:nEM
        pX[c] < pmbe && continue
        assemble!(cb, hb, sb, unknown, um_run("mbe", Xf[c], Pb, Tb), Lb; alpv=Al_pv)
        isempty(unknown) || error("TZ mbe p = $(pX[c]): unclassified phases $(collect(unknown)), add them to NAMLAW or DRY")
        for (k, (i, j)) in enumerate(ix)
            w = clamp((Ptz[i] - Pmbe[1])/(Pmbe[2] - Pmbe[1]), 0.0, 1.0)
            S[i, j, 2, c] = (1 - w)*cb[k] + w*S[i, j, 2, c]; Shyd[i, j, 2, c] = (1 - w)*hb[k] + w*Shyd[i, j, 2, c]
            w < 0.5 && (Shost[i, j, 2, c] = sb[k])
        end
    end
    @printf("sH2O TZ   mbe: %d points per composition below %.2f GPa\n", length(ix), Pmbe[2])
end

# ---- fO2 and native iron over (P, T, R, composition, sector), in chunks of the column-major (i, j, r, k, q) order
if :fO2 in todo
    Xb = [sb24bulk(Xf[k], Rdom[r]) for k in 1:nEM, r in 1:nR]
    CI = CartesianIndices((nP, nT, nR, nEM, 3))
    d  = Initialize_MAGEMin("sb24", verbose=false); t0 = time()
    for c0 in 1:chunk:length(CI)
        sl = CI[c0:min(c0 + chunk - 1, end)]
        oc = multi_point_minimization([10*sectors[t[5]][1][t[1]] for t in sl], [sectors[t[5]][2][t[2]] - 273.15 for t in sl], d,
                                      X=[Xb[t[4], t[3]] for t in sl], Xoxides=Xox24f, name_solvus=true, progressbar=false, sys_in=sys_in)
        for (n, t) in enumerate(sl); F[t] = oc[n].fO2; M[t] = metal_wt_frac(oc[n]); end
        @printf("fO2: %d / %d points, %.0f s\n", c0 + length(sl) - 1, length(CI), time() - t0)
    end
    Finalize_MAGEMin(d)
    for (q, nm) in enumerate(("upper mantle", "transition zone", "lower mantle")), (A, what) in ((F, "fO2"), (M, "metal"))
        bad = [count(!isfinite, view(A, :, :, r, :, q)) for r in 1:nR]
        sum(bad) > 0 && @warn "$(sum(bad)) non-finite $what node(s) in the $nm maps at R = $(round.(Rdom[bad .> 0], digits=5))"
    end
    nfix = sum(enforce_R_monotone!(view(F, :, :, :, :, q)) for q in 1:3)
    nfix > 0 && println("Repaired $nfix non-monotonic fO2 column(s) along R (isotonic).")
    # Non-finite nodes (failed minimisations) take the mean of their finite P-T neighbours, sweep by sweep within each
    #  (R, composition, sector) slice: the fill the shipped tables carry
    nfill = count(!isfinite, F) + count(!isfinite, M); nbr = CartesianIndex.(((-1, 0), (1, 0), (0, -1), (0, 1)))
    for A in (F, M), q in 1:3, k in 1:nEM, r in 1:nR
        Sl = view(A, :, :, r, k, q); bad = findall(!isfinite, Sl)
        while !isempty(bad)
            new = [(v = [Sl[c + d] for d in nbr if checkbounds(Bool, Sl, c + d) && isfinite(Sl[c + d])]; isempty(v) ? NaN : sum(v)/length(v)) for c in bad]
            all(isnan, new) && break
            Sl[bad] .= new; bad = findall(!isfinite, Sl)
        end
    end
    nfill > 0 && println("Filled $nfill non-finite fO2/metal node(s) from their P-T neighbours.")
    M .= max.(0.0, M)
end

# ---- Silica-saturation step in basalt fraction: pa = last p with an olivine-family phase or mw, pb = first p with free silica.
#  Both searches share the p = 0.3 and 1 end members; only nodes that bracket a change are bisected (10 halvings)
if :step in todo
    sect = Tuple((Psec[q], Psec[q+1], Tsec[q]...) for q in 1:3)
    nPc, nTc = step_grid; nK = length(step_R)
    pts = [(q, i, j, k) for q in 1:3 for k in eachindex(step_R) for i in 1:nPc for j in 1:nTc]; N = length(pts)
    Pq  = [10*max(sect[q][1] + (sect[q][2]-sect[q][1])*(i-1)/(nPc-1), 1e-4) for (q,i,j,k) in pts]
    Tq  = [sect[q][3] + (sect[q][4]-sect[q][3])*(j-1)/(nTc-1) - 273.15 for (q,i,j,k) in pts]
    dstep = Initialize_MAGEMin("sb24", verbose=false); t0 = time()
    mini(pv, id) = multi_point_minimization(Pq[id], Tq[id], dstep, Xoxides=Xox24f, sys_in=sys_in, name_solvus=true, progressbar=false,
                                            X=[sb24bulk(pv[m].*XB .+ (1-pv[m]).*XH, step_R[pts[id[m]][4]]) for m in eachindex(id)])
    olo = mini(fill(0.3, N), 1:N); ohi = mini(fill(1.0, N), 1:N)
    stepab = Vector{Vector{Float64}}()
    for (names, grow) in ((("ol","wa","ri","wad","ring","mw"), true), (("coe","st","qtz","seif","apbo"), false))
        has(o) = any(x -> lowercase(x) in names, o.ph)
        ok = [n for n in 1:N if (grow ? (has(olo[n]) && !has(ohi[n])) : (!has(olo[n]) && has(ohi[n])))]
        lo = fill(0.3, length(ok)); hi = fill(1.0, length(ok))
        for _ in 1:10
            isempty(ok) && break
            mid = 0.5 .* (lo .+ hi); o = mini(mid, ok)
            for m in eachindex(ok); (has(o[m]) == grow) ? (lo[m] = mid[m]) : (hi[m] = mid[m]); end
        end
        v = fill(NaN, N); v[ok] .= 0.5 .* (lo .+ hi); push!(stepab, v)
    end
    Finalize_MAGEMin(dstep)
    open(joinpath(outdir, "StagfO2step.dat"), "w") do io
        println(io, nPc, " ", nTc, " ", nK, "\n")
        for q in sect; println(io, join(q, " ")); end
        println(io, join(step_R, " "), "\n")
        for v in stepab, k in 1:nK, q in 1:3
            base = ((q-1)*nK + (k-1))*nPc*nTc
            for i in 1:nPc; println(io, join([@sprintf("%.4f", v[base + (i-1)*nTc + j]) for j in 1:nTc], " ")); end
            println(io)
        end
    end
    @printf("Step table: %.0f s\n", time() - t0)
end

# ---- Solid and melt H2O density corrections
if :dens in todo
    t0 = time()
    SΔρ, MΔρ, ρP, ρT, ρH = PT_H2O_ρ(Psec[1], 50., Tsec[1][1], 3000.; nP=50, nT=50, nH=50,
                                    Hmax_wt=Hmax_wt, melt_H2O_model=melt_H2O_model, melt_floor=melt_floor)
    for (fn, A) in (("StagSdens.dat", SΔρ), ("StagMdens.dat", MΔρ))
        open(joinpath(outdir, fn), "w") do io
            println(io, length(ρP), " ", length(ρT), " ", length(ρH), "\n")
            println(io, ρP[1], " ", ρP[end], " ", ρT[1], " ", ρT[end], " ", round(ρH[1], digits=4), " ", round(ρH[end], digits=4), "\n")
            for slot in 1:2, h in eachindex(ρH)
                for i in eachindex(ρP); for j in eachindex(ρT); print(io, fmtv(A[i, j, h, slot], 4), " "); end; println(io, ""); end
                println(io, "")
            end
        end
    end
    @printf("Density tables: %.0f s\n", time() - t0)
end

# ---- Export. Same layout and number format as the reference writers
header(io) = for (P, T) in sectors; println(io, P[1], " ", P[end], " ", T[1], " ", T[end]); end
if :sH2O_um in todo && :sH2O_deep in todo
    savemaps && serialize(joinpath(outdir, "sH2O_maps.jls"), (sectors=sectors, pX=pX, S=S, Shyd=Shyd, Shost=Shost))
    open(joinpath(outdir, "StagH2O.dat"), "w") do io
        println(io, nP, " ", nT, " ", nEM, " ", join(fO2_Xmix, " "), "\n"); header(io); println(io)
        for q in 1:3, c in vcat([2, 1], 3:nEM)                         # per sector: basalt, harzburgite, mixes (as StagfO2)
            for i in 1:nP; for j in 1:nT; print(io, fmtv(max(1e-10, S[i, j, q, c]), 4), " "); end; println(io, ""); end
            println(io, "")
        end
    end
end
if :fO2 in todo
    for (fn, A, dg) in (("StagfO2.dat", F, 4), ("StagMetal.dat", M, 6))
        open(joinpath(outdir, fn), "w") do io
            println(io, nP, " ", nT, " ", nR, " ", nEM, " ", join(fO2_Xmix, " "), "\n"); header(io)
            for r in 1:nR; print(io, Rdom[r], " "); end; println(io, "\n")
            for r in 1:nR, q in 1:3, k in vcat([2, 1], 3:nEM)          # per R node and sector: basalt, harzburgite, mixes
                for i in 1:nP; for j in 1:nT; print(io, fmtv(A[i, j, r, k, q], dg), " "); end; println(io, ""); end
                println(io, "")
            end
        end
    end
end
