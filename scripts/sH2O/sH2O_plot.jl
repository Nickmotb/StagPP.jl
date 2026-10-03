# Figures of sH2O_math.md into md_figs/: the model's kernels plus the tables of a finished sfO2.jl run.
#  SH2O_RUN (default out/) holds the run; its variants/nodhms and variants/noalpv hold sH2O-only reruns with dhms / Al_pv off
using MAGEMin_C, CairoMakie, Serialization, Printf, Statistics
include(joinpath(@__DIR__, "sfO2_fns.jl"))
run  = get(ENV, "SH2O_RUN", joinpath(@__DIR__, "out")); figs = joinpath(@__DIR__, "md_figs"); mkpath(figs)
# Maps hold every composition (pX); most figures show harzburgite, basalt and the 20 % basalt mix
blk3(m) = (ib = [1, 2, findfirst(≈(0.2), m.pX)]; (sectors=m.sectors, S=m.S[:, :, :, ib], Shyd=m.Shyd[:, :, :, ib], Shost=m.Shost[:, :, :, ib]))
rf   = deserialize(joinpath(run, "sH2O_maps.jls")); rf0 = deserialize(joinpath(run, "variants/nodhms/sH2O_maps.jls"))
r, r0 = blk3(rf), blk3(rf0); r1 = blk3(deserialize(joinpath(run, "variants/noalpv/sH2O_maps.jls")))
Ffit = nam_fits(); A, K = msh_setup()
save_fig(name, f) = (save(joinpath(figs, name), f, px_per_unit=2); println("wrote ", name))

ink = colorant"#0b0b0b"; muted = colorant"#52514e"; faint = colorant"#e8e7e2"
seq = cgrad([colorant"#f3f7fd", colorant"#cde2fb", colorant"#9ec5f4", colorant"#6da7ec", colorant"#3987e5", colorant"#256abf", colorant"#184f95", colorant"#0d366b"])
div = cgrad([colorant"#184f95", colorant"#6da7ec", colorant"#e8e7e2", colorant"#f0a07a", colorant"#b8431a"])
cat3 = (colorant"#2a78d6", colorant"#eb6834", colorant"#1baf7a"); cat4 = (cat3..., colorant"#8a5cd6")
ramp4 = (colorant"#9ec5f4", colorant"#3987e5", colorant"#184f95", colorant"#0b0b0b")
set_theme!(fontsize=12, textcolor=ink, Axis=(xgridvisible=false, ygridvisible=false, titlealign=:left, titlefont=:bold, titlesize=12),
           backgroundcolor=:white)
bnames = ("Harzburgite", "Basalt", "Mix (20 % basalt)"); snames = ("Upper mantle", "Transition zone", "Lower mantle")
# n passes of a [1 2 1] x [1 2 1] filter, for contours without the grid staircase
function smooth(I; n=2)
    for _ in 1:n
        J = similar(I)
        for i in axes(I, 1), j in axes(I, 2)
            sw = 0.0; sv = 0.0
            for di in -1:1, dj in -1:1
                checkbounds(Bool, I, i + di, j + dj) || continue
                k = (2 - abs(di))*(2 - abs(dj)); sw += k; sv += k*I[i+di, j+dj]
            end
            J[i, j] = sv/sw
        end
        I = J
    end
    return I
end
# Phase boundaries: the 0.5 contour of each phase's smoothed indicator
boundaries!(a, P, T, H; color=ink, linewidth=0.9) = for h in unique(H)
    h == "" || contour!(a, collect(P), collect(T), smooth(Float64.(H .== h)), levels=[0.5], color=color, linewidth=linewidth)
end
function labels!(a, P, T, H; nmin=25, fs=10)
    for h in unique(H)
        (h == "" || count(==(h), H) < nmin) && continue
        id = findall(==(h), H); ii = round(Int, median(first.(Tuple.(id)))); jj = round(Int, median(last.(Tuple.(id))))
        H[ii, jj] == h || ((ii, jj) = Tuple(id[argmin([(Tuple(k)[1]-ii)^2 + (Tuple(k)[2]-jj)^2 for k in id])]))
        text!(a, P[ii], T[jj], text=h, align=(:center, :center), fontsize=fs, font=:bold, color=:white, strokecolor=ink, strokewidth=0.6)
    end
end
# Bilinear value of map M[:, :, q, c] at (P, T) on sector q's axes
function at(M, q, c, P, T)
    Pa, Ta = r.sectors[q]; x = clamp((P - Pa[1])/step(Pa) + 1, 1, length(Pa) - 1e-9); y = clamp((T - Ta[1])/step(Ta) + 1, 1, length(Ta) - 1e-9)
    i, j = floor(Int, x), floor(Int, y); u, v = x - i, y - j
    return (1-u)*(1-v)*M[i, j, q, c] + u*(1-v)*M[i+1, j, q, c] + (1-u)*v*M[i, j+1, q, c] + u*v*M[i+1, j+1, q, c]
end
sector(P) = P < 7.0 ? 1 : (P < 25.0 ? 2 : 3)
# Illustrative geotherms (K): mantle adiabat, warm and cold slab tops
geo = (("mantle adiabat", P -> P <= 25 ? 1600 + 12P : 1900 + 600(P - 25)/110),
       ("warm slab", P -> P <= 10 ? 700 + 60P : (P <= 25 ? 1300 + 15(P - 10) : 1525 + 6(P - 25))),
       ("cold slab", P -> P <= 10 ? 500 + 45P : (P <= 25 ? 950 + 15(P - 10) : 1175 + 6(P - 25))))

# ---- Table readers (layout of the sfO2.jl writers)
function read_fo2(fn)
    v = parse.(Float64, split(read(fn, String))); nP, nT, nR, nE = Int.(v[1:4]); o = 4 + nE - 2
    sec = [v[o+4(q-1)+1:o+4q] for q in 1:3]; o += 12; R = v[o+1:o+nR]; o += nR
    X = zeros(nP, nT, nR, nE, 3)
    for rr in 1:nR, q in 1:3, k in vcat([2, 1], 3:nE), i in 1:nP, j in 1:nT; o += 1; X[i, j, rr, k, q] = v[o]; end
    return (sec=sec, R=R, X=X, mix=v[5:4+nE-2])
end
function read_step(fn)
    v = parse.(Float64, split(read(fn, String))); nP, nT, nK = Int.(v[1:3]); sec = [v[3+4(q-1)+1:3+4q] for q in 1:3]; o = 15
    Rk = v[o+1:o+nK]; o += nK; X = zeros(nP, nT, nK, 3, 2)
    for s in 1:2, k in 1:nK, q in 1:3, i in 1:nP, j in 1:nT; o += 1; X[i, j, k, q, s] = v[o]; end
    return (sec=sec, R=Rk, X=X)
end
function read_dens(fn)
    v = parse.(Float64, split(read(fn, String))); nP, nT, nH = Int.(v[1:3]); b = v[4:9]; o = 9; X = zeros(nP, nT, nH, 2)
    for s in 1:2, h in 1:nH, i in 1:nP, j in 1:nT; o += 1; X[i, j, h, s] = v[o]; end
    return (P=LinRange(b[1], b[2], nP), T=LinRange(b[3], b[4], nT), H=LinRange(b[5], b[6], nH), X=X)
end

# ---- S01: which method fills which part of each block's P-T domain
f = Figure(size=(1400, 1000))
lanes = (("ume", colorant"#cfe8d9"), ("mbe", colorant"#fbdcc9"), ("sb24 NAM", colorant"#d6e6fb"), ("blend", colorant"#dcdad3"))
for (row, c) in enumerate((1, 3, 2)), (col, (P0, P1)) in enumerate(((0.0, 25.0), (25.0, 135.0)))
    a = Axis(f[row, col], title=col == 1 ? bnames[c] : "", xlabel=row == 3 ? "Pressure (GPa)" : "", ylabel=col == 1 ? "Temperature (K)" : "")
    box(x0, x1, y0, y1, cl) = (x1 > P0 && x0 < P1) && poly!(a, Rect(max(x0, P0), y0, min(x1, P1) - max(x0, P0), y1 - y0), color=cl, strokewidth=0)
    if c == 2
        box(0, 7, 500, 2000, lanes[2][2]); box(7, 9.5, 700, 4000, lanes[2][2]); box(9.5, 10, 700, 4000, lanes[4][2]); box(10, 135, 700, 4000, lanes[3][2])
    else
        box(0, 5, 500, 2000, lanes[1][2]); box(5, 6, 500, 2000, lanes[4][2]); box(6, 7, 500, 2000, lanes[3][2]); box(7, 135, 700, 4000, lanes[3][2])
    end
    Pd0 = max(P0, c == 2 ? 9.5 : 5.0)
    lines!(a, [Pd0, Pd0, min(P1, 60.0)], [700.0, 1873.0, 1873.0], color=ink, linewidth=2, linestyle=:dash)
    P1 >= 60 && lines!(a, [60.0, 60.0], [700.0, 1873.0], color=ink, linewidth=2, linestyle=:dash)
    P1 > 50 && lines!(a, [50.0, 50.0], [700.0, 1873.0], color=muted, linewidth=1.2, linestyle=:dot)
    lines!(a, [max(5.0, P0), min(60.0, P1)], [773.15, 773.15], color=cat3[1], linewidth=1.2, linestyle=:dash)
    lines!(a, [max(5.0, P0), min(60.0, P1)], [1873.15, 1873.15] .+ 25, color=cat3[1], linewidth=1.2, linestyle=:dash)
    vlines!(a, [7.0, 25.0], color=muted, linewidth=0.8)
    for (nm, g) in geo
        Pg = range(P0, P1, length=100); lines!(a, Pg, g.(Pg), color=muted, linewidth=1)
        col == 2 && text!(a, P1 - 2, g(P1 - 2), text=nm, align=(:right, :bottom), fontsize=9, color=muted)
    end
    xlims!(a, P0, P1); ylims!(a, 450, 4050)
end
Legend(f[4, 1:2], vcat([PolyElement(color=cl) for (nm, cl) in lanes], [LineElement(color=ink, linestyle=:dash, linewidth=2),
       LineElement(color=muted, linestyle=:dot), LineElement(color=cat3[1], linestyle=:dash), LineElement(color=muted)]),
       vcat([nm for (nm, cl) in lanes], ["DHMS domain (sb24 part, P ≤ 60 GPa, T ≤ 1873 K)", "no phase D beyond 50 GPa",
       "a(H2O) map, held flat outside 773-1873 K", "illustrative geotherms (S13)"]), orientation=:horizontal, nbanks=2, framevisible=false, labelsize=10)
save_fig("S01_domain.png", f)

# ---- S02: water fugacity (Pitzer & Sterner) and the DHMS water potential (Brodholt & Wood above 1 GPa)
Pf = range(0.1, 60, length=240); Tf = (800.0, 1200.0, 1600.0, 2000.0)
f = Figure(size=(1400, 480))
a1 = Axis(f[1, 1], yscale=log10, xlabel="Pressure (GPa)", ylabel="f(H2O) (GPa)", title="a  Pitzer & Sterner (1994) fugacity, used by the NAM laws")
a2 = Axis(f[1, 2], xlabel="Pressure (GPa)", ylabel="μ_BW − μ_PS (kJ/mol)", title="b  DHMS water potential against the Pitzer & Sterner one")
a3 = Axis(f[1, 3], xlabel="Pressure (GPa)", ylabel="V (cm³/mol)", title="c  Brodholt & Wood (1993) molar volume")
for (n, T) in enumerate(Tf)
    lines!(a1, Pf, [fH2O(P, T) for P in Pf], color=ramp4[n], linewidth=2, label="$(Int(T)) K")
    lines!(a2, Pf, [(mu_h2o(P, T) - (g_ko(KO_W, 0.0, T) + Rgas*T*log(1e4fH2O(P, T))))/1e3 for P in Pf], color=ramp4[n], linewidth=2)
    lines!(a3, Pf[Pf .>= 1], [bw_V(1e4P, T) for P in Pf[Pf .>= 1]], color=ramp4[n], linewidth=2)
end
axislegend(a1, position=:lt, framevisible=false); hlines!(a2, [0.0], color=muted, linewidth=0.8)
save_fig("S02_water.png", f)

# ---- S03: upper-mantle capacity and the share held by hydrous phases
f = Figure(size=(1400, 820)); P, T = r.sectors[1]
for (col, c) in enumerate((1, 3, 2))
    C = log10.(max.(r.S[:, :, 1, c], 1e-4)); H = copy(r.Shost[:, :, 1, c])
    a = Axis(f[1, col], title="$(bnames[c]): capacity and host", ylabel=col == 1 ? "Temperature (K)" : "")
    hm = heatmap!(a, collect(P), collect(T), C, colormap=seq, colorrange=(-3, 1.3), interpolate=true)
    boundaries!(a, P, T, H; color=(ink, 0.6))
    labels!(a, P, T, H; nmin=30)
    c == 2 || vlines!(a, [5.0, 6.0], color=muted, linestyle=:dot)
    col == 3 && Colorbar(f[1, 4], hm, label="log10 capacity (wt% of solid)", width=12)
    b = Axis(f[2, col], title="share held by hydrous phases", xlabel="Pressure (GPa)", ylabel=col == 1 ? "Temperature (K)" : "")
    hb = heatmap!(b, collect(P), collect(T), r.Shyd[:, :, 1, c] ./ r.S[:, :, 1, c], colormap=seq, colorrange=(0, 1), interpolate=true)
    col == 3 && Colorbar(f[2, 4], hb, label="hydrous / total", width=12)
end
save_fig("S03_upper_mantle.png", f)

# ---- S04: the NAM laws at a(H2O) = 1, against P at four T, with each law's calibration box
lawset = (("olivine (Dong 2020)", :ol, (0.0, 15.0), ((0.01, 14.2), (1200, 2400))), ("opx, pure En (Mierdel 2007)", :opx, (0.0, 15.0), ((0.01, 14.0), (500, 1500))),
          ("opx, Al term (Mierdel 2007)", :opxAl, (0.0, 15.0), ((0.01, 14.0), (1073, 1500))), ("diopside (Bromiley 2004)", :di, (0.0, 12.0), ((0.0, 10.0), (873, 973))),
          ("jadeite (Bromiley & Keppler 2004)", :jd, (0.0, 12.0), ((0.0, 10.0), (873, 973))), ("garnet (Lu & Keppler 1997)", :grt, (0.0, 25.0), ((0.01, 10.0), (700, 1500))),
          ("coesite", :coe, (0.0, 12.0), ((0.0, 10.0), (800, 1073))), ("wadsleyite (fit)", :wa, (12.0, 23.0), (Ffit.wa.box[1], Ffit.wa.box[2])),
          ("ringwoodite (fit); akimotoite = /21", :ri, (16.0, 27.0), (Ffit.ri.box[1], Ffit.ri.box[2])), ("HP clinoenstatite (fit)", :hpx, (6.0, 16.0), (Ffit.hpx.box[1], Ffit.hpx.box[2])),
          ("silica: st | CaCl2-st | α-PbO2 (> 90 GPa here)", :st, (10.0, 135.0), ((0.01, 60.01), (800, 1500))), ("post-perovskite on 30 ppm pv: Al-free | Al", :rNo, (90.0, 135.0), ((90.0, 135.0), (1000, 4000))))
f = Figure(size=(1400, 1050)); Tl = (1000.0, 1400.0, 1800.0, 2200.0)
for (n, (ttl, fld, (Pa, Pb), bx)) in enumerate(lawset)
    a = Axis(f[cld(n, 4), mod1(n, 4)], title=ttl, yscale=log10, xlabel=n > 8 ? "Pressure (GPa)" : "", ylabel=mod1(n, 4) == 1 ? "wt% at a(H2O) = 1" : "")
    vspan!(a, max(bx[1][1], Pa), min(bx[1][2], Pb), color=(faint, 0.8))
    Pl = range(Pa, Pb, length=160)
    for (m, T) in enumerate(Tl)
        L = [nam_laws(P, T, Ffit) for P in Pl]
        y = fld === :st ? [P > 90 ? l.apbo : l.st for (P, l) in zip(Pl, L)] : fld === :rNo ? [3e-3l.rNo for l in L] : [getfield(l, fld) for l in L]
        lines!(a, Pl, max.(y, 1e-6), color=ramp4[m], linewidth=2, label="$(Int(T)) K")
        fld === :ri && lines!(a, Pl, max.([l.ak for l in L], 1e-6), color=ramp4[m], linewidth=1, linestyle=:dash)
        fld === :rNo && lines!(a, Pl, max.([3e-3l.rAl for l in L], 1e-6), color=ramp4[m], linewidth=1, linestyle=:dash)
    end
    text!(a, 0.02, 0.02, text=@sprintf("box %.3g-%.3g GPa, %.0f-%.0f K", bx[1]..., bx[2]...), space=:relative, fontsize=9, color=muted)
end
Legend(f[4, 1:4], [LineElement(color=ramp4[m], linewidth=2) for m in 1:4], ["$(Int(T)) K" for T in Tl], orientation=:horizontal, framevisible=false,
       tellheight=true)
save_fig("S04_nam_laws.png", f)

# ---- S05: bridgmanite from its own Al2O3 (Al_pv) against Litasov et al. 2003, and its effect on the maps
law = bdg_h2o
Xox = ["SiO2", "MgO", "FeO", "CaO", "Al2O3", "Na2O", "Cr2O3", "H2O"]
XH = [43.43, 45.93, 8.34, 0.9, 1.0, 0.01, 0.3, 80.0]; XB = [50.42, 9.77, 7.1, 12.54, 16.8, 2.23, 0.07, 800.0]
Xox24 = ["SiO2", "MgO", "FeO", "CaO", "Al2O3", "Na2O", "Cr2O3", "O"]; Xox24f = oxidized_Xox(Xox24, "Fe_O")
dummy() = (X=zeros(length(Xox24)), Xox=copy(Xox24), mm=get_Xoxmm(Xox24))
sbulk(X) = oxidize_bulk(X[1:end-1], 0.01, dummy(); wt_in=true, wt_out=true, FeFormat="Fe_O", oldXox=Xox[1:end-1])
Xs = (XH, XB, 0.2 .* XB .+ 0.8 .* XH); Pb5 = collect(range(20, 135, length=93)); Tiso = (1600.0, 2200.0)
AL = fill(NaN, length(Pb5), 2, 3); d = Initialize_MAGEMin("sb24", verbose=false)
for c in 1:3, (t, T) in enumerate(Tiso)
    ou = multi_point_minimization(10Pb5, fill(T - 273.15, length(Pb5)), d, X=[sbulk(Xs[c]) for _ in Pb5], Xoxides=Xox24f, sys_in="wt",
                                  name_solvus=true, progressbar=false)
    for (k, o) in enumerate(ou), j in eachindex(o.SS_vec); o.ph[j] == "pv" && (AL[k, t, c] = 100o.SS_vec[j].Comp_wt[findfirst(==("Al2O3"), o.oxides)]); end
end
Finalize_MAGEMin(d)
f = Figure(size=(1400, 980))
a = Axis(f[1, 1], xlabel="Al₂O₃ in bridgmanite (wt%)", ylabel="H₂O in bridgmanite (ppm)", title="a  Bridgmanite law and its calibration data")
for (x0, x1, y0, y1, lab, lp) in ((0.0, 0.3, 80.0, 120.0, "MgSiO₃", (0.5, 170.0)), (4.0, 7.0, 1000.0, 1500.0, "Al-Mg-pv", (7.2, 1400.0)),
                                  (5.0, 6.0, 1400.0, 1800.0, "peridotite pv", (6.2, 1750.0)), (13.0, 17.0, 40.0, 110.0, "MORB pv", (13.0, 200.0)))
    poly!(a, Rect(x0, y0, x1 - x0, y1 - y0), color=(muted, 0.18), strokecolor=muted, strokewidth=1)
    text!(a, lp..., text=lab, fontsize=10, color=muted, align=(:left, :center))
end
xs = range(0, 26, length=300); lines!(a, xs, 1e4 .* law.(xs), color=ink, linewidth=2, label="Al_pv law")
hlines!(a, [30.0], color=muted, linestyle=:dash, linewidth=1.5, label="Al_pv = false, 30 ppm")
for c in 1:3
    v = filter(isfinite, AL[:, :, c]); lo, hi = extrema(v); band!(a, [lo, max(hi, lo + 0.15)], [0.0, 0.0], [2000.0, 2000.0], color=(cat3[c], 0.15))
    text!(a, (lo + min(hi, 26.0))/2, 1930 - 70(c == 3), text=split(bnames[c], " ")[1] * @sprintf(" %.0f-%.0f", lo, hi), color=cat3[c], fontsize=10, align=(:center, :top), font=:bold)
end
xlims!(a, 0, 26); ylims!(a, 0, 2000); axislegend(a, position=:rt, framevisible=false, margin=(0, 0, 40, 0))
b = Axis(f[1, 2], xlabel="Pressure (GPa)", ylabel="Al₂O₃ in bridgmanite (wt%)", title="b  Al in sb24 bridgmanite: garnet dissolves into it")
for c in 1:3, t in 1:2; lines!(b, Pb5, AL[:, t, c], color=cat3[c], linewidth=2, linestyle=t == 1 ? :dash : :solid); end
for c in 1:3; lines!(b, [NaN], [NaN], color=cat3[c], linewidth=2, label=bnames[c]); end
lines!(b, [NaN], [NaN], color=muted, linestyle=:dash, label="1600 K"); lines!(b, [NaN], [NaN], color=muted, label="2200 K")
axislegend(b, position=:lt, framevisible=false)
c5 = Axis(f[2, 1], xlabel="Pressure (GPa)", ylabel="H₂O in bridgmanite (ppm)", yscale=log10, yticks=[10, 30, 100, 300, 1000, 3000], title="c  H₂O held by bridgmanite at saturation")
for c in 1:3, t in 1:2; lines!(c5, Pb5, 1e4 .* law.(AL[:, t, c]), color=cat3[c], linewidth=2, linestyle=t == 1 ? :dash : :solid); end
hlines!(c5, [30.0], color=muted, linestyle=:dot, linewidth=1.5); ylims!(c5, 10, 3000)
e = Axis(f[2, 2], xlabel="Pressure (GPa)", ylabel="Storage capacity (wt%)", yscale=log10, yticks=([1e-3, 1e-2, 1e-1, 1.0], ["0.001", "0.01", "0.1", "1"]),
         title="d  Rock capacity along 2200 K: Al_pv off (thin) and on (thick)")
for c in 1:3
    PP = Float64[]; So = Float64[]; Sn = Float64[]
    for q in 2:3; P, T = r.sectors[q]; j = argmin(abs.(T .- 2200)); append!(PP, P); append!(So, r1.S[:, j, q, c]); append!(Sn, r.S[:, j, q, c]); end
    lines!(e, PP, So, color=(cat3[c], 0.6), linewidth=1); lines!(e, PP, Sn, color=cat3[c], linewidth=2.5)
end
vlines!(e, [25.0], color=muted, linestyle=:dot); ylims!(e, 1e-3, 3)
save_fig("S05_bridgmanite.png", f)

# ---- S06: K&O (2006) benchmark: our MSH minimiser at Mg/Si = 1.4 and the offsets of 50 digitised reaction checkpoints
ix = Dict(n => i for (i, n) in enumerate(KO_NAME))
Gr(n, P, T) = n == "W" ? mu_h2o(P, T) : g_ko(KO[ix[n]], P, T)
RX = Dict(
 3 => ([(1,"phA"),(3,"hen")], [(5,"fo"),(3,"W")]), 4 => ([(1,"phA"),(10,"fo")], [(3,"cHu")]), 5 => ([(1,"cHu"),(1,"hen")], [(5,"fo"),(1,"W")]),
 6 => ([(1,"phA"),(2,"hen")], [(1,"cHu"),(2,"W")]), 11 => ([(5723,"phA"),(22250,"hen"),(17597,"W")], [(26743,"phE")]),
 13 => ([(500,"phE")], [(107,"cHu"),(202,"hen"),(543,"W")]), 14 => ([(2900,"phE")], [(3103,"fo"),(551,"hen"),(3770,"W")]),
 22 => ([(13126,"phA"),(24881,"phD"),(152,"W")], [(52035,"phE")]), 42 => ([(73410,"shB"),(173062,"W")], [(103449,"phA"),(8438,"phD")]),
 43 => ([(16730,"phE")], [(20342,"hyWad"),(738,"hen"),(19308,"W")]), 51 => ([(7989,"shB"),(20422,"phD"),(18964,"W")], [(44630,"phE")]),
 53 => ([(42630,"phE")], [(52834,"hyWad"),(880,"stv"),(49079,"W")]), 55 => ([(16288,"shB"),(38172,"stv")], [(19265,"phE"),(62762,"hyWad")]),
 62 => ([(7724,"shB"),(17913,"stv")], [(41085,"hyWad"),(10518,"W")]), 68 => ([(56400,"hyRin"),(30174,"stv")], [(76704,"aki"),(3102,"shB")]),
 76 => ([(4722,"shB")], [(14309,"hyRin"),(19890,"per"),(7870,"W")]), 78 => ([(11697,"shB"),(53018,"pv")], [(88999,"hyRin"),(13604,"W")]),
 79 => ([(1,"shB"),(7,"stv")], [(10,"pv"),(2,"W")]), 80 => ([(97763,"hyRin")], [(5377,"shB"),(80654,"pv"),(52303,"per")]),
 82 => ([(300,"phD"),(1226,"br")], [(158,"shB"),(1249,"W")]), 83 => ([(1017,"shB"),(7839,"stv")], [(1800,"phD"),(8046,"pv")]),
 84 => ([(89400,"hyRin")], [(88506,"pv"),(82248,"per"),(9834,"W")]), 85 => ([(1,"shB")], [(3,"pv"),(7,"per"),(2,"W")]),
 86 => ([(4145,"shB"),(72538,"phD")], [(127045,"pv"),(90258,"W")]), 87 => ([(11947,"shB")], [(2432,"pv"),(92087,"per"),(21145,"phD")]),
 88 => ([(200,"phD")], [(236,"pv"),(80,"stv"),(226,"W")]), 89 => ([(100,"phD"),(40,"per")], [(158,"pv"),(113,"W")]),
 30 => ([(1,"wd")], [(1,"fo")]), 64 => ([(1,"rg")], [(1,"wd")]), 77 => ([(1,"rg")], [(1,"pv"),(1,"per")]), 75 => ([(1,"aki")], [(1,"pv")]),
 66 => ([(1,"rg"),(1,"stv")], [(2,"aki")]), 65 => ([(1,"wd"),(1,"stv")], [(2,"aki")]), 57 => ([(1,"wd"),(1,"stv")], [(2,"hen")]),
 81 => ([(1,"aki"),(1,"per")], [(1,"rg")]))
dGr(k, P, T, lna) = (s = 0.0; for (c, n) in RX[k][2]; s += c*(Gr(n, P, T) + (n == "W" ? Rgas*T*lna : 0.0)); end;
                     for (c, n) in RX[k][1]; s -= c*(Gr(n, P, T) + (n == "W" ? Rgas*T*lna : 0.0)); end; s/sum(first.(RX[k][1])))
root(g, x0, x1) = (g0 = g(x0); g(x1)*g0 > 0 && return NaN; for _ in 1:80; m = (x0 + x1)/2; gm = g(m); gm*g0 > 0 ? (x0 = m; g0 = gm) : (x1 = m); end; (x0 + x1)/2)
chk = [(87, 30.33, 700, NaN, "Ohtani"), (87, 30.90, 1000, NaN, "Ohtani"), (83, 26.97, 916, NaN, "Fig 5"), (55, 16.85, 1097, NaN, "Fig 5"),
       (4, 7.53, 770, NaN, "Fig 5"), (80, 24.56, 1060, NaN, "Fig 5"), (68, 21.15, 810, NaN, "Fig 5"), (30, 14.57, 1488, NaN, "dry"), (30, 12.90, 520, NaN, "dry"),
       (77, 24.13, 1453, NaN, "dry"), (77, 24.20, 937, NaN, "dry"), (75, 22.73, 1427, NaN, "dry"), (66, 20.05, 1067, NaN, "dry"), (65, 19.21, 1309, NaN, "dry"),
       (57, 16.98, 1446, NaN, "dry"), (64, 18.30, 1076, NaN, "dry"), (81, 23.64, 814, NaN, "dry"), (3, 5.19, 541, 1.0, "J1"), (3, 6.08, 623, 1.0, "Fig 5"),
       (3, 8.32, 792, 1.0, "J2"), (5, 8.32, 792, 1.0, "J2"), (6, 8.32, 792, 1.0, "J2"), (5, 12.27, 1064, 0.68, "J3"), (13, 12.27, 1064, 0.68, "J3"),
       (14, 12.27, 1064, 0.68, "J3"), (14, 13.53, 1208, 0.60, "J4"), (43, 13.53, 1208, 0.60, "J4"), (53, 17.36, 1128, 0.54, "Jb"), (62, 17.36, 1128, 0.54, "Jb"),
       (22, 15.37, 877, 0.85, "J6"), (42, 15.37, 877, 0.85, "J6"), (51, 15.37, 877, 0.85, "J6"), (11, 11.37, 586, 1.0, "J7"), (22, 11.37, 586, 1.0, "J7"),
       (42, 14.49, 546, 1.0, "J8"), (82, 14.49, 546, 1.0, "J8"), (82, 22.41, 732, 1.0, "Fig 5"), (78, 23.33, 1251, 0.08, "J9"), (79, 23.33, 1251, 0.08, "J9"),
       (79, 25.44, 1188, 0.10, "J14"), (86, 25.44, 1188, 0.10, "J14"), (88, 25.44, 1188, 0.10, "J14"), (78, 25.88, 1403, 0.05, "J10"), (84, 25.88, 1403, 0.05, "J10"),
       (85, 25.88, 1403, 0.05, "J10"), (76, 25.88, 1403, 0.05, "J10"), (85, 31.02, 1040, 0.03, "J11"), (86, 31.02, 1040, 0.03, "J11"), (89, 31.02, 1040, 0.03, "J11"),
       (84, 24.33, 1601, 0.012, "J12")]
rows = map(chk) do (k, P0, Tc, a, src)
    T0 = Tc + 273.15; la = isnan(a) ? 0.0 : log(a); lm(P, T) = isnan(a) ? 0.0 : ln_ah2o(P, T, A)
    (k=k, P=P0, T=Tc, a=a, src=src, dP=root(P -> dGr(k, P, T0, la), max(P0 - 6, 1.1), P0 + 6) - P0, dT=root(T -> dGr(k, P0, T, la), T0 - 400, T0 + 400) - T0,
     dPm=root(P -> dGr(k, P, T0, lm(P, T0)), max(P0 - 6, 1.1), P0 + 6) - P0, dTm=root(T -> dGr(k, P0, T, lm(P0, T)), T0 - 400, T0 + 400) - T0)
end
ok2(x) = abs(x.dP) <= 0.2 || abs(x.dT) <= 20; ok2m(x) = abs(x.dPm) <= 0.2 || abs(x.dTm) <= 20
dh = [x for x in rows if !isnan(x.a)]
@printf("S06: dehydration within 0.2 GPa or 20 K: %d/%d at K&O's a, %d/%d with our map; fluid-absent median |dP| %.2f GPa, dry %.2f\n",
        count(ok2, dh), length(dh), count(ok2m, dh), length(dh), median([abs(x.dP) for x in rows if isnan(x.a) && x.src != "dry"]),
        median([abs(x.dP) for x in rows if x.src == "dry"]))
Pg6 = range(5, 35, length=241); Tc6 = range(500, 1600, length=161); W6 = zeros(length(Pg6), length(Tc6)); as6 = fill("", length(Pg6), length(Tc6))
function msh_map!(W, asm, Pv, Tv, MgSi)
    Threads.@threads for i in eachindex(Pv)
        for j in eachindex(Tv)
            P = Pv[i]; T = Tv[j] + 273.15; mu = mu_h2o(P, T) + Rgas*T*ln_ah2o(P, T, A)
            G = ntuple(k -> (k == 4 && P > 50.0) ? 1e12 : g_ko(K[k], P, T) - K[k].h*mu, Val(17)); k1, k2, n1, n2 = msh_min(MgSi, 1.0, G, K)
            w = 18.015(n1*K[k1].h + n2*K[k2].h); m = sum(n*(K[q].m*40.304 + K[q].s*60.084) for (q, n) in ((k1, n1), (k2, n2))) + w
            W[i, j] = 100w/m; asm[i, j] = join(sort([KO_NAME[q] for (q, n) in ((k1, n1), (k2, n2)) if n > 1e-9]), "+")
        end
    end
end
msh_map!(W6, as6, Pg6, Tc6, 1.4)
f = Figure(size=(1400, 1050))
a = Axis(f[1, 1:2], title="a  MgO-SiO2-H2O at Mg/Si = 1.4, water-saturated at the fitted a(H2O), our minimiser", xlabel="Pressure (GPa)", ylabel="Temperature (°C)")
hm = heatmap!(a, Pg6, Tc6, W6, colormap=seq, colorrange=(0, 12), interpolate=true); boundaries!(a, Pg6, Tc6, as6)
labels!(a, Pg6, Tc6, as6; nmin=60)
scatter!(a, [x.P for x in rows if isnan(x.a)], [x.T for x in rows if isnan(x.a)], marker=:diamond, markersize=10, color=:white, strokecolor=ink, strokewidth=1.2,
         label="K&O fluid-absent / dry checkpoint")
scatter!(a, [x.P for x in dh], [x.T for x in dh], markersize=10, color=:white, strokecolor=cat3[2], strokewidth=1.8, label="K&O dehydration checkpoint")
axislegend(a, position=:rt, labelsize=10, framevisible=true)
Colorbar(f[1, 3], hm, label="H2O in the solid (wt%)", width=12)
grp(x) = isnan(x.a) ? (x.src == "dry" ? 2 : 1) : 3
labs = ["R$(x.k) $(round(x.P, digits=1))/$(round(Int, x.T))" for x in rows]
for (col, (fld, fldm, yl)) in enumerate(((:dP, :dPm, "ΔP at the checkpoint T (GPa)"), (:dT, :dTm, "ΔT at the checkpoint P (K)")))
    ax = Axis(f[2, col], ylabel=yl, xticks=(1:length(rows), labs), xticklabelrotation=π/2, xticklabelsize=8, title=col == 1 ? "b  Our boundary minus K&O's" : "")
    hlines!(ax, [0.0], color=muted, linewidth=0.8)
    for g in 1:3
        id = [k for k in eachindex(rows) if grp(rows[k]) == g]
        scatter!(ax, id, [getfield(rows[k], fld) for k in id], color=cat3[g], markersize=9,
                 label=("fluid-absent (incl. Ohtani calibration)", "dry (Fabrichnaya 2004)", "dehydration at K&O's a(H2O)")[g])
    end
    id = [k for k in eachindex(rows) if grp(rows[k]) == 3]
    scatter!(ax, id, [getfield(rows[k], fldm) for k in id], color=:white, strokecolor=cat3[3], strokewidth=1.5, markersize=9, label="dehydration at our a(H2O) map")
    col == 1 && axislegend(ax, position=:lb, labelsize=9, framevisible=true)
end
save_fig("S06_ko_benchmark.png", f)

# ---- S07: the a(H2O) map of the DHMS minimisation
Pa7 = range(5, 64, length=296); Ta7 = range(400, 1700, length=131); La = [log10(exp(ln_ah2o(P, T + 273.15, A))) for P in Pa7, T in Ta7]
f = Figure(size=(1100, 560)); a = Axis(f[1, 1], xlabel="Pressure (GPa)", ylabel="Temperature (°C)", title="a(H2O): thin-plate spline in ln a through the labelled anchors")
hm = heatmap!(a, Pa7, Ta7, La, colormap=Reverse(seq), colorrange=(-2, 0))
contour!(a, Pa7, Ta7, La, levels=log10.([0.02, 0.05, 0.1, 0.2, 0.5, 0.8]), color=ink, linewidth=0.7, labels=true, labelsize=9)
poly!(a, Rect(5, 500, 55, 1100), color=:transparent, strokecolor=muted, strokewidth=1, linestyle=:dash)
nk = size(KO_AH2O, 1); Q = A.Q
scatter!(a, Q[1:nk, 1], Q[1:nk, 2], color=:white, strokecolor=ink, strokewidth=1.2, markersize=9, label="K&O 2006 fits, a = 1 from Komabayashi 2005")
scatter!(a, Q[nk+1:end, 1], Q[nk+1:end, 2], color=:white, strokecolor=cat3[2], strokewidth=1.8, markersize=9, label="Nishi 2014 dehydration lines (phase D, phase H)")
for k in axes(Q, 1); text!(a, Q[k, 1], Q[k, 2], text=string(round(Q[k, 3], sigdigits=2)), offset=(5, 3), fontsize=9, color=ink); end
axislegend(a, position=:lb, labelsize=9, framevisible=true); Colorbar(f[1, 2], hm, label="log10 a(H2O)", width=12)
save_fig("S07_activity.png", f)

# ---- S08: phase H in Nishi et al. (2014) Fig 2 bulks, with their boundaries and observations
nishiA = (([45.4, 48.6], [600.0, 1035.0], "phD = phH + St"), ([30.0, 40.0, 48.6], [1187.0, 1105.0, 1035.0], "phD dehydration"),
          ([48.6, 50.0, 57.9], [1035.0, 1030.0, 1056.0], "phH + St dehydration"))
nishiB = (([28.0, 31.5, 34.0, 39.5], [630.0, 800.0, 900.0, 1000.0], "phH in"), ([40.0, 50.0, 57.9], [997.0, 1030.0, 1056.0], "phH dehydration"))
obsH = [(34.4, 900), (34.9, 900), (40.5, 1000), (41.0, 800), (41.0, 900), (41.0, 700), (42.0, 1000), (45.0, 1000), (50.0, 1000)]
f = Figure(size=(1400, 560))
for (col, (MgSi, ttl, ln_)) in enumerate(((1/1.5, "a  Phase D bulk (Nishi A, Mg/Si = 0.67)", nishiA), (1.0, "b  Phase H bulk (Nishi B, MgSiO4H2)", nishiB)))
    Pp = range(20, 62, length=211); Tp = range(500, 1400, length=181); W8 = zeros(length(Pp), length(Tp)); as8 = fill("", length(Pp), length(Tp))
    msh_map!(W8, as8, Pp, Tp, MgSi)
    ax = Axis(f[1, col], title=ttl, xlabel="Pressure (GPa)", ylabel=col == 1 ? "Temperature (°C)" : "")
    h8 = heatmap!(ax, Pp, Tp, W8, colormap=seq, colorrange=(0, 15), interpolate=true); boundaries!(ax, Pp, Tp, as8)
    labels!(ax, Pp, Tp, as8; nmin=80)
    for (n, (x, y, l)) in enumerate(ln_); scatterlines!(ax, x, y, color=cat3[n], linewidth=2.5, linestyle=:dash, markersize=9, label="Nishi: " * l); end
    if col == 2
        scatter!(ax, first.(obsH), last.(obsH), marker=:rect, color=:white, strokecolor=ink, strokewidth=1.2, markersize=9, label="Nishi: phase H observed")
        scatter!(ax, [50.0], [1650 - 273.15], marker=:utriangle, color=:white, strokecolor=ink, strokewidth=1.2, markersize=11, label="Ohtani 2014: phase H out (DAC)")
        Colorbar(f[1, 3], h8, label="H2O in the solid (wt%)", width=12)
    end
    axislegend(ax, position=:lt, labelsize=9, framevisible=true)
end
save_fig("S08_phaseH.png", f)

# ---- S09: the exact MSH minimisation as a lower convex hull of G* per mole of oxide
f = Figure(size=(1400, 520))
for (col, (P, T)) in enumerate(((20.0, 1100.0), (40.0, 1100.0)))
    mu = mu_h2o(P, T) + Rgas*T*ln_ah2o(P, T, A)
    G = [(k == 4 && P > 50.0) ? NaN : (g_ko(K[k], P, T) - K[k].h*mu)/(K[k].m + K[k].s) for k in eachindex(K)]; x = [K[k].s/(K[k].m + K[k].s) for k in eachindex(K)]
    ref(x) = (1 - x)*G[15] + x*G[16]                                                   # per + stv chord as the zero line
    a = Axis(f[1, col], xlabel="SiO2/(MgO+SiO2)", ylabel="G* − chord(per, stv) (kJ per mol oxide)", title=@sprintf("%s  %.0f GPa, %.0f K", ("a", "b")[col], P, T))
    y = [(G[k] - ref(x[k]))/1e3 for k in eachindex(K)]
    for k in eachindex(K)
        isnan(y[k]) && continue
        scatter!(a, [x[k]], [y[k]], color=KO_DHMS[k] ? cat3[2] : cat3[1], markersize=10)
        text!(a, x[k], y[k], text=KO_NAME[k], offset=(4, 4), fontsize=9, color=ink)
    end
    hull = Int[]; ord = sortperm(x); ord = [k for k in ord if !isnan(y[k])]
    for k in ord
        while length(hull) >= 2
            i, j = hull[end-1], hull[end]; ((x[j] - x[i])*(y[k] - y[i]) - (y[j] - y[i])*(x[k] - x[i]) <= 0) ? pop!(hull) : break
        end
        push!(hull, k)
    end
    lines!(a, x[hull], y[hull], color=ink, linewidth=1.5)
    for (MgSi, lab) in ((1.4, "harzburgite-like 1.4"), (0.5, "Si-rich 0.5")); vlines!(a, [1/(1 + MgSi)], color=muted, linestyle=:dash); text!(a, 1/(1 + MgSi), 0, text=lab, rotation=π/2, fontsize=9, color=muted, align=(:left, :bottom)); end
end
save_fig("S09_msh_hull.png", f)

# ---- S10 / S11 / S12: capacity with host, DHMS gain, share held by hydrous phases (all sectors)
for (fn, M, cm, cr, cl, lab) in (("S10_capacity.png", r.S, seq, (-3, 1.3), "log10 capacity (wt% of solid); dashed: 0.01, 0.1, 1 wt%", true),
                                 ("S11_dhms_gain.png", r.S ./ r0.S, div, (-2.5, 2.5), "log10 (with DHMS / without)", false),
                                 ("S12_hydrous_share.png", r.Shyd ./ r.S, seq, (0, 1), "share of the capacity in hydrous phases", false))
    f = Figure(size=(1400, 1150)); hms = Any[]
    for (row, c) in enumerate((1, 3, 2)), q in 1:3
        P, T = r.sectors[q]; C = cm === div || lab ? log10.(max.(M[:, :, q, c], 1e-4)) : M[:, :, q, c]
        a = Axis(f[row, q], title="$(bnames[c]): $(snames[q])", xlabel=row == 3 ? "Pressure (GPa)" : "", ylabel=q == 1 ? "Temperature (K)" : "")
        push!(hms, heatmap!(a, collect(P), collect(T), C, colormap=cm, colorrange=cr, interpolate=true))
        if lab
            contour!(a, collect(P), collect(T), smooth(C), levels=[-2, -1, 0], color=ink, linewidth=0.7, linestyle=:dash)
            H = copy(r.Shost[:, :, q, c]); boundaries!(a, P, T, H; color=(ink, 0.6))
            labels!(a, P, T, H)
        end
        q == 3 && vlines!(a, [50.0, 60.0], color=muted, linestyle=:dash, linewidth=1)
        q == 1 && c != 2 && vlines!(a, [5.0, 6.0], color=muted, linestyle=:dot, linewidth=1)
        q == 2 && c == 2 && vlines!(a, [9.5, 10.0], color=muted, linestyle=:dot, linewidth=1)
    end
    Colorbar(f[1:3, 4], hms[1], label=cl, width=12)
    save_fig(fn, f)
end

# ---- S13: capacity along illustrative geotherms, with and without DHMS
f = Figure(size=(1400, 820)); Pgeo = collect(range(0.5, 134.5, length=540))
for (col, (nm, g)) in enumerate(geo)
    a = Axis(f[1, col], yscale=log10, title=nm, xlabel="Pressure (GPa)", ylabel=col == 1 ? "Capacity (wt% of solid)" : "")
    b = Axis(f[2, col], xlabel="Pressure (GPa)", ylabel=col == 1 ? "Temperature (K)" : "", title="the geotherm")
    for c in (1, 3, 2)
        lines!(a, Pgeo, [max(at(r.S, sector(P), c, P, g(P)), 1e-4) for P in Pgeo], color=cat3[c], linewidth=2.2, label=bnames[c])
        lines!(a, Pgeo, [max(at(r0.S, sector(P), c, P, g(P)), 1e-4) for P in Pgeo], color=cat3[c], linewidth=1, linestyle=:dash)
    end
    vlines!(a, [7.0, 25.0], color=muted, linewidth=0.8); ylims!(a, 1e-3, 30)
    lines!(b, Pgeo, g.(Pgeo), color=ink, linewidth=2); vlines!(b, [7.0, 25.0], color=muted, linewidth=0.8)
    col == 1 && axislegend(a, position=:rt, labelsize=9, framevisible=true)
end
save_fig("S13_geotherms.png", f)

# ---- S14: capacity against T at fixed P, with (solid) and without (dashed) DHMS
Ps14 = (6.5, 10.0, 13.0, 16.0, 19.0, 22.0, 28.0, 40.0, 55.0); pc = cgrad([colorant"#9ec5f4", colorant"#0d366b"], length(Ps14), categorical=true)
f = Figure(size=(1400, 480))
for (col, c) in enumerate((1, 3, 2))
    a = Axis(f[1, col], yscale=log10, title=bnames[c], xlabel="Temperature (K)", ylabel=col == 1 ? "Capacity (wt% of solid)" : "")
    for (n, Pv) in enumerate(Ps14)
        q = sector(Pv); P, T = r.sectors[q]; i = argmin(abs.(collect(P) .- Pv))
        lines!(a, collect(T), max.(r.S[i, :, q, c], 1e-4), color=pc[n], linewidth=2, label=@sprintf("%.1f GPa", P[i]))
        lines!(a, collect(T), max.(r0.S[i, :, q, c], 1e-4), color=pc[n], linewidth=1, linestyle=:dash)
    end
    xlims!(a, 500, 2200); ylims!(a, 1e-3, 30); col == 3 && axislegend(a, position=:rt, labelsize=9, nbanks=2, framevisible=true)
end
save_fig("S14_T_profiles.png", f)

# ---- S15: seams and hand-overs, capacity against P through 0-30 GPa at fixed T
f = Figure(size=(1400, 480)); Tsm = (800.0, 1000.0, 1300.0, 1700.0)
for (col, c) in enumerate((1, 3, 2))
    a = Axis(f[1, col], yscale=log10, title=bnames[c], xlabel="Pressure (GPa)", ylabel=col == 1 ? "Capacity (wt% of solid)" : "")
    vspan!(a, (c == 2 ? (9.5, 10.0) : (5.0, 6.0))..., color=(faint, 0.8))
    for (n, T) in enumerate(Tsm)
        for q in 1:3
            P, Ta = r.sectors[q]; (T < Ta[1] || T > Ta[end]) && continue
            j = argmin(abs.(collect(Ta) .- T)); sel = collect(P) .<= 30
            lines!(a, collect(P)[sel], max.(r.S[sel, j, q, c], 1e-4), color=ramp4[n], linewidth=2, label=q == 1 || (q == 2 && T < 1000) ? "$(Int(T)) K" : nothing)
        end
    end
    vlines!(a, [7.0, 25.0], color=muted, linewidth=0.8)
    xlims!(a, 0, 30); ylims!(a, 1e-3, 30); col == 1 && axislegend(a, position=:rt, labelsize=9, framevisible=true, unique=true)
end
save_fig("S15_seams.png", f)

# ---- S16: fO2 and native iron (StagfO2 / StagMetal)
fo = read_fo2(joinpath(run, "StagfO2.dat")); me = read_fo2(joinpath(run, "StagMetal.dat"))
ir(R) = argmin(abs.(fo.R .- R)); Rp = (ir(-0.5), ir(0.0), ir(0.03), ir(0.3))
f = Figure(size=(1400, 1150)); hms = Any[]
for (row, k) in enumerate((1, 2)), q in 1:3
    P = LinRange(fo.sec[q][1], fo.sec[q][2], size(fo.X, 1)); T = LinRange(fo.sec[q][3], fo.sec[q][4], size(fo.X, 2))
    a = Axis(f[row, q], title="$(bnames[k]): log10 fO2 at R = $(round(fo.R[Rp[3]], digits=3)), $(snames[q])", xlabel=row == 2 ? "Pressure (GPa)" : "",
             ylabel=q == 1 ? "Temperature (K)" : "")
    push!(hms, heatmap!(a, collect(P), collect(T), fo.X[:, :, Rp[3], k, q], colormap=seq, colorrange=(-25, 5)))
    contour!(a, collect(P), collect(T), fo.X[:, :, Rp[3], k, q], levels=-24:4:4, color=ink, linewidth=0.6, labels=true, labelsize=8)
end
Colorbar(f[1:2, 4], hms[1], label="log10 fO2", width=12)
for (col, (nm, q, P0, T0)) in enumerate((("UM 3 GPa, 1600 K", 1, 3.0, 1600.0), ("TZ 18 GPa, 1900 K", 2, 18.0, 1900.0), ("LM 60 GPa, 2300 K", 3, 60.0, 2300.0)))
    Pq = LinRange(fo.sec[q][1], fo.sec[q][2], size(fo.X, 1)); Tq = LinRange(fo.sec[q][3], fo.sec[q][4], size(fo.X, 2))
    i = argmin(abs.(Pq .- P0)); j = argmin(abs.(Tq .- T0))
    a = Axis(f[3, col], title="fO2 and metal against R, $nm", xlabel="R = Fe3+/ΣFe (signed; < 0: metal)", ylabel=col == 1 ? "log10 fO2" : "")
    b = Axis(f[3, col], yaxisposition=:right, ylabel=col == 3 ? "native Fe (wt fraction, dashed)" : ""); hidexdecorations!(b); hidespines!(b)
    for (m, k) in enumerate((1, 2, 3, 6))
        lines!(a, fo.R, fo.X[i, j, :, k, q], color=cat4[m], linewidth=2, label=k <= 2 ? bnames[k] : @sprintf("mix %.2f", fo.mix[k-2]))
        lines!(b, fo.R, me.X[i, j, :, k, q], color=cat4[m], linewidth=1, linestyle=:dash)
    end
    col == 1 && axislegend(a, position=:lt, labelsize=9, framevisible=true)
end
save_fig("S16_fO2_metal.png", f)

# ---- S17: silica-saturation step (StagfO2step): basalt fraction at olivine-out (pa) and silica-in (pb)
st = read_step(joinpath(run, "StagfO2step.dat")); kR = argmin(abs.(st.R .- 0.01))
f = Figure(size=(1400, 820)); hms = Any[]
for (row, (s, nm)) in enumerate(((1, "pa: last basalt fraction with olivine-family or mw"), (2, "pb: first basalt fraction with free silica"))), q in 1:3
    P = LinRange(st.sec[q][1], st.sec[q][2], size(st.X, 1)); T = LinRange(st.sec[q][3], st.sec[q][4], size(st.X, 2))
    a = Axis(f[row, q], title="$(snames[q]): $nm", xlabel=row == 2 ? "Pressure (GPa)" : "", ylabel=q == 1 ? "Temperature (K)" : "", titlesize=10)
    push!(hms, heatmap!(a, collect(P), collect(T), st.X[:, :, kR, q, s], colormap=seq, colorrange=(0.3, 1.0), nan_color=faint))
end
Colorbar(f[1:2, 4], hms[1], label="basalt fraction (R = $(st.R[kR]))", width=12)
save_fig("S17_step.png", f)

# ---- S19: capacity against basalt fraction at single nodes: the ten table compositions and the silica step
f = Figure(size=(1400, 820)); o = sortperm(rf.pX)
for (n, (nm, q, P0, T0)) in enumerate((("UM 3 GPa, 800 K", 1, 3.0, 800.0), ("UM 3 GPa, 1500 K", 1, 3.0, 1500.0), ("TZ 12 GPa, 1000 K", 2, 12.0, 1000.0),
                                       ("TZ 18 GPa, 1700 K", 2, 18.0, 1700.0), ("LM 45 GPa, 1200 K", 3, 45.0, 1200.0), ("LM 90 GPa, 2500 K", 3, 90.0, 2500.0)))
    P, T = rf.sectors[q]; i = argmin(abs.(collect(P) .- P0)); j = argmin(abs.(collect(T) .- T0))
    a = Axis(f[cld(n, 3), mod1(n, 3)], yscale=log10, title=@sprintf("%s (node %.1f GPa, %.0f K)", nm, P[i], T[j]), xlabel=n > 3 ? "basalt fraction p" : "",
             ylabel=mod1(n, 3) == 1 ? "Capacity (wt% of solid)" : "")
    Ps = LinRange(st.sec[q][1], st.sec[q][2], size(st.X, 1)); Ts = LinRange(st.sec[q][3], st.sec[q][4], size(st.X, 2))
    pa, pb = st.X[argmin(abs.(Ps .- P0)), argmin(abs.(Ts .- T0)), kR, q, :]
    isfinite(pa) && isfinite(pb) && vspan!(a, pa, pb, color=(faint, 0.9))
    lines!(a, rf.pX[o], max.(rf.S[i, j, q, o], 1e-4), color=cat3[1], linewidth=2, label="with DHMS")
    scatter!(a, rf.pX[o], max.(rf.S[i, j, q, o], 1e-4), color=cat3[1], markersize=9)
    lines!(a, rf.pX[o], max.(rf0.S[i, j, q, o], 1e-4), color=cat3[2], linewidth=1.2, linestyle=:dash, label="without")
    xlims!(a, 0, 1); n == 1 && axislegend(a, position=:rt, labelsize=9, framevisible=true)
end
save_fig("S19_composition.png", f)

# ---- S18: hydrous density corrections (StagSdens / StagMdens) and the melt water volume
sd = read_dens(joinpath(run, "StagSdens.dat")); md = read_dens(joinpath(run, "StagMdens.dat")); ih = argmin(abs.(sd.H .- 4.0))
f = Figure(size=(1400, 820)); hms = Any[]
for (col, (D, nm, s)) in enumerate(((sd, "solid, harzburgite", 2), (sd, "solid, basalt", 1), (md, "melt, harzburgite", 2)))
    a = Axis(f[1, col], title=@sprintf("%s: %% density drop at %.1f wt%% H2O", nm, D.H[ih]), xlabel="Pressure (GPa)", ylabel=col == 1 ? "Temperature (K)" : "")
    push!(hms, heatmap!(a, collect(D.P), collect(D.T), D.X[:, :, ih, s], colormap=seq, colorrange=(0, 15)))
end
Colorbar(f[1, 4], hms[1], label="% density decrease against dry", width=12)
Xl = ["SiO2", "MgO", "FeO", "CaO", "Al2O3", "Na2O", "Cr2O3", "O", "H2O"]
Xh = filter_ox_list([43.43, 45.93, 8.34, 0.9, 1.0, 0.01, 0.3, 0.0, 0.0], Xl, Xl, wt_in=true, wt_out=false); Xh ./= sum(Xh)
MMv = get_Xoxmm(Xl); MA = sum(Xh .* MMv); Pv = range(1e-4, 50, length=120); w = 4.0
x = (w/MMv[end])/((w/MMv[end]) + (100 - w)/MA); Xw = (1 - x) .* Xh; Xw[end] = x; MH = (1 - x)*MA + x*MMv[end]
a = Axis(f[2, 1:2], xlabel="Pressure (GPa)", ylabel="partial molar V of H2O in melt (cm³/mol)", title="melt water volume: ideal (used) and Jing & Karato hard sphere, 4 wt%, harzburgite melt")
for (n, T) in enumerate((1800.0, 2400.0, 3000.0))
    V1 = (MH/Jing_Karato_silicate_melt_density(1e-4, T, Xw, Xl) - (1 - x)*MA/Jing_Karato_silicate_melt_density(1e-4, T, Xh, Xl))/x
    lines!(a, Pv, [melt_H2O_partial_volume(P, V1) for P in Pv], color=ramp4[n+1], linewidth=2, label="ideal, $(Int(T)) K")
    lines!(a, Pv, [(MH/Jing_Karato_silicate_melt_density(P, T, Xw, Xl) - (1 - x)*MA/Jing_Karato_silicate_melt_density(P, T, Xh, Xl))/x for P in Pv],
           color=ramp4[n+1], linewidth=1.2, linestyle=:dash, label="JK, $(Int(T)) K")
end
hlines!(a, [0.0], color=muted, linewidth=0.8); axislegend(a, position=:rt, labelsize=9, nbanks=2, framevisible=true)
b = Axis(f[2, 3], xlabel="H2O (wt%)", ylabel="% density decrease", title="all three at 20 GPa, 2000 K")
i20 = argmin(abs.(sd.P .- 20)); j20 = argmin(abs.(sd.T .- 2000))
for (m, (D, s, nm)) in enumerate(((sd, 2, "solid harz"), (sd, 1, "solid basalt"), (md, 2, "melt harz"), (md, 1, "melt basalt")))
    lines!(b, collect(D.H), D.X[i20, j20, :, s], color=cat4[m], linewidth=2, label=nm)
end
axislegend(b, position=:lt, labelsize=9, framevisible=true)
save_fig("S18_density.png", f)
println("all figures written to ", figs)
