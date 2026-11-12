# Worked examples of the three fluid EoS on the ideality pages of TOex_solver_flow.tex (pages 6a-6c)
using Printf
src = read(joinpath(homedir(),"Desktop/StagPP.jl/src/TOex.jl"), String)
a = findfirst("# --- Zhang & Duan (2009) equation of state", src).start
b = findfirst("# ---- Section 3b", src).start
module E; const R = 8.31446261815324; end
Base.include_string(E, src[a:b-1])
using .E
x = zeros(7); x[2] = 0.9; x[3] = 0.1               # H2O 0.9, CO2 0.1
println("=== ZD09  X_H2O 0.9 X_CO2 0.1, 3 GPa, 1600 K")
P, T = 3.0e4, 1600.0
ε=0.0; σ=0.0; εi=zeros(7); σi=zeros(7)
for i in 1:7, j in 1:7
    global ε += x[i]*x[j]*E.ZD_K1[i,j]*√(E.ZD_eps[i]*E.ZD_eps[j]); global σ += x[i]*x[j]*E.ZD_K2[i,j]*(E.ZD_sig[i]+E.ZD_sig[j])/2
    εi[i] += E.ZD_K1[i,j]*x[j]*√(E.ZD_eps[i]*E.ZD_eps[j]); σi[i] += E.ZD_K2[i,j]*x[j]*(E.ZD_sig[i]+E.ZD_sig[j])/2
end
Tm = 154*T/ε; Pm = 3.06258806161*σ^3*P/ε; Vm = E.zd_findroot(Tm, Pm)
c = E.ZD_c; g = [c[1]+c[2]/Tm^2+c[3]/Tm^3, c[4]+c[5]/Tm^2+c[6]/Tm^3, c[7]+c[8]/Tm^2+c[9]/Tm^3, c[10]+c[11]/Tm^2+c[12]/Tm^3]
Z = 1 + g[1]/Vm + g[2]/Vm^2 + g[3]/Vm^4 + g[4]/Vm^5 + c[13]/Tm^3/Vm^2*(c[14]+c[15]/Vm^2)*exp(-c[15]/Vm^2)
V = 1000Vm*(σ/3.691)^3
@printf("eps %.2f K  sigma %.4f A  Tm %.3f  Pm %.4e  Vm %.5f L/mol  Z %.4f  V %.3f cm3/mol\n", ε, σ, Tm, Pm, Vm, Z, V)
@printf("terms of Z-1: g1/V %.4f  g2/V2 %.4f  g3/V4 %.4f  g4/V5 %.4f  exp %.4f\n", g[1]/Vm, g[2]/Vm^2, g[3]/Vm^4, g[4]/Vm^5, c[13]/Tm^3/Vm^2*(c[14]+c[15]/Vm^2)*exp(-c[15]/Vm^2))
vol, fug = E.zd_fugacity(x, P, T)
S1 = (c[1]+c[2]/Tm^2+c[3]/Tm^3)/Vm + (c[4]+c[5]/Tm^2+c[6]/Tm^3)/(2Vm^2) + (c[7]+c[8]/Tm^2+c[9]/Tm^3)/(4Vm^4) + (c[10]+c[11]/Tm^2+c[12]/Tm^3)/(5Vm^5) + c[13]/(2c[15]*Tm^3)*(c[14]+1-(c[14]+1+c[15]/Vm^2)*exp(-c[15]/Vm^2))
S2 = (2c[2]/Tm^2+3c[3]/Tm^3)/Vm + (2c[5]/Tm^2+3c[6]/Tm^3)/(2Vm^2) + (2c[8]/Tm^2+3c[9]/Tm^3)/(4Vm^4) + (2c[11]/Tm^2+3c[12]/Tm^3)/(5Vm^5) + 3c[13]/(2c[15]*Tm^3)*(c[14]+1-(c[14]+1-c[15]/Vm^2)*exp(-c[15]/Vm^2))
@printf("Z-1-lnZ %.4f  S1 %.4f  S2 %.4f\n", Z-1-log(Z), S1, S2)
for (i,nm) in ((2,"H2O"),(3,"CO2"))
    @printf("%s: eps_i'/eps %.4f sigma_i'/sigma %.4f  -2S2(1-e) %.4f  6(1-Z)(1-s) %.4f  ln phi %.4f  phi %.3f\n", nm, εi[i]/ε, σi[i]/σ, -2S2*(1-εi[i]/ε), 6(1-Z)*(1-σi[i]/σ), log(fug[i]/(x[i]*P)), fug[i]/(x[i]*P))
end
@printf("check V from zd_fugacity %.3f\n", vol)

println("=== MRK  same composition, 0.3 GPa, 600 K")
P, T = 3000.0, 600.0
av, bv = E.mrk_ab(T); awc = E.mrk_awc(T)
@printf("a_H2O %.4e  a_CO2 %.4e  sqrt(a_w a_c) %.4e  a_wc %.4e   b_H2O %.2f b_CO2 %.2f\n", av[2], av[3], sqrt(av[2]*av[3]), awc, bv[2], bv[3])
am = x[2]^2*av[2] + x[3]^2*av[3] + x[2]*x[3]*awc; bm = x[2]*bv[2] + x[3]*bv[3]
rt = E.RK_R*T; st = √T; c1 = -rt/P
r = E.cubic_roots(c1, c1*bm + am/(st*P) - bm^2, -am*bm/(st*P))
@printf("a_mix %.4e  b_mix %.3f  RT/P %.3f  roots %s\n", am, bm, rt/P, string(round.(r, digits=3)))
v, lnφ = E.mrk_lnφ(x, P, T)
d2 = log((v+bm)/v)/(rt*st*bm); d3 = am*d2/bm - am/((v+bm)*rt*st*bm) + 1/(v-bm); d6 = log(rt/(v-bm)) - log(P)
@printf("chosen V %.3f cm3/mol  ln(RT/(P(V-b))) %.4f\n", v, d6)
for (i,nm) in ((2,"H2O"),(3,"CO2"))
    s = sum(x[j]*(i in (2,3) && j in (2,3) && i!=j ? awc/2 : √(av[i]*av[j])) for j in 1:7)
    @printf("%s: b_i d3 %.4f  -2 sum x_j a_ij d2 %.4f  ln phi_MRK %.4f  phi %.4f\n", nm, bv[i]*d3, -2s*d2, lnφ[i], exp(lnφ[i]))
end
println("=== PS94 pure H2O and CO2, 0.3 GPa, 600 K; hybrid correction")
for (i,k,nm) in ((2,1,"H2O"),(3,2,"CO2"))
    vps, lps = E.ps94_pure(P, T, k)
    e = zeros(7); e[i] = 1.0; vm, lm = E.mrk_lnφ(e, P, T)
    @printf("%s pure: PS94 V %.3f ln phi %.4f | MRK pure V %.3f ln phi %.4f | delta %.4f\n", nm, vps, lps, vm, lm[i], lps - lm[i])
end
vh, fh = E.hyb_fugacity(x, P, T)
@printf("hybrid mixture: V %.3f  ln phi H2O %.4f  CO2 %.4f\n", vh, log(fh[2]/(x[2]*P)), log(fh[3]/(x[3]*P)))
vz, fz = E.zd_fugacity(x, P, T)
@printf("ZD09 at the same point: V %.3f  ln phi H2O %.4f CO2 %.4f\n", vz, log(fz[2]/(x[2]*P)), log(fz[3]/(x[3]*P)))
println("=== blend at 723 K, 0.3 GPa")
T = 723.0
vh, fh = E.hyb_fugacity(x, P, T); vz, fz = E.zd_fugacity(x, P, T); vb, fb = E.fluid_fugacity(x, P, T)
s = (T-673)/100; w = s^2*(3-2s)
@printf("w %.3f  ln phi H2O: hyb %.4f ZD %.4f blend %.4f | V hyb %.3f ZD %.3f blend %.3f\n", w, log(fh[2]/(x[2]*P)), log(fz[2]/(x[2]*P)), log(fb[2]/(x[2]*P)), vh, vz, vb)
println("=== pure water, ZD09 vs hybrid at 0.3 GPa 500 K (why the hybrid exists)")
e = zeros(7); e[2] = 1.0
vz, fz = E.zd_fugacity(e, 3000.0, 500.0); vh, fh = E.hyb_fugacity(e, 3000.0, 500.0)
@printf("ZD09 V %.3f phi %.3e | hybrid V %.3f phi %.4e\n", vz, fz[2]/3000, vh, fh[2]/3000)
for (P,T) in ((20.0,500.0),(30.0,500.0),(100.0,500.0))
    av,bv = E.mrk_ab(T); A,B = av[2],bv[2]; rt=E.RK_R*T; st=√T; c1=-rt/P
    r = filter(x->x>B, E.cubic_roots(c1, c1*B + A/(st*P) - B^2, -A*B/(st*P)))
    if length(r) >= 2
        vmin,vmax = first(r), last(r)
        pdv = P*(vmax-vmin) - rt*log((vmax-B)/(vmin-B)) - A/(B*st)*log((vmax+B)/(vmin+B)*vmin/vmax)
        @printf("P %5.1f bar T %.0f K: roots %s  G(vap)-G(liq) %.2f J/mol*10 -> %s\n", P, T, string(round.(r,digits=2)), pdv, pdv>0 ? "liquid (Vmin)" : "vapour (Vmax)")
    else
        @printf("P %5.1f bar T %.0f K: single root %s\n", P, T, string(round.(r,digits=2)))
    end
end
