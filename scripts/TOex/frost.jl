using MAGEMin_C, Printf
const mmFeO, mmFe, mmO = 71.85, 55.85, 16.0
Xox = ["SiO2","MgO","FeO","CaO","Al2O3","Na2O","Cr2O3","O"]
oxFe= ["SiO2","MgO","Fe" ,"CaO","Al2O3","Na2O","Cr2O3","O"]
XH  = [43.43, 45.93, 8.34, 0.90, 1.00, 0.01, 0.30, 0.0]   # harzburgite (depleted)
XF  = [45.00, 38.30, 8.20, 3.50, 4.50, 0.35, 0.40, 0.0]   # fertile lherzolite, Al-richer
isFe(p) = lowercase(p) in ("feg","fea","fee","iron","fe")

# TOex.jl's own buffer expressions
fmq(T,Pb) = -25096.3/T + 8.735 + 0.11*(Pb-1)/T
iw(T,Pb)  = -27489.0/T + 6.702 + 0.055*(Pb-1)/T

# bulk in sb24's Fe+O basis at a GIVEN Fe3+/SumFe, held FIXED while P varies
function bulk_at(X, r)
    nFe = X[3]/mmFeO
    x = copy(X); x[3] = nFe*mmFe; x[8] = nFe*(1 + 0.5r)*mmO
    return 1e2 .* x ./ sum(x)
end

d = Initialize_MAGEMin("sb24", verbose=false)
for (nm, X) in (("harzburgite (Al2O3 1.0)", XH), ("fertile lherzolite (Al2O3 4.5)", XF)), r in (0.03,)
    println("\n", "="^104)
    @printf("%s   bulk Fe3+/SumFe held FIXED at %.3f, adiabat T = 1550 + 12P\n", nm, r)
    println("="^104)
    @printf("%6s %6s %10s %9s %9s %11s %10s   %s\n","P GPa","T K","fO2","dFMQ","dIW","metal wt%","Fe0/SFe","assemblage")
    for P in (1.0,2.0,3.0,5.0,7.0,8.0,10.0,12.0,14.0,16.0,20.0,25.0)
        T  = 1550.0 + 12.0P
        o  = single_point_minimization(10P, T-273.15, d, X=bulk_at(X,r), Xoxides=oxFe, sys_in="wt")
        mf = sum(o.ph_frac_wt[i] for i in eachindex(o.ph) if isFe(o.ph[i]); init=0.0)
        r0 = (mf/mmFe)/((X[3]/sum(X))/mmFeO)
        @printf("%6.1f %6.0f %10.3f %+9.3f %+9.3f %11.5f %10.5f   %s\n",
                P, T, o.fO2, o.fO2-fmq(T,1e4P), o.fO2-iw(T,1e4P), 1e2mf, r0, join(o.ph," "))
    end
end
Finalize_MAGEMin(d)
