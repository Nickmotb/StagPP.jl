# Julia reference for the Fortran port: coh_state_u and solid_h2o_curve at a fixed state set
include(joinpath(homedir(), "Desktop/StagPP.jl/scripts/TOex/solid_scratch/SP_cairo.jl")); using .StagPP, Printf
tb = StagPP.load_StagH2Oa(joinpath(homedir(), "Desktop/StagPP.jl/scripts/sH2O/out/aw16/StagH2Oa.dat"))
Tad(P) = P <= 23.5 ? 1600 + 12P : 1600 + 12*23.5 + 5*(P - 23.5)
Tsl(P) = min(700 + 60P, Tad(P) - 300)
st = open("states.txt", "w"); out = open("j_coh.txt", "w"); cv = open("j_curve.txt", "w")
k = 0; t = 0.0
for P in (0.5, 1.0, 3.0, 6.0, 10.0, 15.0, 25.0, 60.0), Tk in (:ad, :slab), p in (0.0, 1.0), ϕ in (0.0, 0.02),
    H in (1e-5, 1e-4, 1e-3, 5e-3), d in (-3.0, -1.0, 0.0, 1.0)
    global k += 1
    T = Tk === :ad ? Tad(P) : Tsl(P); TC = 1e-4; TOex = 1e-3
    (; X, XB, mmXox) = StagPP.endmember_composition(p)
    (; Mf, Ms, Mc, Mh, TOₑₓ, molfracXB, nB_per_mass) = StagPP.mass_budget(1e3, ϕ, TC, TOex, H, XB, mmXox)
    Csat = StagPP.solid_h2o_curve(tb, P, T, p)
    sol = (; Ms, C = (la = log.(tb.aw), lC = log.(max.(Csat, 1e-12))), φw0 = P <= 50.0 ? StagPP.pure_water_φ(P, T) : 1.0)
    f = StagPP.logfO2_EMOG(T, 1e4P) + d
    global t += @elapsed cs = StagPP.coh_state_u(f, T, P, Mh, Mc, Mf*nB_per_mass, TOₑₓ; Mf=Mf, molfracXB=molfracXB, sol=sol)
    @printf(st, "%d %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n", k, f, P, T, p, ϕ, TOex, TC, H)
    @printf(out, "%d", k)
    for v in (cs.aC, cs.n_fluid, cs.XH2O, cs.XCO2, cs.XCH4, cs.XCO, cs.XH2, cs.nC_fluid, cs.nC_melt, cs.nC_solid,
              cs.nC_carb, cs.molH_fl, cs.nH_solid, cs.nH2_melt, cs.nH_H2O_melt, cs.a_solid, cs.fH2, cs.fₑₓ, cs.S_fluid,
              Float64(cs.it_ideality), cs.vol_fluid, cs.res_ideality)
        @printf(out, " %.17g", v)
    end
    println(out)
    @printf(cv, "%d %.17g %.17g %.17g %.17g\n", k, Csat[1], Csat[6], Csat[11], Csat[16])
end
close(st); close(out); close(cv)
@printf("%d states, Julia coh_state_u %.3f ms per call\n", k, 1e3t/k)
