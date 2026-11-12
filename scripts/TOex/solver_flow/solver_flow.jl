# Rebuilds the inner levels of coh_state_u (carbon x, hydrogen y, fluid amount n) for one converged state, φ frozen,
# checks them against partition_Oₑₓ, and plots each budget against its unknown for TOex_solver_flow.pdf
include(joinpath(@__DIR__, "..", "solid_scratch", "SP_cairo.jl")); using .StagPP, MAGEMin_C, CairoMakie, Serialization, Printf
data = Initialize_MAGEMin("sb24", verbose=false); tb = load_StagH2Oa()
P, T, TOw, Cw, Hw = 3.0, 1636.0, 0.02, 0.05, 0.1                # harzburgite, no melt, wt% budgets
o = partition_Oₑₓ(P, T, 0.0, 0.0, 1e-2TOw, 1e-2Cw, 1e-2Hw; data=data, verbose=false, out=true, sH2Oa=tb)
Pbar, f = 1e4P, o.fO2; fO₂ = 10.0^f
K1, K2, K3, K4 = StagPP.reaction_Ks(T, Pbar)
# φ of the converged fluid (ZD09 order CH4, H2O, CO2, H2, CO, O2, C2H6)
vx = zeros(7); vx[1], vx[2], vx[3], vx[4], vx[5] = o.XCH4, o.XH2O, o.XCO2, o.XH2, o.XCO
_, fug = StagPP.fluid_fugacity(vx, Pbar, T); φ = [vx[i] > 0 ? fug[i]/(vx[i]*Pbar) : 1.0 for i in 1:7]
NC, NH = 10Cw/StagPP.mm.C, 2*10Hw/StagPP.mm.H2O                  # moles C and H in the 1 kg cell
n_sol = o.molH_fluid/(2o.XH2O + 4o.XCH4 + 2o.XH2)
fE = StagPP.logfO2_EMOG(T, Pbar); x_cap = log(min(1.0, 10.0^(fE - f)))
Ccurve = StagPP.solid_h2o_curve(tb, P, T, 0.0); Csol = (la = log.(tb.aw), lC = log.(max.(Ccurve, 1e-12))); φw0 = StagPP.pure_water_φ(P, T)
sp(x, y) = (aC = exp(x); fH2 = exp(y); s = sqrt(fO₂);
            (H2O = K3*fH2*s/(φ[2]*Pbar), CO2 = aC*K1*fO₂/(φ[3]*Pbar), CH4 = aC*K4*fH2^2/(φ[1]*Pbar), CO = aC*K2*s/(φ[5]*Pbar), H2 = fH2/(φ[4]*Pbar)))
G(x) = NC*exp(min(x - x_cap, 0.0)/(0.01*log(10)))               # graphite/carbonate sink
Cfl(x, y, n) = (s = sp(x, y); n*(s.CO2 + s.CO + s.CH4))
a_of(y) = min(1.0, K3*exp(y)*sqrt(fO₂)/(φw0*Pbar))
Hsol(y) = 2e-2*1000*StagPP.solid_h2o_wt(a_of(y), Csol)/StagPP.mm.H2O
Hfl(x, y, n) = (s = sp(x, y); n*(2s.H2O + 4s.CH4 + 2s.H2))
lres(d, t) = d <= 0 ? -1e3 : log(d/t)
x_of(y, n) = StagPP.safe_root(x -> lres(Cfl(x, y, n) + G(x), NC), x_cap - 700, x_cap; atol=1e-12).x
y_ref = log(Pbar/(1/φ[4] + K3*sqrt(fO₂)/φ[2]))
y_of(n) = StagPP.safe_root(y -> lres(Hsol(y) + Hfl(x_of(y, n), y, n), NH), y_ref - 140, y_ref + 20; atol=1e-12).x
S_of(n) = (y = y_of(n); s = sp(x_of(y, n), y); s.H2O + s.CO2 + s.CH4 + s.CO + s.H2)
n_r = exp(StagPP.safe_root(u -> log(S_of(exp(u))), log(NC + NH/2) - 70, log(NC + NH/2); atol=1e-12).x)
y_r = y_of(n_r); x_r = x_of(y_r, n_r)
@printf("solver : ΔFMQ %+.3f  n %.6e  aC %.6f  fH2 %.6e bar  shares C gph %.3f fl %.3f  H solid %.3f fl %.3f\n", o.FMQ, n_sol, o.aC, o.fH2, o.fC_gph, o.fC_fluid, o.fH_solid, o.fH_fluid)
@printf("replica: n %.6e  aC %.6f  fH2 %.6e bar  shares C gph %.3f fl %.3f  H solid %.3f fl %.3f\n", n_r, exp(x_r), exp(y_r),
        G(x_r)/NC, Cfl(x_r, y_r, n_r)/NC, Hsol(y_r)/NH, Hfl(x_r, y_r, n_r)/NH)
@printf("x_cap %.4f (log10 aC cap %.3f), fE %.3f, f %.3f, y_ref %.3f\n", x_cap, x_cap/log(10), fE, f, y_ref)

# Figures: (a) carbon vs log10 aC at the solved (y, n); (b) hydrogen vs log10 fH2 at the solved n, x re-solved per y;
# (c) S vs n; (d) oxygen demand vs ΔFMQ from the behaviour sweep (TOex at which each fO2 is the answer)
ink = colorant"#0b0b0b"; muted = colorant"#52514e"; cB = colorant"#2a78d6"; cS = colorant"#4d4d4d"; cG = colorant"#1baf7a"; cR = colorant"#d64a3b"
set_theme!(fontsize=15, textcolor=ink, Axis=(xgridvisible=false, ygridvisible=false, titlealign=:left, titlefont=:bold, titlesize=15))
fg = Figure(size=(1500, 1050))
ga = fg[1, 1] = GridLayout()
a1 = Axis(ga[1, 1], title="(a) carbon: solve x = ln a_C   (y, n fixed)", xlabel="log₁₀ a_C", ylabel="carbon demand / budget", yscale=log10)
lx = range(x_cap/log(10) - 3, x_cap/log(10), length=600); xs = lx .* log(10)
lines!(a1, lx, max.([Cfl(x, y_r, n_r) for x in xs] ./ NC, 1e-6), color=cB, linewidth=2.5, label="fluid  n·(X_CO₂+X_CO+X_CH₄)")
lines!(a1, lx, max.([G(x) for x in xs] ./ NC, 1e-6), color=cS, linewidth=2.5, label="graphite sink G(x)")
lines!(a1, lx, [(Cfl(x, y_r, n_r) + G(x)) for x in xs] ./ NC, color=ink, linewidth=3, linestyle=:dash, label="total")
hlines!(a1, 1.0, color=cR, linewidth=1.5); scatter!(a1, [x_r/log(10)], [1.0], color=cR, markersize=14)
ylims!(a1, 1e-3, 3); axislegend(a1, position=:lt, framevisible=false, labelsize=13)
ins = Axis(ga[1, 2], title="zoom: last 0.06 decade", xlabel="log₁₀ a_C", ylabel="demand / budget (linear)", titlesize=13, xticks=[-0.06, -0.03, 0.0])
colsize!(ga, 2, Relative(0.3))
lz = range(x_cap/log(10) - 0.06, x_cap/log(10), length=300); xz = lz .* log(10)
lines!(ins, lz, [Cfl(x, y_r, n_r) for x in xz] ./ NC, color=cB, linewidth=2); lines!(ins, lz, G.(xz) ./ NC, color=cS, linewidth=2)
lines!(ins, lz, [(Cfl(x, y_r, n_r) + G(x)) for x in xz] ./ NC, color=ink, linewidth=2.5, linestyle=:dash)
hlines!(ins, 1.0, color=cR); scatter!(ins, [x_r/log(10)], [1.0], color=cR, markersize=10)
a2 = Axis(fg[1, 2], title="(b) hydrogen: solve y = ln f_H₂   (n fixed, x re-solved at every y)", xlabel="log₁₀ f_H₂ (bar)", ylabel="hydrogen demand / budget", yscale=log10)
ly = range(y_r/log(10) - 3, y_r/log(10) + 1.0, length=300); ys = ly .* log(10)
xy = [x_of(y, n_r) for y in ys]
lines!(a2, ly, max.(Hsol.(ys) ./ NH, 1e-6), color=cG, linewidth=2.5, label="solid  C(a(y))")
lines!(a2, ly, max.([Hfl(xy[i], ys[i], n_r) for i in eachindex(ys)] ./ NH, 1e-6), color=cB, linewidth=2.5, label="fluid  n·(2X_H₂O+4X_CH₄+2X_H₂)")
lines!(a2, ly, [(Hsol(ys[i]) + Hfl(xy[i], ys[i], n_r)) for i in eachindex(ys)] ./ NH, color=ink, linewidth=3, linestyle=:dash, label="total")
hlines!(a2, 1.0, color=cR, linewidth=1.5); scatter!(a2, [y_r/log(10)], [1.0], color=cR, markersize=14)
ylims!(a2, 1e-3, 10); axislegend(a2, position=:lt, framevisible=false, labelsize=13)
a3 = Axis(fg[2, 1], title="(c) saturation: solve n   (x, y re-solved at every n)", xlabel="fluid amount n (mol per kg)", ylabel="S(n) = Σ Xᵢ", xscale=log10, yscale=log10)
ns = exp10.(range(log10(n_r) - 4, log10(NC + NH/2), length=160))
lines!(a3, ns, S_of.(ns), color=ink, linewidth=2.5); hlines!(a3, 1.0, color=cR, linewidth=1.5); scatter!(a3, [n_r], [1.0], color=cR, markersize=14)
a4 = Axis(fg[2, 2], title="(d) oxygen: solve fO₂   (everything above re-solved at every fO₂)", xlabel="ΔFMQ", ylabel="oxygen demand (wt%)",
          yscale=Makie.pseudolog10)
R0 = sort(filter(z -> z.P == 3.0 && z.ϕ == 0.0, deserialize(joinpath(@__DIR__, "..", "behaviour", "runs", "fo2.jls"))), by = z -> z.TO)
lines!(a4, [z.FMQ for z in R0], [z.TO for z in R0], color=ink, linewidth=2.5); hlines!(a4, TOw, color=cR, linewidth=1.5)
scatter!(a4, [o.FMQ], [TOw], color=cR, markersize=14)
Label(fg[0, 1:2], @sprintf("One state, each level against its own unknown: harzburgite, %g GPa, %g K, no melt, %g wt%% C, %g wt%% H₂O, TOex %g wt%% (red: budget and root)", P, T, Cw, Hw, TOw),
      font=:bold, fontsize=17)
save(joinpath(@__DIR__, "figs", "levels.png"), fg, px_per_unit=2); Finalize_MAGEMin(data); println("done")
