# Demand plots for TOex_solver_flow.tex: the oxygen budget each fO2 needs, and the carbon and hydrogen
#  budgets split among their sinks against fO2, for three melt fractions at one P-T point.
# Run: julia --project=<StagPP.jl> demand_plots.jl  (writes figs/demand_{O,C,H}.png next to this file)

isdefined(Main, :TOexLite) || @eval using StagPP     # a Makie-free loader may stand in for StagPP
const M = isdefined(Main, :TOexLite) ? Main.TOexLite : StagPP
using CairoMakie, MAGEMin_C, Printf

# ---- the cell: 3 GPa, 1700 K, 20 % basalt, 0.1 wt% carbon, 0.3 wt% water; melt fraction varies ----
const P, T, p   = 3.0, 1700.0, 0.2
const TC, H     = 1e-3, 3e-3                       # mass fractions
const TOex0     = 1e-3                             # oxygen budget of the reference state (mass fraction)
const PHIS      = (0.0, 0.01, 0.1)                 # melt fractions shown
const COLS      = ("#0072B2", "#E69F00", "#009E73")   # Okabe-Ito, colour-blind safe
const STYLES    = (:solid, :dash, :dot)            # second encoding, so colour is never the only cue
const LABELS    = ("no melt", "1 % melt", "10 % melt")
const OUT       = joinpath(@__DIR__, "figs")
isdir(OUT) || mkpath(OUT)
data = Initialize_MAGEMin("sb24", verbose=false)

const MAGC = Dict{Any,Any}()                          # one MAGEMin multi_point per (P, T, bulk), reused by every solve
part(ϕ, TOex) = M.partition_Oₑₓ(P, T, p, ϕ, TOex, TC, H; verbose=false, data=data, out=true, mag_cache=MAGC)

# ---- one sweep through the oxygen budget per melt fraction, shared by all three plots ----
budgets = 10 .^ range(-5, log10(5e-3), length=30)
println("loaded, sweeping"); flush(stdout)
function solve_logged(ϕ, b)
    t = @elapsed o = part(ϕ, b)
    @printf("ϕ=%.2f TOex=%.2e wt%%  ΔFMQ=%6.2f  FeR_solid=%.3f  solid=%.2f carb=%.2f coh=%.2f  conv=%s  %.2fs\n",
            ϕ, 100b, o.FMQ, o.FeR_solid, o.sOₑₓ, o.carbO, o.cohsink, o.converged, t); flush(stdout)
    return o
end
sweep = [sort!(filter(o -> o.converged, [solve_logged(ϕ, b) for b in budgets]), by = o -> o.FMQ) for ϕ in PHIS]

# ---- 1. oxygen: the budget each fO2 needs (the oxygen demand at the solved state) ----
figO = Figure(size=(520, 370), fontsize=13)
axO  = Axis(figO[1, 1], xlabel="ΔFMQ", ylabel="excess-O budget [wt%]", yscale=log10,
            title="Oxygen demand against fO₂", xgridvisible=false, ygridcolor=(:gray, 0.25))
for (k, r) in enumerate(sweep)
    lines!(axO, [o.FMQ for o in r], [1e2*o.TOₑₓ_g/1e3 for o in r]; color=COLS[k], linestyle=STYLES[k], linewidth=2, label=LABELS[k])
end
Legend(figO[2, 1], axO; orientation=:horizontal, framevisible=false)   # below the plot, clear of the lines
save(joinpath(OUT, "demand_O.png"), figO; px_per_unit=3)

# ---- 2-3. carbon and hydrogen: each budget's split among its sinks, swept through fO2 by the oxygen budget ----
const SCOLS = ("#0072B2", "#E69F00", "#009E73", "#CC79A7")   # Okabe-Ito, one colour per sink
const SSTY  = (:solid, :dash, :dot, :dashdot)                 # second encoding per sink
function shares_fig(title, ylab, sinks, fname)
    fig = Figure(size=(1050, 380), fontsize=13); ax1 = nothing
    for (k, r) in enumerate(sweep)
        ax = Axis(fig[1, k], xlabel="ΔFMQ", ylabel=k == 1 ? ylab : "", title=LABELS[k],
                  xgridvisible=false, ygridcolor=(:gray, 0.25))
        ylims!(ax, -0.02, 1.02); k == 1 && (ax1 = ax)
        x = [o.FMQ for o in r]
        for (j, (lab, g)) in enumerate(sinks)
            lines!(ax, x, g.(r); color=SCOLS[j], linestyle=SSTY[j], linewidth=2, label=lab)
        end
    end
    Label(fig[0, 1:length(sweep)], title, font=:bold)
    Legend(fig[2, 1:length(sweep)], ax1; orientation=:horizontal, framevisible=false)
    save(joinpath(OUT, fname), fig; px_per_unit=3)
end

shares_fig("Carbon demand by sink against fO₂", "share of the carbon budget",
           (("melt CO₂", o -> o.fC_melt), ("fluid", o -> o.fC_fluid),
            ("graphite / diamond", o -> o.fC_gph), ("carbonate", o -> o.fC_carb)), "demand_C.png")
h2m(o) = o.molH_total > 0 ? 2o.nH2_melt/o.molH_total : 0.0         # melt H2, as a share of the H held
shares_fig("Hydrogen demand by sink against fO₂", "share of the hydrogen budget",
           (("melt H₂", h2m), ("melt H₂O + OH", o -> o.fH_melt - h2m(o)),
            ("solid", o -> o.fH_solid), ("fluid", o -> o.fH_fluid)), "demand_H.png")

Finalize_MAGEMin(data)
println("demand plots written to ", OUT)
