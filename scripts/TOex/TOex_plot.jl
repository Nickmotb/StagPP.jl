using StagPP, GLMakie, ColorSchemes, Printf

# ===============================================================================================
# TOex_plot.jl — donut + bar grid for a single TOex budget slice (mat[:,:,:,k])
# ===============================================================================================
#
# Layout: x = j (the 4 profiles: default / more melt / more carbon / more water)
#         y = i (the 4 PT levels, depth increasing downward)
#
# piechart=true   each cell is a two-ring donut with two stacked bars in the hole:
#                   outer ring        TOₑₓ sink shares      (sOₑₓ, carbO, mOₑₓ, carbonO, cohsink)
#                   inner ring        carbon partitioning   (fC_gph, fC_carb, fC_melt, fC_fluid)
#                   left list         COH fluid speciation  (XH2O, XCO2, XCH4, XCO, XH2)
#                   right list        solid iron speciation (Fe³⁺, Fe²⁺, Fe⁰diss, Fe⁰disp)
#                 fO₂ / ΔFMQ / ΔIW / melt XCO₂ sit in a white box above each donut.
# piechart=false  each cell is a four-column table of colour swatches and numbers, so the values
#                 can be read straight off the plot: % of TOₑₓ, mol% of the COH fluid, % of TC,
#                 % of solid ΣFe.
#
# RINGS vs LISTS. Only the two budget shares get rings: they are the shares whose wedges must be
# compared against each other, and a ring reads a part-of-whole best. The two speciations are
# printed lists instead — at 1/20 of the cell area a wedge or a bar segment cannot be measured by
# eye, so the number is the only honest way to show it, and the hole was free anyway.
#
# COLOUR DESIGN. carbO/carbonO/cohsink (outer ring) and fC_carb/fC_melt/fC_fluid (inner ring) are
# literally the SAME physical reservoir — solid carbonate, melt-dissolved carbon, and the COH
# fluid — read against two different budgets (O and C). They share a colour across the two rings
# on purpose, so a wedge can be tracked from the outer ring straight through to the inner one.
# graphite/diamond has no oxygen-budget counterpart (0 excess O per C) and gets its own colour —
# black, since that is what graphite is. Each list carries its own independent palette, because a
# mole fraction and a per-iron share are not budget shares and must not read as one.
# All palettes are hand-picked, not continuously sampled from a 10-stop colormap — sampling
# a qualitative scheme at only 5 of its stops interpolates BETWEEN categorical colours, which is
# what made "melt Fe³⁺" and "XCO₂" land on near-identical blends in an earlier version.
#
# The IRON list is the one exception to independence: its Fe³⁺ swatch is the outer ring's "solid
# Fe³⁺" reservoir seen per-iron instead of per-oxygen, so it keeps that colour. The two metallic
# swatches are one reservoir reached by two routes, so they are one hue at two lightnesses; Fe²⁺
# takes the only hue left unused. SOLID iron only — melt Fe³⁺ stays on the ring, where it is a sink.
#
# A grey "unplaced" entry appears on the OUTER ring whenever the sink shares don't sum to 1, which
# happens only in the :capacity_exceeded regime (TOex.jl classifies !converged first) — hiding it
# by renormalising would mask exactly the failure mode TOex.jl is built to surface. The carbon ring
# has no analogous case: carbon mass balance is satisfied exactly at every solve (TOex_math.txt
# §3), so fC_gph+fC_carb+fC_melt+fC_fluid is always either 1 (TC>0) or identically 0 (TC=0, drawn
# as an empty ring). Non-converged cells get a dashed red outline.
#
# The COH fluid share on the OUTER ring is SIGNED (TOex.jl:956): a reducing fluid speciates H₂O to
# CH₄/H₂ and hands that oxygen back to the rock, so its share goes negative while the rock shares
# exceed 1 — the five still sum to 1. That ring is therefore normalised by Σ|share|, which makes it
# the total oxygen *traffic*: every wedge is a real slice, the circle closes at exactly 360°, and
# the fluid stays in the pie whichever way it flows. A negative wedge is a SOURCE and is drawn
# tinted with a solid outline instead of a flat fill. Nothing else here is signed.

# Slot layout of `mat`'s third dimension, as filled by TOex_test.jl.
const IDX = (conv=1, fO2=2, fmq=3, sinks=4:8, coh=9:13, xco2=14, carbon=15:18, iron=19:22,
             nfe=23, iw=24)

# --- Outer ring: TOₑₓ sinks — sOₑₓ, carbO, mOₑₓ, carbonO, cohsink ---
const clr_sFe   = colorant"#1f77b4"   # solid Fe³⁺    — O-budget only
const clr_carb  = colorant"#bcbd22"   # solid carbonate      \  shared with the carbon ring:
const clr_mCO2  = colorant"#d62728"   # melt CO₂ / melt C     > these are the same reservoir,
const clr_fluid = colorant"#17becf"   # COH fluid net O / TC /  read against two different budgets
const clr_mFe   = colorant"#9467bd"   # melt Fe³⁺     — O-budget only
const clr_gph   = colorant"#2b2b2b"   # graphite/diamond — C-budget only (0 excess O per C)
const unplaced_color = RGBAf(0.75, 0.75, 0.75, 1.0)

const toex_colors   = [clr_sFe, clr_carb, clr_mFe, clr_mCO2, clr_fluid]   # sOₑₓ,carbO,mOₑₓ,carbonO,cohsink
const carbon_colors = [clr_gph, clr_carb, clr_mCO2, colorant"#1f77b4"]            # fC_gph,fC_carb,fC_melt,fC_fluid

# --- Left list: COH speciation — an independent palette (a mole fraction, not a budget share) ---
const coh_colors = [colorant"#17becf", colorant"#ff7f0e", colorant"#8c564b",
                     colorant"#e377c2", colorant"#7f7f7f"]                 # XH2O,XCO2,XCH4,XCO,XH2

# --- Right list: solid iron speciation — a share of ΣFe, not of a budget ---
const clr_FeIII = clr_sFe                    # the TOₑₓ ring's solid Fe³⁺, counted per iron
const clr_FeII  = colorant"#2ca02c"          # green — the one hue unused elsewhere in the figure
const clr_Fe0d  = colorant"#5c5346"          # dark iron-brown  \ one metal, two routes: same
const clr_Fe0p  = colorant"#9a8b73"          # light iron-brown / family, two lightnesses
const iron_colors = [clr_FeIII, clr_FeII, clr_Fe0d, clr_Fe0p]

const toex_names   = ["solid Fe³⁺", "solid carbonate", "melt Fe³⁺", "melt CO₂", "COH fluid"]
const coh_names    = ["H₂O", "CO₂", "CH₄", "CO", "H₂"]
const carbon_names = ["graphite/diamond", "solid carbonate", "melt", "COH fluid"]
const iron_names   = ["Fe³⁺ (Fe₂O₃)", "Fe²⁺ (FeO)", "Fe⁰ — dissociation", "Fe⁰ — disproportionation"]

# Fixed point where it is readable, scientific where it would otherwise print as 0.00 or overflow.
function fmtval(x)
    a = abs(x)
    a == 0 && return "0"
    (a < 1e-2 || a >= 1e4) && return @sprintf("%.1e", x)
    return @sprintf("%.2f", x)
end

# Percent for the in-hole lists, floored instead of exponential: at this size 1.1e-10% and 0%
# say the same thing, and an exponent costs the width two real digits would use.
function fmtpct(x)
    a = abs(x)
    a == 0    && return "0%"
    a < 1e-2  && return "<0.01%"
    a >= 99.995 && return @sprintf("%.0f%%", x)
    return @sprintf("%.2f%%", x)
end

# ---------------------------------------------------------------------------
# Geometry: one filled annulus sector (or a plain wedge when r_in == 0)
# ---------------------------------------------------------------------------
function wedge_points(center::Point2f, r_in::Real, r_out::Real, θ0::Real, θ1::Real; n=48)
    outer = [center + Point2f(r_out*cos(θ), r_out*sin(θ)) for θ in range(θ0, θ1, length=n)]
    r_in <= 0 && return vcat([center], outer)
    inner = [center + Point2f(r_in*cos(θ), r_in*sin(θ)) for θ in range(θ1, θ0, length=n)]
    return vcat(outer, inner)
end

circle_points(center::Point2f, r::Real; n=64) = [center + Point2f(r*cos(θ), r*sin(θ)) for θ in range(0, 2π, length=n)]

# ---------------------------------------------------------------------------
# One ring: `fracs` drawn clockwise from 12 o'clock between radii r_in/r_out
# ---------------------------------------------------------------------------
function draw_ring!(ax, center::Point2f, r_in::Real, r_out::Real, fracs::AbstractVector,
                     colors::AbstractVector; start_angle=π/2, mincut=1e-6)
    # Σ|f|, not Σf: a negative share must still occupy arc, or the ring cannot close. Sweeping a
    # raw signed share runs backwards and self-overlaps, which is what tore the earlier plots up.
    total = sum(abs, fracs)
    total <= 0 && return nothing
    θ = start_angle
    for (f, c) in zip(fracs, colors)
        frac = abs(f)/total
        θ1 = θ - 2π*frac                      # clockwise
        if frac >= mincut
            src = f < 0                       # a source, not a sink (never true for carbon fracs)
            poly!(ax, wedge_points(center, r_in, r_out, θ, θ1),
                  color       = src ? (c, 0.30) : c,
                  strokecolor = src ? c : :white,
                  strokewidth = src ? 0.8 : 0.25)
        end
        θ = θ1
    end
end

# Solid Fe speciation as shares of ΣFe. Fe²⁺ is the remainder, so the four close on 1 by
# construction; an ironless cell returns zeros and draws as an empty list, like TC=0 does.
function iron_entries(FeR, FeM, FeM_diss, FeM_disp, nFe)
    nFe <= 0 && return [0.0, 0.0, 0.0, 0.0]
    return [FeR, max(0.0, 1 - FeR - FeM), FeM_diss, FeM_disp]
end

# Sink shares plus the shortfall, which is nonzero only when the solve did not converge.
function sink_entries(toex_fracs; mincut=1e-6)
    vals, cols = collect(toex_fracs), collect(toex_colors)
    unplaced = max(0.0, 1.0 - sum(toex_fracs))
    if unplaced > mincut
        push!(vals, unplaced); push!(cols, unplaced_color)
    end
    return vals, cols
end

# ---------------------------------------------------------------------------
# One vertical swatch + percent list, top-down, headed by `label`
# Values are printed raw: both lists already sum to 1, so a shortfall should show rather than hide.
# ---------------------------------------------------------------------------
function draw_list!(ax, x, y, fracs, colors; label="", dy=0.04, sq=0.026, fontsize=6.5)
    sum(abs, fracs) <= 0 && return nothing
    label != "" && text!(ax, Point2f(x, y), text=label, align=(:left, :center),
                         fontsize=fontsize, font=:bold, color=(:black, 0.6))
    for (n, (f, c)) in enumerate(zip(fracs, colors))
        yr = y - n*dy
        poly!(ax, Rect2f(x, yr - sq/2, sq, sq), color=c, strokecolor=:black, strokewidth=0.2)
        text!(ax, Point2f(x + sq + 0.007, yr), text=fmtpct(100f),
              align=(:left, :center), fontsize=fontsize)
    end
end

# ---------------------------------------------------------------------------
# Pie mode: TOₑₓ ring + carbon ring, COH and iron lists in the hole, label box above
# The box sits ~3x closer to its own donut than to the row above, so attribution is unambiguous.
# ---------------------------------------------------------------------------
function draw_glyph!(ax, center::Point2f; toex_fracs, coh_fracs, carbon_fracs, iron_fracs,
                      fO2, ΔFMQ, ΔIW, XCO2melt, converged=true,
                      r_c=(0.205, 0.275), r_o=(0.29, 0.37), ypie=-0.105,
                      lx=0.138, ly=0.10, ldy=0.04, fontsize=7.5, listfont=6.5, mincut=1e-6)
    cd = Point2f(center[1], center[2] + ypie)

    vals, cols = sink_entries(toex_fracs; mincut)
    draw_ring!(ax, cd, r_o[1], r_o[2], vals, cols; mincut)
    draw_ring!(ax, cd, r_c[1], r_c[2], carbon_fracs, carbon_colors; mincut)

    # Six rows of two columns fit the hole: the widest corner lands at r=0.18 against a carbon
    # ring starting at 0.205. The 0.152 column pitch leaves 0.033 clear between the two.
    draw_list!(ax, cd[1]-lx,       cd[2]+ly, coh_fracs,  coh_colors;  label="COH", dy=ldy, fontsize=listfont)
    draw_list!(ax, cd[1]-lx+0.152, cd[2]+ly, iron_fracs, iron_colors; label="Fe",  dy=ldy, fontsize=listfont)

    bx, by, bwid, bhgt = center[1] - 0.33, center[2] + 0.30, 0.66, 0.145
    poly!(ax, Rect2f(bx, by, bwid, bhgt), color=:white, strokecolor=(:black, 0.35), strokewidth=0.6)
    text!(ax, Point2f(center[1], by + 0.100), align=(:center, :center), fontsize=fontsize,
          text=@sprintf("fO₂ = %.2f     XCO₂ = %s", fO2, fmtval(XCO2melt)))
    text!(ax, Point2f(center[1], by + 0.043), align=(:center, :center), fontsize=fontsize,
          text=@sprintf("ΔFMQ = %+.2f     ΔIW = %+.2f", ΔFMQ, ΔIW))

    !converged && lines!(ax, circle_points(cd, r_o[2]*1.05), color=:red,
                         linewidth=2, linestyle=:dash)
end

# ---------------------------------------------------------------------------
# Table mode: four columns of swatch + number, readable straight off the plot
# ---------------------------------------------------------------------------
function draw_table!(ax, center::Point2f; toex_fracs, coh_fracs, carbon_fracs, iron_fracs,
                      fO2, ΔFMQ, ΔIW, XCO2melt, converged=true, sq=0.05, fontsize=10,
                      headersize=12, mincut=1e-6, panel=true)
    cx, cy = center[1], center[2]

    # Each solve sits on its own card, widened to take the fourth column. Cards are pulled in from
    # the full 1.0 grid pitch so the whitespace lands between simulations, not inside one.
    panel && poly!(ax, Rect2f(cx - 0.46, cy - 0.36, 0.92, 0.78), color=:white,
                   strokecolor=(:black, 0.25), strokewidth=0.5)

    text!(ax, Point2f(cx, cy + 0.37), align=(:center, :top), fontsize=headersize, font=:bold,
          text=@sprintf("fO₂=%.2f   ΔFMQ=%+.2f   ΔIW=%+.2f", fO2, ΔFMQ, ΔIW))
    text!(ax, Point2f(cx, cy + 0.28), align=(:center, :top), fontsize=headersize,
          text="melt XCO₂ = " * fmtval(XCO2melt))

    vals, cols = sink_entries(toex_fracs; mincut)
    swatch!(x, y, c) = poly!(ax, Rect2f(x, y - sq/2, sq, sq), color=c,
                             strokecolor=:black, strokewidth=0.25)
    # Start-to-start column pitch, not cell-wide; tightened from 0.32 to seat the iron column.
    xs = cx .+ (-0.445, -0.22, 0.005, 0.23)
    y0, dy = cy + 0.14, 0.075

    # One-word headers: with four columns the legend alone no longer says which is which.
    for (x, h) in zip(xs, ("TOₑₓ", "COH", "C", "Fe"))
        text!(ax, Point2f(x, y0 + 0.05), text=h, align=(:left, :center),
              fontsize=fontsize-1, font=:bold, color=(:black, 0.55))
    end

    cols_of = (cols, coh_colors, carbon_colors, iron_colors)
    vals_of = (vals, coh_fracs, carbon_fracs, iron_fracs)
    for r in 1:maximum(length, vals_of)
        y = y0 - (r-1)*dy
        for (x, v, c) in zip(xs, vals_of, cols_of)
            r > length(v) && continue
            swatch!(x, y, c[r])
            text!(ax, Point2f(x + sq + 0.016, y), text=fmtval(100*v[r]) * "%",
                  align=(:left, :center), fontsize=fontsize)
        end
    end

    if !converged
        box = [Point2f(cx-0.46, cy-0.36), Point2f(cx+0.46, cy-0.36), Point2f(cx+0.46, cy+0.42),
               Point2f(cx-0.46, cy+0.42), Point2f(cx-0.46, cy-0.36)]
        lines!(ax, box, color=:red, linewidth=1.5, linestyle=:dash)
    end
end

# ---------------------------------------------------------------------------
# Full grid for one TOex budget k
# ---------------------------------------------------------------------------
function plot_toex_grid(mat, PT, TOex, ϕ, p, TC, H; k=1, piechart=true,
                         glyph_kwargs=NamedTuple(), savein="")
    nrow, ncol = size(mat, 1), size(mat, 2)
    row_labels = ["$(round(Int, GPa2km(PT[i,1]))) km\n$(Int(PT[i,2])) K" for i in 1:nrow]

    # Row i is drawn at nrow+1-i so depth increases downward. Done by placement rather than
    # yreversed=true, which would mirror the axis and flip every pie to counter-clockwise.
    ypos(i) = nrow + 1 - i

    col_labels = ["default", "↑ melt fraction ($(1e2(ϕ[2]))%)", "↑ carbon($(TC[2]) wt%)", "↑ water($(H[2]) wt%)"]

    fig = Figure(size = piechart ? (1560, 1180) : (1980, 1340))
    ax = Axis(fig[1, 1], xticks=(1:ncol, col_labels), yticks=(1:nrow, reverse(row_labels)),
              aspect = piechart ? DataAspect() : nothing, xgridvisible=false, ygridvisible=false,
              backgroundcolor= :gray90,
              title="--- Default: ϕ=$(1e2ϕ[1])% | Xsolid=$p | Carbon=$(TC[1]) wt% | Water=$(H[1]) wt% | TOₑₓ=$(TOex[k]) wt% ---")

    for i in 1:nrow, j in 1:ncol
        args = (toex_fracs = mat[i, j, IDX.sinks, k], coh_fracs = mat[i, j, IDX.coh, k],
                carbon_fracs = mat[i, j, IDX.carbon, k],
                iron_fracs = iron_entries(mat[i, j, IDX.iron, k]..., mat[i, j, IDX.nfe, k]),
                fO2 = mat[i, j, IDX.fO2, k], ΔFMQ = mat[i, j, IDX.fmq, k],
                ΔIW = mat[i, j, IDX.iw, k],
                XCO2melt = mat[i, j, IDX.xco2, k], converged = mat[i, j, IDX.conv, k] == 1)
        c = Point2f(j, ypos(i))
        piechart ? draw_glyph!(ax, c; args..., glyph_kwargs...) :
                   draw_table!(ax, c; args..., glyph_kwargs...)
    end
    xlims!(ax, 0.4, ncol + 0.6); ylims!(ax, 0.4, nrow + 0.6)

    # Legends in a 2x2 block: rings on the top row, bars on the bottom, left/right matching where
    # each bar sits in the glyph. Patches and text enlarged — the old single column was unreadable.
    legkw = (patchsize=(20, 20), labelsize=12, titlesize=13, rowgap=3, padding=(10, 10, 10, 10),
             valign=:top)

    sink_elems = [PolyElement(color=c) for c in toex_colors]
    push!(sink_elems, PolyElement(color=(toex_colors[5], 0.30), strokecolor=toex_colors[5], strokewidth=0.8))
    push!(sink_elems, PolyElement(color=unplaced_color))
    ttl = piechart ? ("Outer ring — TOₑₓ sinks (% of TOₑₓ)", "Inner ring — carbon (% of TC)",
                      "Left list — COH speciation (mol%)", "Right list — solid iron (% of ΣFe)") :
                     ("TOₑₓ sinks (% of TOₑₓ)", "C — carbon partitioning (% of TC)",
                      "COH — fluid speciation (mol%)", "Fe — solid iron (% of ΣFe)")

    lg = GridLayout(fig[1, 2])
    Legend(lg[1, 1], sink_elems, vcat(toex_names, "COH fluid (source)", "unplaced"), ttl[1]; legkw...)
    Legend(lg[1, 2], [PolyElement(color=c) for c in carbon_colors], carbon_names, ttl[2]; legkw...)
    Legend(lg[2, 1], [PolyElement(color=c) for c in coh_colors],    coh_names,    ttl[3]; legkw...)
    Legend(lg[2, 2], [PolyElement(color=c) for c in iron_colors],   iron_names,   ttl[4]; legkw...)
    # Content-sized rows plus an empty row to absorb the slack. Without the spacer the two legend
    # rows split the figure height between them and the pairs drift half a page apart.
    lg[3, 1:2] = GridLayout()
    colgap!(lg, 8); rowgap!(lg, 6)
    rowsize!(lg, 1, Auto(false)); rowsize!(lg, 2, Auto(false)); rowsize!(lg, 3, Auto(true))
    (savein!="") && (GLMakie.save("./scripts/TOex/$savein.png", fig, px_per_unit=3.0))
    fig
end

# ---------------------------------------------------------------------------
# Example call — run TOex_test.jl first so `mat`, `PT`, `TOex` exist, then:
#     fig = plot_toex_grid(mat, PT, TOex, ϕ, p, TC, H; k=1, piechart=true)
#     display(fig)
# ---------------------------------------------------------------------------

fig = plot_toex_grid(mat, PT, TOex, ϕ, p, TC, H; k=1, piechart=true, savein="")
display(fig)
