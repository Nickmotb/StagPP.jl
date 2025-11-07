# time.dat field --> String
LTN = Dict(
    "time"  => L"Evolution\;time\;[\mathrm{Gyr}]",
    "F_top" => L"Surface\;heat\;flux\;[\mathrm{mW/m^2}]",
    "F_bot" => L"CMB\;heat\;flux\;[\mathrm{mW/m^2}]",
    "Tmean" => L"Mean\;mantle\;temperature\;[\mathrm{K}]",
    "Vrms"  => L"RMS\;velocity\;[\mathrm{cm/yr}]",
    "eta_amean" => L"Mean\;viscosity\;[\mathrm{Pa\cdot s}]",
    "eta_mean"  => L"Mean\;viscosity\;[\mathrm{Pa\cdot s}]",
    "erupta" => L"Erupted\;material\;[\mathrm{kg}]",
    "erupt_rate" => L"Eruption\;rate\;[\mathrm{kg/s}]",
    "Tcmb" => L"CMB\;temperature\;[\mathrm{K}]",
    "Tsurf" => L"Surface\;temperature\;[\mathrm{K}]",
    "Tpotl" => L"Potential\;temperature\;[\mathrm{K}]",
    "cH2O_mean" => L"Mean\;mantle\;H_2O\;content\;[\mathrm{wt%}]",
    "IngassedH2O" => L"Total\;ingassed\;H_2O\;[\mathrm{kg}]",
    "OutgassedH2O" => L"Total\;outgassed\;H_2O\;[\mathrm{kg}]",
    "SaturationOutgassH2O" => L"Saturation\;outgassed\;H_2O\;[\mathrm{kg}]",
    "EruptedH2O" => L"Erupted\;H_2O\;[\mathrm{kg}]",
    "SurfOceanMass3D" => L"Surface\;ocean\;mass\;[\mathrm{kg}]",
    "TotH2OMassOceanPlusMantle" => L"Total\;system\;H_2O\;[\mathrm{kg}]",
    "Mob" => L"Lithosphere\;mobility\;[\mathrm{v_{surf}/v_{rms}}]",
    "Vsurf" => L"Surface\;velocity\;[\mathrm{cm/yr}]",
)

# rprof.dat field --> String
LRN = Dict(
    "Tmean" => L"Mean\;mantle\;temperature\;[\mathrm{K}]", "logTmean" => L"Mean\;mantle\;temperature\;[log_{10}K]",
    "Vrms"  => L"RMS\;velocity\;[\mathrm{cm/yr}]", "logvrms" => L"RMS\;velocity\;[log_{10}(\mathrm{cm/yr})]",
    "eta_amean" => L"Mean\;viscosity\;[\mathrm{Pa\cdot s}]", "logeta_amean" => L"Mean\;viscosity\;[log_{10}(\mathrm{Pa\cdot s})]",
    "eta_mean"  => L"Mean\;viscosity\;[\mathrm{Pa\cdot s}]", "logeta_mean"  => L"Mean\;viscosity\;[log_{10}(\mathrm{Pa\cdot s})]",
    "Water" => L"H_2O\;content\;[\mathrm{wt%}]", "logWater" => L"H_2O\;content\;[log_{10}(\mathrm{wt%})]",
    "Wsol" => L"H_2O\;storage\;capacity\;[\mathrm{wt%}]", "logWsol" => L"H_2O\;storage\;capacity\;[log_{10}(\mathrm{wt%})]",
    "fO2" => L"fO_2\;[log_{10}FMQ]", "logfO2" => L"fO_2\;[log_{10}FMQ]",
    "H2Ofree" => L"Free\;H_2O\;[\mathrm{wt%}]", "logH2Ofree" => L"Free\;H_2O\;[log_{10}(\mathrm{wt%})]",
    "rhomean" => L"Mean\;density\;[\mathrm{kg/m^3}]", "logrhomean" => L"Mean\;density\;[log_{10}(\mathrm{kg/m^3})]",
    "bsmean" => L"Mean\;basalt\;fraction", "logbsmean" => L"Mean\;basalt\;fraction\;[\mathrm{log_{10}}]",
    "etalog" => L"Viscosity\;[\mathrm{Pa\cdot s}]", "logetalog" => L"Viscosity\;[log_{10}(\mathrm{Pa\cdot s})]",
)

# Colormap to field dictionary
CTF = Dict(
    "Temperature" => :vik100,
    "T_residual"  => :vik100,
    "Viscosity"   => :davos,
    "Water"       => :Blues,
    "Water_solubility" => Reverse(:grays),
    "fO2"         => :plasma,
    "Density"     => :turbo,
    "Basalt"      => Reverse(:RdGy),
    "Velocity"    => :viridis,
    "Melt_fraction" => :default,
    "dmelt/dt"      => :default,
)

# =============================
# ==== Post-Processing API ====
# =============================

"""
    Plot the evolution of a field over time.
    
    \t Basic usage: \t time_vs_field(Dblock::DataBlock, field::String; kwargs...)

    Optional arguments (kwargs):

      • -- Canvas arguments --

        - fsize::Tuple{Int64,Int64} \t-->\t Figure size in pixels [default: (800, 600)]
        - xlabelsize::Int64 \t\t-->\t X-axis label font size [default: 25]
        - ylabelsize::Int64 \t\t-->\t Y-axis label font size [default: 25]
        - xgrid::Bool \t\t\t-->\t Show x-grid [default: true]
        - ygrid::Bool \t\t\t-->\t Show y-grid [default: true]
        - xticklabelsize::Int64 \t\t-->\t X-axis tick label font size [default: 17]
        - yticklabelsize::Int64 \t\t-->\t Y-axis tick label font size [default: 17]
        - xticksize::Int64 \t\t-->\t X-axis tick size [default: 10]
        - yticksize::Int64 \t\t-->\t Y-axis tick size [default: 10]
        - xlabelpadding::Int64 \t\t-->\t X-axis label padding [default: 15]
        - ylabelpadding::Int64 \t\t-->\t Y-axis label padding [default: 15]
        - yreversed::Bool \t\t-->\t Reverse y-axis [default: false]
        - savein::String \t\t-->\t Save figure in file [default: "" (no saving)]

      • -- Data arguments --

        - subsample::Int64 \t\t-->\t Subsample factor for plotting [default: 0 (automatic)]
        - tstart::Union{Float64, Nothing} -->\t Start time [default: nothing (first time)]
        - tend::Union{Float64, Nothing} \t-->\t End time [default: nothing (last time)]
        - color::Symbol \t\t\t-->\t Line/scatter color [default: :blue]
        - scatter::Bool \t\t\t-->\t Enables scatters [default: false]
            - markersize::Int64 \t\t-->\t Marker size [default: 10]
            - marker::Symbol \t\t-->\t Marker type [default: :circle]
        - line::Bool \t\t\t-->\t Enables lines [default: true]
            - linewidth::Float64 \t\t-->\t Line width [default: 2.5]
            - linestyle::Symbol \t\t-->\t Line style [default: :solid]
        - mov_avg::Bool \t\t\t-->\t Enables moving average [default: false]
            - mov_avg_window::Int64 \t-->\t Moving average window [default: 0 (automatic)]

"""
function time_vs_field(Dblock::DataBlock, field::String; fsize=(800, 600), subsample=0, color=:blue, xlabelsize=25, ylabelsize=25, xgrid=true, ygrid=true,
                            xticklabelsize=17, yticklabelsize=17, xticksize=10, yticksize=10, xlabelpadding=15, ylabelpadding=15, yreversed=false, 
                            mov_avg=false, mov_avg_window=0, savein="", tstart=nothing, tend=nothing,
                            scatter=true,
                                markersize=10, 
                                marker=:circle,
                            line=true,
                                linewidth=2.5, 
                                linestyle=:solid,
                            fig = nothing, fpos = (1,1), disp=true
                            )
    
    # Field in plates.dat flag
    in_plates = (field=="Mob") || (field=="Vsurf")
    
    # Get encoding
    idxT, idxR, idxP = data_encoding(Dblock.timeheader, Dblock.rprofheader, Dblock.platesheader)
    in_plates ? (@assert haskey(idxP, field) "Field $field not found in plates header") : (@assert haskey(idxT, field) "Field $field not found in time header")

    # Automatic subsample
    (subsample==0 && !in_plates) && (subsample = max(1, Int(size(Dblock.timedata, 1) ÷ 500)))
    (subsample==0) && (subsample = max(1, Int(size(Dblock.platesdata, 1) ÷ 500)))

    # Vectors
    xvec = in_plates ? Dblock.platesdata[1:subsample:end, idxP["time"]] : Dblock.timedata[1:subsample:end,idxT["time"]]
    yvec = in_plates ? Dblock.platesdata[1:subsample:end, idxP[field]] : Dblock.timedata[1:subsample:end, idxT[field]]
    isnothing(tstart) && (tstart  = first(xvec))
    isnothing(tend) && (tend    = last(xvec))

    # Define mov_avg window automatically if not given
    (mov_avg && mov_avg_window == 0) && (mov_avg_window = max(3, Int(size(yvec, 1) ÷ 50)))

    # Plot
    isnothing(fig) && (fig = Figure(size = fsize))
    ax = Axis(fig[fpos[1], fpos[2]], xlabel = L"Time\;[Gyr]", ylabel = LTN[field], xlabelsize=xlabelsize, ylabelsize=ylabelsize, xgridvisible=xgrid, ygridvisible=ygrid, yreversed=yreversed,
               xticklabelsize=xticklabelsize, yticklabelsize=yticklabelsize, xticksize=xticksize, yticksize=yticksize, xlabelpadding=xlabelpadding, ylabelpadding=ylabelpadding)
    scatter && scatter!(ax, xvec, yvec, color=color, markersize=markersize, marker=marker, alpha = mov_avg ? 0.3 : 1.0)
    line && lines!(ax, xvec, yvec, color=color, linewidth=linewidth, linestyle=linestyle, alpha = mov_avg ? 0.3 : 1.0)
    mov_avg && lines!(ax, xvec, EasyFit.movavg(yvec, mov_avg_window).x, color=color, linewidth=linewidth, linestyle=:solid)
    tstartidx = findfirst(xvec .>= tstart); (isnothing(tstartidx)) && error("tstart=$tstart beyond data range")
    tendidx   = findlast(xvec .<= tend); (isnothing(tendidx)) && error("tend=$tend beyond data range")
    localmin, localmax = minimum(yvec[tstartidx:tendidx]), maximum(yvec[tstartidx:tendidx])
    Δ = 0.2(localmax - localmin)
    second_axis!(fig, fpos, xvec, yvec, field, ylabelsize, yticklabelsize, ylabelpadding, yreversed, true, (tstart, tend), (localmin-Δ, localmax+Δ))
    xlims!(ax, tstart, tend)
    (yreversed ? ylims!(ax, localmax+Δ, localmin-Δ) : ylims!(ax, localmin-Δ, localmax+Δ))
    disp && display(fig)
    (savein != "") && save(savein*".png", fig)
end

"""
    Plot the correlation between two time-dependent fields.
    
    \t Basic usage: \t field_vs_field(Dblock::DataBlock, fieldx::String, fieldy::String; kwargs...)

    Optional arguments (kwargs):

      • -- Canvas arguments --

        - fsize::Tuple{Int64,Int64} \t-->\t Figure size in pixels [default: (800, 600)]
        - xlabelsize::Int64 \t\t-->\t X-axis label font size [default: 25]
        - ylabelsize::Int64 \t\t-->\t Y-axis label font size [default: 25]
        - xticklabelsize::Int64 \t\t-->\t X-axis tick label font size [default: 17]
        - yticklabelsize::Int64 \t\t-->\t Y-axis tick label font size [default: 17]
        - xticksize::Int64 \t\t-->\t X-axis tick size [default: 10]
        - yticksize::Int64 \t\t-->\t Y-axis tick size [default: 10]
        - xlabelpadding::Int64 \t\t-->\t X-axis label padding [default: 15]
        - ylabelpadding::Int64 \t\t-->\t Y-axis label padding [default: 15]
        - savein::String \t\t-->\t Save figure in file [default: "" (no saving)]

      • -- Data arguments --

        - subsample_size::Int64 \t\t-->\t Number of points to plot [default: 1000]
        - tstart::Union{Float64, Nothing} \t-->\t Start time [default: nothing (first time)]
        - tend::Union{Float64, Nothing} \t-->\t End time [default: nothing (last time)]
        - color::Symbol \t\t\t-->\t Line/scatter color [default: :blue]
        - scatter::Bool \t\t\t-->\t Enables scatters [default: true]
            - markersize::Int64 \t\t-->\t Marker size [default: 10]
            - marker::Symbol \t\t-->\t Marker type [default: :circle]
        - line::Bool \t\t\t-->\t Enables lines [default: true]
            - linewidth::Float64 \t\t--
"""
function field_vs_field(Dblock::DataBlock, fieldx::String, fieldy::String; fsize=(800, 600), xlabelsize=25, ylabelsize=25,
                            xticklabelsize=17, yticklabelsize=17, xticksize=10, yticksize=10, xlabelpadding=15, ylabelpadding=15,
                            mov_avg=false, mov_avg_window=0,subsample_size=1000, color=:blue, tstart=nothing, tend=nothing,
                            scatter=true,
                                markersize=10,
                                marker=:circle,
                            line=true,
                                linewidth=2.5,
                                linestyle=:solid,
                            fig=nothing, fpos=(1,1), disp=true, savein=""
                            )

    # Field in plates.dat flag
    in_plates1 = (fieldx=="Mob") || (fieldx=="Vsurf")
    in_plates2 = (fieldy=="Mob") || (fieldy=="Vsurf")

    # Get encoding
    idxT, idxR, idxP = data_encoding(Dblock.timeheader, Dblock.rprofheader, Dblock.platesheader)
    @assert (haskey(idxT, fieldx) || haskey(idxP, fieldx)) "Field $fieldx not found in time header"
    @assert (haskey(idxT, fieldy) || haskey(idxP, fieldy)) "Field $fieldy not found in time header"

    # Automatic subsample
    xvec = in_plates1 ? Dblock.platesdata[:, idxP[fieldx]] : Dblock.timedata[:, idxT[fieldx]]
    yvec = in_plates2 ? Dblock.platesdata[:, idxP[fieldy]] : Dblock.timedata[:, idxT[fieldy]]
    tvec1 = in_plates1 ? Dblock.platesdata[:, idxP["time"]] : Dblock.timedata[:, idxT["time"]]
    tvec2 = in_plates2 ? Dblock.platesdata[:, idxP["time"]] : Dblock.timedata[:, idxT["time"]]

    # Time window
    isnothing(tstart) && (tstart = first(tvec1))
    isnothing(tend) && (tend = last(tvec1))
    tstartidx = findfirst(tvec1 .>= tstart); (isnothing(tstartidx)) && error("tstart=$tstart beyond data range")
    tendidx   = findlast(tvec1 .<= tend); (isnothing(tendidx)) && error("tend=$tend beyond data range")

    # Vectors
    itpx, itpy = interpolate((sort(tvec1),), xvec, Gridded(Linear())), interpolate((sort(tvec2),), yvec, Gridded(Linear()))
    subsample = max(1, Int(size(tvec1[tstartidx:tendidx], 1) ÷ subsample_size))
    tp = tvec1[tstartidx:subsample:tendidx]
    xp, yp = itpx(tp), itpy(tp)

    # Define mov_avg window automatically if not given
    (mov_avg && mov_avg_window == 0) && (mov_avg_window = max(3, Int(size(yp, 1) ÷ 50)))

    # Plot
    isnothing(fig) && (fig = Figure(size = fsize))
    ax = Axis(fig[fpos[1], fpos[2]], xlabel = LTN[fieldx], ylabel = LTN[fieldy], xlabelsize=xlabelsize, ylabelsize=ylabelsize,
               xticklabelsize=xticklabelsize, yticklabelsize=yticklabelsize, xticksize=xticksize, yticksize=yticksize,
               xlabelpadding=xlabelpadding, ylabelpadding=ylabelpadding)
    scatter && scatter!(ax, xp, yp, color=color, markersize=markersize, marker=marker, alpha = mov_avg ? 0.3 : 1.0)
    line && lines!(ax, xp, yp, color=color, linewidth=linewidth, linestyle=linestyle, alpha = mov_avg ? 0.3 : 1.0)
    mov_avg && lines!(ax, xp, EasyFit.movavg(yp, mov_avg_window).x, color=color, linewidth=linewidth, linestyle=:solid)

    # Set zoom
    localminx, localmaxx, localminy, localmaxy = 0.95minimum(xp), 1.05maximum(xp), 0.95minimum(yp), 1.05maximum(yp)
    xlims!(ax, localminx, localmaxx); ylims!(ax, localminy, localmaxy)

    # Display and save
    disp && display(fig)
    (savein != "") && save(savein*".png", fig)
end

"""
    Plot the evolution of a radial profile field over time.
    
    \t Basic usage: \t rprof_vs_field(Dblock::DataBlock, field::String; kwargs...)

    Optional arguments (kwargs):
        - fsize::Tuple{Int64,Int64} \t-->\t Figure size in pixels [default: (800, 600)]
        - Paxis::Bool \t\t\t-->\t Adds a pressure axis on the right [default: false]
        - cmap::Symbol \t\t\t-->\t Colormap [default: :vik100]
        - cmap_reverse::Bool \t\t-->\t Reverses colormap [default: false]
        - logscale::Bool \t\t\t-->\t Plots log₁₀ of the field [default: false]
        - tstart::Union{Float64, Nothing} \t-->\t Start time [default: nothing (first time)]
        - tend::Union{Float64, Nothing} \t-->\t End time [default: nothing (last time)]

"""
function rprof_vs_field(Dblock::DataBlock, field::String; fsize=(900, 600), cmap=:vik100, logscale=false, cmap_reverse=false, colorrange=(nothing, nothing), Paxis=true, tstart=nothing, tend=nothing,
                            fig=nothing, fpos=(1,1), disp=true, savein="", interpolate=false, np=500)

    !isnothing(fig) && (Paxis = false)

    # Get encoding
    idxT, idxR, idxP = data_encoding(Dblock.timeheader, Dblock.rprofheader, Dblock.platesheader)
    @assert haskey(idxR, field) "Field $field not found in rprof header"

    # Vectors
    yp = logscale ? log10.(Dblock.rprofdata[:,:,idxR[field]]) : Dblock.rprofdata[:,:,idxR[field]]
    isnothing(tstart) && (tstart = first(Dblock.rproftime))
    isnothing(tend) && (tend = last(Dblock.rproftime))

    # Interpolate to regular grid if requested
    if interpolate
        itp = Interpolations.interpolate((Dblock.rprofdata[:,1,idxR["r"]], Dblock.rproftime), yp, Gridded(Linear()))
        rgrid = LinRange(first(Dblock.rprofdata[:,1,idxR["r"]]), last(Dblock.rprofdata[:,1,idxR["r"]]), np)
        tgrid = LinRange(Dblock.rproftime[1], Dblock.rproftime[end], np)
        yp = [itp(r, t) for t in tgrid, r in rgrid]'
    end

    # Automatic colorrange if not given
    isnothing(colorrange[1]) && (colorrange = (minimum(yp), colorrange[2]))
    isnothing(colorrange[2]) && (colorrange = (colorrange[1], maximum(yp)))

    # Plot
    isnothing(fig) && (fig = Figure(size = fsize))
    ax = Axis(fig[Paxis ? 2 : fpos[1], fpos[2]], xlabel = L"Time\;[Gyr]", ylabel = L"Radius\;[km]", xlabelsize=25, ylabelsize=25, xticklabelsize=17, yticklabelsize=17, xticksize=10, yticksize=10, xlabelpadding=15, ylabelpadding=15)
    logscale && (field="log"*field)
    cmap_reverse && (cmap = Reverse(cmap))
    hm = heatmap!(ax, interpolate ? tgrid : Dblock.rproftime, interpolate ? 1e-3rgrid : 1e-3Dblock.rprofdata[:,idxR["r"],1], yp', colormap=cmap, interpolate=interpolate, colorrange=colorrange)
    Colorbar(fig[fpos[1], Paxis ? 1 : fpos[2]+1], hm; label = LRN[field], labelsize=25, ticklabelsize=15, vertical= Paxis ? false : true, labelpadding=13)
    Paxis && second_axis!(fig, fpos, Dblock.rproftime, 1e-3abs.(Dblock.rprofdata[:,idxR["r"],1] .- Dblock.rprofdata[end,idxR["r"],1]), "Pressure", 25, 17, 15, true, false, (tstart, tend), (nothing, nothing))
    xlims!(ax, tstart, tend)
    disp && display(fig)
    (savein != "") && save(savein*".png", fig)

end

"""
    Plot the mantle water content profile at a given time.
    
    \t Basic usage: \t mantle_water(Dblock::DataBlock, time::Float64; kwargs...)

    Optional arguments (kwargs):

      • -- Canvas arguments --

        - fsize::Tuple{Int64,Int64} \t-->\t Figure size in pixels [default: (900, 600)]
        - xlabelsize::Int64 \t\t-->\t X-axis label font size [default: 25]
        - ylabelsize::Int64 \t\t-->\t Y-axis label font size [default: 25]
        - xticklabelsize::Int64 \t\t-->\t X-axis tick label font size [default: 17]
        - yticklabelsize::Int64 \t\t-->\t Y-axis tick label font size [default: 17]
        - xticksize::Int64 \t\t-->\t X-axis tick size [default: 10]
        - yticksize::Int64 \t\t-->\t Y-axis tick size [default: 10]
        - xlabelpadding::Int64 \t\t-->\t X-axis label padding [default: 15]
        - ylabelpadding::Int64 \t\t-->\t Y-axis label padding [default: 15]
        - savein::String \t\t-->\t Save figure in file [default: "" (no saving)]
"""
function mantle_water(Dblock::DataBlock, time::Float64; fig=nothing, fpos=(1,1), fsize=(900,600), disp=true, savein="",
                        xlabelsize=25, ylabelsize=25, xticklabelsize=17, yticklabelsize=17, xticksize=10, yticksize=10, xlabelpadding=15, ylabelpadding=15)

    # Checks and indexing setup
    @assert Dblock.metadata.H₂O_tracked "Mantle water content tracking not enabled in simulation $(Dblock.metadata.Sname)"
    idxT, idxR, idxP = data_encoding(Dblock)
    mintime = min(time, Dblock.rproftime[end], Dblock.timedata[end,idxT["time"]])
    timeidx = findfirst(Dblock.rproftime .>= mintime)
    timeidxT = findfirst(Dblock.timedata[:,idxT["time"]] .>= mintime)

    # Initialize figure 
    isnothing(fig) && (fig = Figure(size = fsize))

    # Gather Data
    r = Dblock.rprofdata[:, timeidx, idxR["r"]]
    midpoint_r = 0.5(r[1:end-1] .+ r[2:end])
    midpoint_r² = 0.5(midpoint_r[1:end-1] .+ midpoint_r[2:end])
    ρ = Dblock.rprofdata[:, timeidx, idxR["rhomean"]]
    ∂ρ_∂r = diff(ρ)./diff(r)
    ∂²ρ_∂r² = diff(∂ρ_∂r)./diff(midpoint_r)
    water = Dblock.rprofdata[:, timeidx, idxR["Water"]]
    s = Dblock.rprofdata[:, timeidx, idxR["Wsol"]]
    ∂M = Dblock.rprofdata[:, timeidx, idxR["dM"]]
    H2Okg = water./om.*∂M.*1e-2 # wt% to fraction
    sH2Okg = s./om.*∂M.*1e-2

    # Find radial phase boundaries
    pv_ring, wad_ol, lith_ol = idx_ph_transitions(Dblock; timeidx=timeidx)

    # Ocean mass
    omass = Dblock.timedata[timeidxT,idxT["SurfOceanMass3D"]]/om

    # Divide by sectors
    lm_H2Okg = sum(H2Okg[1:pv_ring])
    tz_H2Okg = sum(H2Okg[pv_ring+1:wad_ol])
    um_H2Okg = sum(H2Okg[wad_ol+1:lith_ol])
    crust_H2Okg = sum(H2Okg[lith_ol+1:end])

    ax = Axis(fig[fpos[1], fpos[2]], xlabel = L"Radius\;[km]", ylabel = L"Cumulative\;H_2O\;mass\;[OM]", xlabelsize=xlabelsize, ylabelsize=ylabelsize,
               xticklabelsize=xticklabelsize, yticklabelsize=yticklabelsize, xticksize=xticksize, yticksize=yticksize,
               xlabelpadding=xlabelpadding, ylabelpadding=ylabelpadding)

    # Sector bands
    # -- LM
        band!(ax, 1e-3r[1:pv_ring], zeros(length(r[1:pv_ring])), lm_H2Okg*ones(length(r[1:pv_ring])), color=:orange, alpha=0.3)
        band!(ax, 1e-3r[pv_ring:end], zeros(length(r[pv_ring:end])), lm_H2Okg*ones(length(r[pv_ring:end])), color=:orange, alpha=0.1)
    # -- TZ
        band!(ax, 1e-3r[1:pv_ring], lm_H2Okg*ones(length(r[1:pv_ring])), (lm_H2Okg+tz_H2Okg)*ones(length(r[1:pv_ring])), color=:green, alpha=0.1)
        band!(ax, 1e-3r[pv_ring:wad_ol], lm_H2Okg*ones(length(r[pv_ring:wad_ol])), (lm_H2Okg+tz_H2Okg)*ones(length(r[pv_ring:wad_ol])), color=:green, alpha=0.3)
        band!(ax, 1e-3r[wad_ol:end], lm_H2Okg*ones(length(r[wad_ol:end])), (lm_H2Okg+tz_H2Okg)*ones(length(r[wad_ol:end])), color=:green, alpha=0.1)
    # -- UM
        band!(ax, 1e-3r[1:wad_ol], (lm_H2Okg+tz_H2Okg)*ones(length(r[1:wad_ol])), (lm_H2Okg+tz_H2Okg+um_H2Okg)*ones(length(r[1:wad_ol])), color=:blue, alpha=0.1)
        band!(ax, 1e-3r[wad_ol:lith_ol], (lm_H2Okg+tz_H2Okg)*ones(length(r[wad_ol:lith_ol])), (lm_H2Okg+tz_H2Okg+um_H2Okg)*ones(length(r[wad_ol:lith_ol])), color=:blue, alpha=0.3)
        band!(ax, 1e-3r[lith_ol:end], (lm_H2Okg+tz_H2Okg)*ones(length(r[lith_ol:end])), (lm_H2Okg+tz_H2Okg+um_H2Okg)*ones(length(r[lith_ol:end])), color=:blue, alpha=0.1)
    # -- Crust
        band!(ax, 1e-3r[1:lith_ol], (lm_H2Okg+tz_H2Okg+um_H2Okg)*ones(length(r[1:lith_ol])), (lm_H2Okg+tz_H2Okg+um_H2Okg+crust_H2Okg)*ones(length(r[1:lith_ol])), color=:brown, alpha=0.1)
        band!(ax, 1e-3r[lith_ol:end], (lm_H2Okg+tz_H2Okg+um_H2Okg)*ones(length(r[lith_ol:end])), (lm_H2Okg+tz_H2Okg+um_H2Okg+crust_H2Okg)*ones(length(r[lith_ol:end])), color=:brown, alpha=0.3)

    # Text
    f = 0.5
    text!(ax, 200, f*(lm_H2Okg), text="$(round(lm_H2Okg, digits=3)) OM (Lower mantle)", color=:black, align = (:left, :center), fontsize=20)
    text!(ax, 200, lm_H2Okg + f*tz_H2Okg, text="$(round(tz_H2Okg, digits=3)) OM (MTZ)", color=:black, align = (:left, :center), fontsize=20)
    text!(ax, 200, lm_H2Okg + tz_H2Okg + f*um_H2Okg, text="$(round(um_H2Okg, digits=3)) OM (Upper mantle)", color=:black, align = (:left, :center), fontsize=20)
    text!(ax, 200, lm_H2Okg + tz_H2Okg + um_H2Okg + f*crust_H2Okg, text="$(round(crust_H2Okg, digits=3)) OM (Crust)", color=:black, align = (:left, :center), fontsize=20)

    # Lines
    lines!(ax, 1e-3r, cumsum(H2Okg), color=:black, linewidth=3.0, alpha=1.0)
    # lines!(ax, 1e-3r, cumsum(sH2Okg), color=:blue, linestyle=:dash, linewidth=2.5, alpha=0.3, label = "Mantle water capacity")
    scatter!(ax, 0.0, 0.0, color=:black, alpha=0.0, label = "Ocean mass at $(round(time, digits=2)) Gyr : $(round(omass, digits=2)) OM")
    # scatter!(ax, 0.0, 0.0, color=:black, alpha=0.0, label = "Parameter | Integrated Total H₂O mass at $(round(time, digits=2)) Gyr : $(round(Dblock.metadata.totH₂O/om, digits=2)) | $(round(sum(H2Okg)+omass, digits=2)) OM")
    scatter!(ax, 0.0, 0.0, color=:black, alpha=0.0, label = "Integrated mantle H₂O mass at $(round(time, digits=2)) Gyr : $(round(sum(H2Okg), digits=2)) OM")
    axislegend(ax, position=:rb, framevisible=true, fontsize=15, padding=10, rowgap=10)

    # Display and save
    disp && display(fig)
    (savein != "") && save(savein*".png", fig)
end

"""
    Plot a snapshot of a 2D StagYY at a given time.
    
    \t Basic usage: \t snapshot(Dblock::DataBlock, stime::Float64, field::Union{String, Vector{String}}; kwargs...)

    Optional arguments (kwargs):

      • -- Canvas arguments --

        - fsize::Tuple{Int64,Int64} \t-->\t Figure size in pixels [default: (800, 800)]
        - fpos::Tuple{Int64,Int64} \t-->\t Figure position in a grid layout [default: (1, 1)]
        - cmaparray::Union{Nothing, Vector{ColorScheme}} \t-->\t Colormap array for multiple fields [default: nothing]
        - logscale::Bool \t\t\t-->\t Plots log₁₀ of the field [default: false]
"""
function snapshot(Dblock, stime, field; fig=nothing, fpos=(1,1), fsize=(800,800), cmaparray=nothing, disp=true, logscale=false)

    # Read selected VTk file
    pvd_fname = joinpath(Dblock.metadata.outdir, Dblock.metadata.sroot*".pvd")
    # Check file existance
    @assert isfile(pvd_fname) "VTK PVD file $pvd_fname not found. Make sure to convert binary output into readable VTK using StagYY utilities."
    # Read
    vtkname, vtktime = pick_VTK_file_at_time(min(stime, Dblock.metadata.tend), pvd_fname)
    vtk = readVTK(vtkname)

    # Get colormap for fields
    typeof(field) == String && (field = [field])
    isnothing(cmaparray) ? (cmap = [CTF[f] for f in field]) : (cmap = cmaparray)

    # Expand cell colors into point colors
    c = reshape(vtk[field[1]], Dblock.metadata.nx, Dblock.metadata.ny, Dblock.metadata.nz)  # Cell value array
    p = zeros(Float32, Dblock.metadata.nx+1, Dblock.metadata.ny+1, Dblock.metadata.nz+1)    # Empty point value array
    Threads.@threads for k in 1:Dblock.metadata.nz+1
        for j in 1:Dblock.metadata.ny+1
            @inbounds for i in 1:Dblock.metadata.nx+1
                ii = min(i, Dblock.metadata.nx)
                jj = min(j, Dblock.metadata.ny)
                kk = min(k, Dblock.metadata.nz)
                p[i,j,k] = c[ii,jj,kk]
            end
        end
    end

    # Initialize figure 
    isnothing(fig) && (fig = Figure(size = fsize))
    ax = Axis(fig[fpos[1], fpos[2]], xgridvisible=false, ygridvisible=false)
    hidespines!(ax); hidexdecorations!(ax); hideydecorations!(ax);
    plot!(ax, vtk["Points"][1,:], vtk["Points"][2,:], color= logscale ? log10.(vec(p)) : vec(p), colormap=cmap[1], markersize=8)
    
    # Time print
    text!(ax, 0.0, 0.0, text="Time = $(round(1e-3vtktime, digits=2)) Gyr", color=:black, align = (:center, :center), fontsize=20)

    # display
    disp && display(fig)
end

function IOplot(Dblock; fig=nothing, fpos=(1,1), fsize=(800,800), disp=true)

    # Checks
    @assert Dblock.metadata.H₂O_tracked "Mantle water content tracking not enabled in simulation $(Dblock.metadata.Sname)"

    # Get encoding
    idxT, idxR, idxP = data_encoding(Dblock.timeheader, Dblock.rprofheader, Dblock.platesheader)

    # Vectors
    time = Dblock.timedata[:,idxT["time"]]
    rawingas, rawoutgas = Dblock.timedata[:,idxT["IngassedH2O"]], Dblock.timedata[:,idxT["OutgassedH2O"]]
    mavg_ingas, mavg_outgas = EasyFit.movavg(rawingas, 10).x, EasyFit.movavg(rawoutgas, 10).x

    # Initialise figure
    isnothing(fig) && (fig = Figure(size = fsize))
    ax = Axis(fig[fpos[1], fpos[2]], xlabel = L"Time\;[\mathrm{Gyr}]", ylabel = L"H_2O\;mass\;[\mathrm{OM}]", xlabelsize=25, ylabelsize=25, xticklabelsize=17, yticklabelsize=17, xticksize=10, yticksize=10, xlabelpadding=15, ylabelpadding=15)
    lines!(ax, time, rawingas, color=:blue, linewidth=2.5, linestyle=:solid, label="Ingassed H₂O")
    display(fig)


end

# =================================================================================
# =======  Water Storage Capacity (sᴴ²ᴼ) and Oxygen fugacity profile (fO₂)  =======
# =================================================================================

"""
    Solve for mantle water storage capacity (sᴴ²ᴼ) and oxygen fugacity (fO₂) over P-T-X space.

    \t Basic usage: \t solve_sH2O_fO2(nP::Int64, nT::Int64; kwargs...)

    Optional arguments (kwargs):

      • -- General arguments --

        - s::Bool \t\t\t-->\t Solve for sᴴ²ᴼ [default: true]
        - fO2::Bool \t\t\t-->\t Solve for fO₂ [default: true]
        - DBswitchP::Float64 \t\t-->\t Pressure (GPa) to switch between upper mantle and transition zone databases [default: 7.0]
        - plt::Bool \t\t\t-->\t Plot results [default: false]
        - disp_prog::Bool \t\t-->\t Display progress bars [default: true]
        - DHMS::Bool \t\t\t-->\t Explore DHMS paths for transition zone calculations [default: true]
        - Clist::Vector{String} \t-->\t List of oxides in the system [default: ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "Cr2O3", "O", "H2O"]]
        - XB::Vector{Float64} \t\t-->\t Bulk composition in wt% for sᴴ²ᴼ calculations [default: Upper mantle composition]
        - XH::Vector{Float64} \t\t-->\t Bulk composition in wt% for fO₂ calculations [default: Mid-ocean ridge basalt composition]
        - Prange::Tuple{Float64, Float64} \t-->\t Pressure range in GPa [default: (0.1, 130.0)]
        - Trange::Tuple{Float64, Float64} \t-->\t Temperature range in K [default: (500.0, 4000.0)]
        - verbose::Bool \t\t-->\t Enable verbose output [default: true]

      • -- Plot arguments -- (plot_sᴴ²ᴼ(...) call)

        - cmap::Symbol \t\t\t-->\t Colormap for plots [default: :vik100]
        - interp::Bool \t\t\t-->\t Interpolate results for smoother plots [default: false]
        - cmap_reverse::Bool \t\t-->\t Reverse colormap [default: false]
        - logscale::Bool \t\t\t-->\t Plot log₁₀ of sᴴ²ᴼ or fO₂ [default: true]
        - phase_out::Vector{String} \t-->\t List of phases to exclude from calculations [default: ["chl"]]
"""
function solve_sH2O_fO2(nP::Int64, nT::Int64;
                        s=true, fO2=true, DBswitchP=7.0, plt=false, disp_prog=true, DHMS=true,
                        Clist=["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "Cr2O3", "O", "H2O"],
                        XB=[49.33, 15.31, 10.82, 7.41, 10.33, 0.19, 2.53, 1.46, 0.0, 0.0, 100.0],
                        XH=[45.5, 2.59, 4.05, 35.22, 7.26, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0],
                        Prange=(0.1, 130.0), Trange=(500.0, 4000.0), verbose=true,
                        cmap=:vik100, interp=false, cmap_reverse=false, logscale=true, phase_out=["chl"],
                        test_path=false
                        )

    # Checks
    (!s && !fO2) && return
    @assert length(Clist) == length(XB) "Length of Clist and XB must match."
    @assert length(Clist) == length(XH) "Length of Clist and XH must match."

    # Vetorization size
    nPnT = nP*nT

    if s
        # Memory allocations
        # --- Axis Vectors
            Pum, Tum = LinRange(Prange[1], DBswitchP, nP), LinRange(Trange[1], 2000., nT)
            Ptz, Ttz = LinRange(DBswitchP, 25., nP), LinRange(1000., Trange[2], nT)
            Plm, Tlm = LinRange(25., Prange[2], nP), LinRange(1000., Trange[2], nT)
        # --- Others
            Pv, Tv = zeros(Float64, 2nPnT), zeros(Float64, 2nPnT)
            Xv = vcat(map(Vector, eachrow(repeat(XH', outer=nPnT))), map(Vector, eachrow(repeat(XB', outer=nPnT))))
            tnP, tnP05 = Int64(ceil(1.1max(nP, nT))), Int64(ceil(0.5max(nP, nT)))
        # --- Maps
            um, tz, lm = zeros(Float64, nP*nT, 2), zeros(Float64, nP*nT, 2), zeros(Float64, nP*nT, 2) # later reshaped as 'T, P, reshape(um, nP, :)'

        # Upper mantle mesh vectorization (MAGEMin parallelization)
            mesh_vectorization!(Pum, Tum, nP, nT, Pv, Tv)
        # Minimizer call + assembly
            verbose && println("Calculating upper mantle sᴴ²ᴼ...")
            data    = Initialize_MAGEMin("um", verbose=false, buffer="aH2O");
            rm_list = remove_phases(phase_out, "um")
            outHB   = multi_point_minimization(10Pv, Tv.-273.15, data, X=Xv, Xoxides=Clist, sys_in="wt", name_solvus=true, B=ones(length(Pv)), progressbar=disp_prog, rm_list=rm_list) # kbar and K
            Finalize_MAGEMin(data);
            sᴴ²ᴼ_assembler!(um, outHB, nPnT)
            um = cat(reshape(um[:,1], nP, nT), reshape(um[:,2], nP, nT), dims=3)

        # Mineral-bound sᴴ²ᴼ assembly
            verbose && println("Calculating Mineral-bound sᴴ²ᴼ curves...")
            min_s = min_sᴴ²ᴼ_assembler(tnP);

        # Transition zone mesh vectorization
            mesh_vectorization!(Ptz, Ttz, nP, nT, Pv, Tv)
        # Minimizer call + assembly
            verbose && println("Calculating transition zone sᴴ²ᴼ...")
            data    = Initialize_MAGEMin("mtl", verbose=false);
            outHB   = multi_point_minimization(10Pv, Tv.-273.15, data, X=Xv, Xoxides=Clist, sys_in="wt", name_solvus=true, progressbar=disp_prog) # kbar and K
            Finalize_MAGEMin(data);

            DHMS && verbose && println("Exploring DHMS paths...")
            ppaths, paths  =  DHMS ? path_solve(XH, Clist, phase_out, (@view Pv[1:nPnT]), (@view Tv[1:nPnT]), DBswitchP, (@view outHB[1:nPnT]), test_path; npaths=max(50,tnP05), ns=max(100,tnP), Pend=Prange[2]) : (0.0, 0.0)
            ∫sᴴ²ᴼ!(tz, (@view Pv[1:nPnT]), (@view Tv[1:nPnT]), min_s, outHB, DHMS, ppaths, paths, nPnT)
            tz = cat(reshape(tz[:,1], nP, nT), reshape(tz[:,2], nP, nT), dims=3)
    
        # Lower mantle mesh vectorization
            mesh_vectorization!(Plm, Tlm, nP, nT, Pv, Tv)
        # Minimizer call + assembly
            verbose && println("Calculating lower mantle sᴴ²ᴼ...")
            data    = Initialize_MAGEMin("sb21", verbose=false);
            outHB   = multi_point_minimization(10Pv, Tv.-273.15, data, X=Xv, Xoxides=Clist, sys_in="wt", name_solvus=true, progressbar=disp_prog) # kbar and K
            ∫sᴴ²ᴼ!(lm, (@view Pv[1:nPnT]), (@view Tv[1:nPnT]), min_s, outHB, DHMS, ppaths, paths, nPnT)
            lm = cat(reshape(lm[:,1], nP, nT), reshape(lm[:,2], nP, nT), dims=3)
            Finalize_MAGEMin(data);
        # Return structure
        smap = sᴴ²ᴼ( um, tz, lm, Pum, Tum, Ptz, Ttz, Plm, Tlm )

        # Export
        fmap = zeros(Float64, nP, nT) # Dummy
        write_output(smap, fmap, s=s, fO2=fO2)

        # Plot
        plt && plot_sᴴ²ᴼ(smap; cmap = cmap, interp = interp, cmap_reverse = cmap_reverse, logscale = logscale)

        return smap

        end

end

# ======================
# ===== Minimizer ======
# ======================

"""
    Plot phase diagram minimization maps for a given mantle sector and composition.
    
    \t Basic usage: \t minmap(sector::String, em::String; kwargs...)

    Required arguments:

      • sector::String \t-->\t Mantle sector to model ('um', 'tz', 'lm')

      • em::String \t\t-->\t Earth mantle composition ('XH', 'XB', 'Custom')

    Optional arguments (kwargs):

        - nP::Int64 \t\t-->\t Number of pressure points [default: 50]
        - nT::Int64 \t\t-->\t Number of temperature points [default: 50]
        - Clist::Array{String,1} \t-->\t List of oxides in composition [default: ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "Cr2O3", "O", "H2O"]]
        - XB::Array{Float64,1} \t-->\t StagYY's enriched endmember [default: [49.33, 15.31, 10.82, 7.41, 10.33, 0.19, 2.53, 1.46, 0.0, 0.0, 100.0]]
        - XH::Array{Float64,1} \t-->\t Bulk composition for HARZBURGITE (wt%) [default: [45.5, 2.59, 4.05, 35.22, 7.26, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0]]
        - ccomp::Array{Float64,1} \t-->\t Custom bulk composition (wt%) [default: zeros(Float64, length(Clist))]
        - DBswitchP::Float64 \t-->\t Pressure (GPa) of upper mantle to transition zone switch [default: 7.0]
        - interp::Bool \t\t-->\t Interpolates heatmap [default: false]
        - cmap::Symbol \t\t-->\t Colormap [default: :BuPu]
        - Prange ::Tuple{Float64,Float64} \t-->\t Pressure range (GPa) [default: (0.1, 130.0)]
        - Trange ::Tuple{Float64,Float64} \t-->\t Temperature range (K) [default: (500.0, 4000.0)]
        - ncols::Int64 \t\t-->\t Number of columns in figure [default: 4]
        - savein::String \t-->\t Save figure in file [default: "" (no saving)]
        - phase_out::Array{String,1} \t-->\t List of phases to exclude from minimization [default: ["chl"]]
        - H2Osat::Bool \t-->\t Considers H2O saturation in upper mantle calculations [default: true]
"""
function minmap(sector, em::String; nP=50, nT=50,
                Clist=["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "Cr2O3", "O", "H2O"],
                XB=[49.33, 15.31, 10.82, 7.41, 10.33, 0.19, 2.53, 1.46, 0.0, 0.0, 0.0],
                XH=[45.5, 2.59, 4.05, 35.22, 7.26, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                ccomp=zeros(Float64, length(Clist)),DBswitchP=7.0, interp=false, cmap=:BuPu,
                Prange=(0.1, 130.0), Trange=(500.0, 4000.0), ncols=4, savein="", phase_out=["chl"], H2Osat=true
                )

    # Check
    @assert em in ["XH", "XB", "Custom"] "Invalid Earth mantle: $em. Choose from 'XH', 'XB', 'Custom'."
    # Composition
    X = em=="XB" ? XB : em=="XH" ? XH : ccomp
    H2Osat && (X[end] = 100.0) # Set H2O to 100 for saturation calculations
    # Vectors
    Pum, Tum = LinRange(Prange[1], DBswitchP, nP), LinRange(Trange[1], 2000., nT)
    Ptz, Ttz = LinRange(DBswitchP, 25., nP), LinRange(1200., Trange[2], nT)
    Plm, Tlm = LinRange(25., Prange[2], nP), LinRange(1200., Trange[2], nT)
    Pv, Tv = zeros(Float64, nP*nT), zeros(Float64, nP*nT)
    Xv = map(Vector, eachrow(repeat(X', outer=length(Pv))))
    Pv .= repeat(sector=="um" ? Pum : sector=="tz" ? Ptz : Plm, outer=nT); Tv .= repeat(sector=="um" ? Tum : sector=="tz" ? Ttz : Tlm, inner=nP)
    # Minimizer
    if sector == "um"
        data = H2Osat ? Initialize_MAGEMin("um", verbose=false, buffer="aH2O") : Initialize_MAGEMin("um", verbose=false);
        rm_list = remove_phases(phase_out, "um")
        out = multi_point_minimization(10Pv, Tv.-273.15, data, X=Xv, Xoxides=Clist, B=ones(length(Pv)), sys_in="wt", name_solvus=true, rm_list=rm_list)
    else
        data = Initialize_MAGEMin(sector=="tz" ? "mtl" : "sb21", verbose=false);
        out = multi_point_minimization(10Pv, Tv.-273.15, data, X=Xv, Xoxides=Clist, sys_in="wt", name_solvus=true)
    end
    Finalize_MAGEMin(data);
    # Phase iterator
    observed_phases = Dict{String, Array{Float64, 2}}()
    for n in eachindex(out)
        # Vectorized to mesh indexing
        i, j = get_ip_it_nrows(nP, n)
        # Phase iterator
        for (idx, phase) in enumerate(out[n].ph)
            # Initialise phase if not present
            !haskey(observed_phases, phase) && (observed_phases[phase] = zeros(Float64, nP, nT))
            # Assign solution phase value
            observed_phases[phase][i, j] = out[n].ph_frac_wt[idx]
        end
    end
    # plot
    fig = Figure(size=(1000,1000)); i=1; j=1
    xlabsz, ylabsz, titlesz, xticklabsz, yticklabsz, xticksz, yticksz = 20, 20, 22, 16, 16, 12, 12
    for (ph, map) in observed_phases
        ax = Axis(fig[i, j], title=ph, ylabel=L"Pressure\;[GPa]", xlabel=L"Temperature\;[K]", xlabelsize=xlabsz, ylabelsize=ylabsz,
                    xticklabelsize=xticklabsz, yticklabelsize=yticklabsz, xticksize=xticksz, yticksize=yticksz, titlesize=titlesz, yreversed=true)
        hm = heatmap!(ax, sector=="um" ? Tum : sector=="tz" ? Ttz : Tlm, sector=="um" ? Pum : sector=="tz" ? Ptz : Plm, map', colormap=cmap, interpolate=interp)
        Colorbar(fig[i, j+1], hm; labelsize=20, ticklabelsize=15, vertical=true, labelpadding=13)
        j+=2; (j==2ncols+1) && (i+=1; j=1)
    end
    !isempty(savein) && CairoMakie.save("./" * savein * ".png", fig)
    display(fig)
end


"""
    Solve a single P-T point for given composition and return the minimization output.
    
    \t Basic usage: \t solve_point(P::Float64, T::Float64, em::String; kwargs...)

    Optional arguments (kwargs):

        - Clist::Array{String,1} \t-->\t List of oxides in the composition [default: ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "Cr2O3", "O", "H2O"]]
        - XB::Array{Float64,1} \t\t-->\t StagYY's default enriched endmember [default: [49.33, 15.31, 10.82, 7.41, 10.33, 0.19, 2.53, 1.46, 0.0, 0.0, 100.0]]
        - XH::Array{Float64,1} \t\t-->\t StagYY's default depleted endmember [default: [45.5, 2.59, 4.05, 35.22, 7.26, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0]]
        - ccomp::Array{Float64,1} \t-->\t Custom composition in wt% [default: zeros(Float64, length(Clist))]
        - DBswitchP::Float64 \t\t-->\t Pressure (GPa) at which to switch from upper mantle to transition zone database [default: 7.0 GPa]
        - phase_out::Array{String,1} \t-->\t List of phases to remove from the minimization [default: ["chl"]]
        - H2Osat::Bool \t\t-->\t Considers H2O saturation in upper mantle calculations [default: true]
"""
function solve_point(P, T, em;
                    Clist=["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "Cr2O3", "O", "H2O"],
                    XB=[49.33, 15.31, 10.82, 7.41, 10.33, 0.19, 2.53, 1.46, 0.0, 0.0, 100.0],
                    XH=[45.5, 2.59, 4.05, 35.22, 7.26, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0],
                    ccomp=zeros(Float64, length(Clist)),DBswitchP=7.0,phase_out=["chl"],H2Osat=true
                    )

    # Checks
    @assert length(XB) == length(Clist) "Length of Clist and XB must match."
    @assert length(XH) == length(Clist) "Length of Clist and XH must match."
    @assert em in ["XH", "XB", "Custom"] "Endmember must be either: XB, XH or Custom"
    @assert (em!="Custom" || (em=="Custom" && sum(ccomp)!=0.0)) "Please provide a custom composition if chosen"

    # Minimize
    X = em=="XB" ? XB : em=="XH" ? XH : ccomp
    H2Osat && (X[end] = 100.0) # Set H2O to 100 for saturation calculations
    if P <= DBswitchP
    data = H2Osat ? Initialize_MAGEMin("um", verbose=false, buffer="aH2O") : Initialize_MAGEMin("um", verbose=false);
    rm_list = remove_phases(phase_out, "um")
    out = single_point_minimization(10P, T-273.15, data, X=X, Xoxides=Clist, B=1.0, sys_in="wt", name_solvus=true, rm_list=rm_list)
    else
        data = Initialize_MAGEMin(P<=25. ? "mtl" : "sb21", verbose=false);
        out = single_point_minimization(10P, T-273.15, data, X=X, Xoxides=Clist, sys_in="wt", name_solvus=true)
    end
    Finalize_MAGEMin(data);

    return out
end