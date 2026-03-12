# time.dat field --> String
LTN = Dict(
    "time"  => L"Evolution\;time\;[\mathrm{Gyr}]",
    "F_top" => L"Surface\;heat\;flux\;[\mathrm{mW/m^2}]", "F_top_norm" => L"Normalized\;surface\;heat\;flux",
    "F_bot" => L"CMB\;heat\;flux\;[\mathrm{mW/m^2}]", "F_bot_norm" => L"Normalized\;CMB\;heat\;flux",
    "Tmean" => L"Mean\;mantle\;temperature\;[\mathrm{K}]", "Tmean_norm" => L"Normalized\;mean\;mantle\;temperature",
    "Vrms"  => L"RMS\;velocity\;[\mathrm{cm/yr}]", "Vrms_norm"  => L"Normalized\;RMS\;velocity",
    "eta_amean" => L"Mean\;viscosity\;[\mathrm{Pa\cdot s}]", "eta_amean_norm" => L"Normalized\;mean\;viscosity",
    "eta_mean"  => L"Mean\;viscosity\;[\mathrm{Pa\cdot s}]", "eta_mean_norm"  => L"Normalized\;mean\;viscosity",
    "erupta" => L"Recycled\;material\;[\mathrm{Mantle\;masses}]", "erupta_norm" => L"Normalized\;recycled\;material",
    "erupt_rate" => L"Eruption\;rate\;[\mathrm{kg/s}]", "erupt_rate_norm" => L"Normalized\;eruption\;rate",
    "Tcmb" => L"CMB\;temperature\;[\mathrm{K}]", "Tcmb_norm" => L"Normalized\;CMB\;temperature",
    "Tsurf" => L"Surface\;temperature\;[\mathrm{K}]", "Tsurf_norm" => L"Normalized\;surface\;temperature",
    "Tpotl" => L"Potential\;temperature\;[\mathrm{K}]", "Tpotl_norm" => L"Normalized\;potential\;temperature",
    "cH2O_mean" => L"Mean\;mantle\;H_2O\;content\;[\mathrm{wt%}]", "cH2O_mean_norm" => L"Normalized\;mean\;mantle\;H_2O\;content",
    "IngassedH2O" => L"Total\;ingassed\;H_2O\;[\mathrm{kg}]", "IngassedH2O_norm" => L"Normalized\;total\;ingassed\;H_2O",
    "OutgassedH2O" => L"Total\;outgassed\;H_2O\;[\mathrm{kg}]", "OutgassedH2O_norm" => L"Normalized\;total\;outgassed\;H_2O",
    "SaturationOutgasH2O" => L"Saturation\;outgassed\;H_2O\;[\mathrm{kg}]", "SaturationOutgassH2O_norm" => L"Normalized\;saturation\;outgassed\;H_2O",
    "EruptedH2O" => L"Erupted\;H_2O\;[\mathrm{kg}]", "EruptedH2O_norm" => L"Normalized\;erupted\;H_2O",
    "eH2O/e" => L"erupted\;H_2O\;/\;erupta\;[\mathrm{wt%}]", "eH2O/e_norm" => L"Normalized\;erupted\;H_2O\;/\;erupta",
    "I/O_H2O" => L"Ingassed\;/\;Outgassed\;H_2O", "I/O_H2O_norm" => L"Normalized\;ingassed\;/\;outgassed\;H_2O",
    "SurfOceanMass3D" => L"Surface\;ocean\;mass\;[\mathrm{OM}]", "SurfOceanMass3D_norm" => L"Normalized\;surface\;ocean\;mass",
    "TotH2OMassOceanPlusMantle" => L"Total\;system\;H_2O\;[\mathrm{kg}]", "TotH2OMassOceanPlusMantle_norm" => L"Normalized\;total\;system\;H_2O",
    "Mob" => L"Lithosphere\;mobility\;[\mathrm{v_{surf}/v_{rms}}]", "Mob_norm" => L"Normalized\;lithosphere\;mobility",
    "Vsurf" => L"Surface\;velocity\;[\mathrm{cm/yr}]", "Vsurf_norm" => L"Normalized\;surface\;velocity",

)

# rprof.dat field --> String
LRN = Dict(
    "Tmean" => L"Mean\;mantle\;temperature\;[\mathrm{K}]", "logTmean" => L"Mean\;mantle\;temperature\;[log_{10}K]",
    "Vrms"  => L"RMS\;velocity\;[\mathrm{cm/yr}]", "logvrms" => L"RMS\;velocity\;[log_{10}(\mathrm{cm/yr})]",
    "eta_amean" => L"Mean\;viscosity\;[\mathrm{Pa\cdot s}]", "logeta_amean" => L"Mean\;viscosity\;[log_{10}(\mathrm{Pa\cdot s})]",
    "eta_mean"  => L"Mean\;viscosity\;[\mathrm{Pa\cdot s}]", "logeta_mean"  => L"Mean\;viscosity\;[log_{10}(\mathrm{Pa\cdot s})]",
    "Water" => L"H_2O\;content\;[\mathrm{wt%}]", "logWater" => L"H_2O\;content\;[log_{10}(\mathrm{wt%})]",
    "boundH2O" => L"Bound\;H_2O\;[\mathrm{wt%}]", "logbH2O" => L"Bound\;H_2O\;[log_{10}(\mathrm{wt%})]",
    "freeH2O" => L"Free\;H_2O\;[\mathrm{wt%}]", "logfH2O" => L"Free\;H_2O\;[log_{10}(\mathrm{wt%})]",
    "Wsol" => L"H_2O\;storage\;capacity\;[\mathrm{wt%}]", "logWsol" => L"H_2O\;storage\;capacity\;[log_{10}(\mathrm{wt%})]",
    "satH2O" => L"H_2O\;saturation\;[\mathrm{%}]", "logsatH2O" => L"H_2O\;saturation\;[log_{10}(\mathrm{%})]",
    "fO2" => L"fO_2\;[log_{10}FMQ]", "logfO2" => L"fO_2\;[log_{10}FMQ]",
    "H2Ofree" => L"Free\;H_2O\;[\mathrm{wt%}]", "logH2Ofree" => L"Free\;H_2O\;[log_{10}(\mathrm{wt%})]",
    "rhomean" => L"Mean\;density\;[\mathrm{kg/m^3}]", "logrhomean" => L"Mean\;density\;[log_{10}(\mathrm{kg/m^3})]",
    "bsmean" => L"Mean\;basalt\;fraction", "logbsmean" => L"Mean\;basalt\;fraction\;[\mathrm{log_{10}}]",
    "etalog" => L"Viscosity\;[\mathrm{Pa\cdot s}]", "logetalog" => L"Viscosity\;[log_{10}(\mathrm{Pa\cdot s})]",
    "H2Odarcy" => L"H_2O\;darcy\;velocity\;[\mathrm{cm/yr}]", "logH2Odarcy" => L"H_2O\;darcy\;velocity\;[log_{10}(\mathrm{cm/yr})]",
    "fmeltmean" => L"Mean\;melt\;fraction", "logfmeltmean" => L"Mean\;melt\;fraction\;[log_{10}]",
    "Tmean_um" => L"Mean\;mantle\;temperature\;UM\;[\mathrm{K}]", "logTmean_um" => L"Mean\;mantle\;temperature\;UM\;[log_{10}K]",
    "Tmean_tz" => L"Mean\;mantle\;temperature\;TZ\;[\mathrm{K}]", "logTmean_tz" => L"Mean\;mantle\;temperature\;TZ\;[log_{10}K]",
    "Tmean_lm" => L"Mean\;mantle\;temperature\;LM\;[\mathrm{K}]", "logTmean_lm" => L"Mean\;mantle\;temperature\;LM\;[log_{10}K]",
    "Tmean_crust" => L"Mean\;mantle\;temperature\;Crust\;[\mathrm{K}]", "logTmean_crust" => L"Mean\;mantle\;temperature\;Crust\;[log_{10}K]",
    "Vrms_um"  => L"RMS\;velocity\;UM\;[\mathrm{cm/yr}]", "logvrms_um" => L"RMS\;velocity\;UM\;[log_{10}(\mathrm{cm/yr})]",
    "Vrms_tz"  => L"RMS\;velocity\;TZ\;[\mathrm{cm/yr}]", "logvrms_tz" => L"RMS\;velocity\;TZ\;[log_{10}(\mathrm{cm/yr})]",
    "Vrms_lm"  => L"RMS\;velocity\;LM\;[\mathrm{cm/yr}]", "logvrms_lm" => L"RMS\;velocity\;LM\;[log_{10}(\mathrm{cm/yr})]",
    "Vrms_crust"  => L"RMS\;velocity\;Crust\;[\mathrm{cm/yr}]", "logvrms_crust" => L"RMS\;velocity\;Crust\;[log_{10}(\mathrm{cm/yr})]",
    "fmeltmean_um" => L"Mean\;melt\;fraction\;UM", "logfmeltmean_um" => L"Mean\;melt\;fraction\;UM\;[log_{10}]",
    "fmeltmean_tz" => L"Mean\;melt\;fraction\;TZ", "logfmeltmean_tz" => L"Mean\;melt\;fraction\;TZ\;[log_{10}]",
    "fmeltmean_lm" => L"Mean\;melt\;fraction\;LM", "logfmeltmean_lm" => L"Mean\;melt\;fraction\;LM\;[log_{10}]",
    "fmeltmean_crust" => L"Mean\;melt\;fraction\;Crust", "logfmeltmean_crust" => L"Mean\;melt\;fraction\;Crust\;[log_{10}]",
    "etalog_um" => L"Viscosity\;UM\;[\mathrm{Pa\cdot s}]", "logetalog_um" => L"Viscosity\;UM\;[log_{10}(\mathrm{Pa\cdot s})]",
    "etalog_tz" => L"Viscosity\;TZ\;[\mathrm{Pa\cdot s}]", "logetalog_tz" => L"Viscosity\;TZ\;[log_{10}(\mathrm{Pa\cdot s})]",
    "etalog_lm" => L"Viscosity\;LM\;[\mathrm{Pa\cdot s}]", "logetalog_lm" => L"Viscosity\;LM\;[log_{10}(\mathrm{Pa\cdot s})]",
    "etalog_crust" => L"Viscosity\;Crust\;[\mathrm{Pa\cdot s}]", "logetalog_crust" => L"Viscosity\;Crust\;[log_{10}(\mathrm{Pa\cdot s})]",
    "Water_um" => L"H_2O\;content\;UM\;[\mathrm{wt%}]", "logWater_um" => L"H_2O\;content\;UM\;[log_{10}(\mathrm{wt%})]",
    "Water_tz" => L"H_2O\;content\;TZ\;[\mathrm{wt%}]", "logWater_tz" => L"H_2O\;content\;TZ\;[log_{10}(\mathrm{wt%})]",
    "Water_lm" => L"H_2O\;content\;LM\;[\mathrm{wt%}]", "logWater_lm" => L"H_2O\;content\;LM\;[log_{10}(\mathrm{wt%})]",
    "Water_crust" => L"H_2O\;content\;Crust\;[\mathrm{wt%}]", "logWater_crust" => L"H_2O\;content\;Crust\;[log_{10}(\mathrm{wt%})]",
    "boundH2O_um" => L"Bound\;H_2O\;UM\;[\mathrm{wt%}]", "logbH2O_um" => L"Bound\;H_2O\;UM\;[log_{10}(\mathrm{wt%})]",
    "boundH2O_tz" => L"Bound\;H_2O\;TZ\;[\mathrm{wt%}]", "logbH2O_tz" => L"Bound\;H_2O\;TZ\;[\mathrm{wt%}]", 
    "boundH2O_lm" => L"Bound\;H_2O\;LM\;[\mathrm{wt%}]", "logbH2O_lm" => L"Bound\;H_2O\;LM\;[\mathrm{wt%}]", 
    "boundH2O_crust" => L"Bound\;H_2O\;Crust\;[\mathrm{wt%}]", "logbH2O_crust" => L"Bound\;H_2O\;Crust\;[\mathrm{wt%}]",
    "freeH2O_um" => L"Free\;H_2O\;UM\;[\mathrm{wt%}]", "logfH2O_um" => L"Free\;H_2O\;UM\;[log_{10}(\mathrm{wt%})]",
    "freeH2O_tz" => L"Free\;H_2O\;TZ\;[\mathrm{wt%}]", "logfH2O_tz" => L"Free\;H_2O\;TZ\;[log_{10}(\mathrm{wt%})]",
    "freeH2O_lm" => L"Free\;H_2O\;LM\;[\mathrm{wt%}]", "logfH2O_lm" => L"Free\;H_2O\;LM\;[\mathrm{wt%}]", 
    "freeH2O_crust" => L"Free\;H_2O\;Crust\;[\mathrm{wt%}]", "logfH2O_crust" => L"Free\;H_2O\;Crust\;[\mathrm{wt%}]",
    "Wsol_um" => L"H_2O\;storage\;capacity\;UM\;[\mathrm{wt%}]", "logWsol_um" => L"H_2O\;storage\;capacity\;UM\;[log_{10}(\mathrm{wt%})]",
    "Wsol_tz" => L"H_2O\;storage\;capacity\;TZ\;[\mathrm{wt%}]", "logWsol_tz" => L"H_2O\;storage\;capacity\;TZ\;[log_{10}(\mathrm{wt%})]",
    "Wsol_lm" => L"H_2O\;storage\;capacity\;LM\;[\mathrm{wt%}]", "logWsol_lm" => L"H_2O\;storage\;capacity\;LM\;[log_{10}(\mathrm{wt%})]",
    "Wsol_crust" => L"H_2O\;storage\;capacity\;Crust\;[\mathrm{wt%}]", "logWsol_crust" => L"H_2O\;storage\;capacity\;Crust\;[log_{10}(\mathrm{wt%})]",
    "satH2O_um" => L"H_2O\;saturation\;UM\;[\mathrm{%}]", "logsatH2O_um" => L"H_2O\;saturation\;UM\;[log_{10}(\mathrm{%})]",
    "satH2O_tz" => L"H_2O\;saturation\;TZ\;[\mathrm{%}]", "logsatH2O_tz" => L"H_2O\;saturation\;TZ\;[log_{10}(\mathrm{%})]",
    "satH2O_lm" => L"H_2O\;saturation\;LM\;[\mathrm{%}]", "logsatH2O_lm" => L"H_2O\;saturation\;LM\;[log_{10}(\mathrm{%})]",
    "satH2O_crust" => L"H_2O\;saturation\;Crust\;[\mathrm{%}]", "logsatH2O_crust" => L"H_2O\;saturation\;Crust\;[log_{10}(\mathrm{%})]",
)

# Colormap to field dictionary
CTF = Dict(
    "Temperature" => :vik100,
    "T: Residual"  => :vik100,
    "Viscosity"   => :davos,
    "C:Water"       => :Blues,
    "Water_solubility" => Reverse(:grays),
    "fO2"         => :plasma,
    "Density"     => :turbo,
    "C:Basalt"      => Reverse(:RdGy),
    "Velocity"    => :viridis,
    "Melt_fraction" => :default,
    "dmelt/dt"      => :default,
)

# ======================================
# ====  General Post-Processing API ====
# ======================================

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
function time_vs_field(Dblock, field::String; fsize=(800, 600), subsample=0, color=:blue, xlabelsize=25, ylabelsize=25, xgrid=true, ygrid=true,
                            xticklabelsize=17, yticklabelsize=17, xticksize=10, yticksize=10, xlabelpadding=15, ylabelpadding=15, yreversed=false, 
                            mov_avg=false, mov_avg_window=0, savein="", tstart=nothing, tend=nothing,
                            scatter=true,
                                markersize=10, 
                                marker=:circle,
                            line=true,
                                linewidth=2.5, 
                                linestyle=:solid,
                            fig = nothing, fpos = (1,1), disp=true, logscale=false,
                            )
    
    if Dblock isa Vector{DataBlock}
        time_vs_field_multiple(Dblock, field; fsize=fsize, subsample=subsample, xlabelsize=xlabelsize, ylabelsize=ylabelsize, xgrid=xgrid, ygrid=ygrid,
                            xticklabelsize=xticklabelsize, yticklabelsize=yticklabelsize, xticksize=xticksize, yticksize=yticksize, xlabelpadding=xlabelpadding, ylabelpadding=ylabelpadding, yreversed=yreversed, 
                            mov_avg=mov_avg, mov_avg_window=mov_avg_window, savein=savein, tstart=tstart, tend=tend,
                            scatter=scatter,
                                markersize=markersize,
                                marker=marker,
                            line=line,
                                linewidth=linewidth, 
                                linestyle=linestyle,
                            fig = fig, fpos = fpos, disp=disp,
                            logscale=logscale,
                            clrsmap = :vik100)
        return
    end

    # Field to use
    field2use = replace(field, "_lm"=>"", "_tz"=>"", "_um"=>"", "_crust"=>"")

    # Field in plates.dat flag
    in_plates = (field2use=="Mob") || (field2use=="Vsurf") && !contains(field, "_lm") && !contains(field, "_tz") && !contains(field, "_um") && !contains(field, "_crust")
    in_time  = (field2use in Dblock.timeheader) && !contains(field, "_lm") && !contains(field, "_tz") && !contains(field, "_um") && !contains(field, "_crust")
    in_rprof = (field2use in Dblock.rprofheader) && !in_time && !in_plates
    @assert in_plates || in_time || in_rprof "Field $field2use not found in any header"
    
    # Get encoding
    idxT, idxR, idxP = data_encoding(Dblock.timeheader, Dblock.rprofheader, Dblock.platesheader)

    # Find radial phase boundaries
    pv_ring, wad_ol, lith_ol = idx_ph_transitions(Dblock; timeidx=1)
    rprof_range = if contains(field, "_lm")
        1:pv_ring
    elseif contains(field, "_tz")
        pv_ring+1:wad_ol
    elseif contains(field, "_um")
        wad_ol+1:lith_ol
    elseif contains(field, "_crust")
        lith_ol+1:size(Dblock.rprofdata, 1)
    else
        1:size(Dblock.rprofdata, 1)
    end

    # Automatic subsample
    (subsample==0 && !in_plates && !in_rprof) && (subsample = max(1, Int(size(Dblock.timedata, 1) ÷ 500)))
    (subsample==0 && in_rprof) && (subsample = max(1, Int(length(Dblock.rproftime) ÷ 500)))
    (subsample==0) && (subsample = max(1, Int(size(Dblock.platesdata, 1) ÷ 500)))

    # Vectors
    xvec = in_plates ? Dblock.platesdata[1:subsample:end, idxP["time"]] : in_time ? Dblock.timedata[1:subsample:end,idxT["time"]] : Dblock.rproftime[1:subsample:end]
    yvec = in_plates ? Dblock.platesdata[1:subsample:end, idxP[field2use]] : in_time ? Dblock.timedata[1:subsample:end, idxT[field2use]] : 
            vec(map(col -> mean(skipmissing(col)), eachcol(replace(Dblock.rprofdata[rprof_range, 1:subsample:end, idxR[field2use]], NaN=>missing))))
    isnothing(tstart) && (tstart  = first(xvec))
    isnothing(tend) && (tend    = last(xvec))

    # Check for log
    logscale && (yvec .= max.(yvec, 1e-8))

    # Define mov_avg window automatically if not given
    (mov_avg && mov_avg_window == 0) && (mov_avg_window = max(3, Int(size(yvec, 1) ÷ 50)))

    # Plot
    isnothing(fig) && (fig = Figure(size = fsize))
    ax = Axis(fig[fpos[1], fpos[2]], xlabel = L"Time\;[Gyr]", ylabel = in_rprof ? LRN[field] : LTN[field], xlabelsize=xlabelsize, ylabelsize=ylabelsize, xgridvisible=xgrid, ygridvisible=ygrid, yreversed=yreversed,
               xticklabelsize=xticklabelsize, yticklabelsize=yticklabelsize, xticksize=xticksize, yticksize=yticksize, xlabelpadding=xlabelpadding, ylabelpadding=ylabelpadding, yscale= logscale ? log10 : identity)
    scatter && scatter!(ax, xvec, yvec, color=color, markersize=markersize, marker=marker, alpha = mov_avg ? 0.2 : 1.0)
    line && lines!(ax, xvec, yvec, color=color, linewidth=linewidth, linestyle=linestyle, alpha = mov_avg ? 0.2 : 1.0)
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

function time_vs_field_multiple(Dblock, field::String; fsize=(800, 600), subsample=0, xlabelsize=25, ylabelsize=25, xgrid=true, ygrid=true,
                            xticklabelsize=17, yticklabelsize=17, xticksize=10, yticksize=10, xlabelpadding=15, ylabelpadding=15, yreversed=false, 
                            mov_avg=false, mov_avg_window=0, savein="", tstart=nothing, tend=nothing,
                            scatter=true,
                                markersize=10, 
                                marker=:circle,
                            line=true,
                                linewidth=2.5, 
                                linestyle=:solid,
                            fig = nothing, fpos = (1,1), disp=true,
                            clrsmap = :vik100, logscale=false,
                            )

    # field to use
    field2use = replace(field, "_lm"=>"", "_tz"=>"", "_um"=>"", "_crust"=>"")
    
    # Field in plates.dat flag
    in_plates   = (field2use=="Mob") || (field2use=="Vsurf") && !contains(field, "_lm") && !contains(field, "_tz") && !contains(field, "_um") && !contains(field, "_crust")
    in_time     = (field2use in Dblock[1].timeheader) && !contains(field, "_lm") && !contains(field, "_tz") && !contains(field, "_um") && !contains(field, "_crust")
    in_rprof    = (field2use in Dblock[1].rprofheader) && !in_time && !in_plates
    @assert in_plates || in_time || in_rprof "Field $field2use not found in any header"

    # Plot
    isnothing(fig) && (fig = Figure(size = fsize))
    ax = Axis(fig[fpos[1], fpos[2]], xlabel = L"Time\;[Gyr]", ylabel = in_rprof ? LRN[field] : LTN[field], xlabelsize=xlabelsize, ylabelsize=ylabelsize, xgridvisible=xgrid, ygridvisible=ygrid, yreversed=yreversed,
               xticklabelsize=xticklabelsize, yticklabelsize=yticklabelsize, xticksize=xticksize, yticksize=yticksize, xlabelpadding=xlabelpadding, ylabelpadding=ylabelpadding, yscale=logscale ? log10 : identity)
    clrs = cpalette(clrsmap, length(Dblock))
    localmin, localmax = Inf, -Inf
    isnothing(tstart) && (ts = Inf)
    isnothing(tend) && (te = -Inf)
    for (b, blck) in enumerate(Dblock)
        # Get encoding
        idxT, idxR, idxP = data_encoding(blck.timeheader, blck.rprofheader, blck.platesheader)

        # Find radial phase boundaries
        pv_ring, wad_ol, lith_ol = idx_ph_transitions(blck; timeidx=1)
        rprof_range = if contains(field, "_lm")
            1:pv_ring
        elseif contains(field, "_tz")
            pv_ring+1:wad_ol
        elseif contains(field, "_um")
            wad_ol+1:lith_ol
        elseif contains(field, "_crust")
            lith_ol+1:size(blck.rprofdata, 1)
        else
            1:size(blck.rprofdata, 1)
        end

        # Automatic subsample
        (subsample==0 && !in_plates && !in_rprof) && (subsample = max(1, Int(size(blck.timedata, 1) ÷ 500)))
        (subsample==0 && in_plates) && (subsample = max(1, Int(size(blck.platesdata, 1) ÷ 500)))
        (subsample==0 && in_rprof) && (subsample = max(1, Int(length(blck.rproftime) ÷ 500)))

        # Vectors
        xvec = in_plates ? blck.platesdata[1:subsample:end, idxP["time"]] : in_time ? blck.timedata[1:subsample:end,idxT["time"]] : blck.rproftime[1:subsample:end]
        yvec = in_plates ? blck.platesdata[1:subsample:end, idxP[field2use]] : in_time ? blck.timedata[1:subsample:end, idxT[field2use]] : 
                    vec(map(col -> mean(skipmissing(col)), eachcol(replace(blck.rprofdata[rprof_range, 1:subsample:end, idxR[field2use]], NaN=>missing))))
        
        # Time axis cuts
        isnothing(tstart) && (ts  = min(first(xvec), ts); idxt1 = findfirst(xvec .>= ts))
        isnothing(tend) && (te    = max(last(xvec), te); idxt2 = findlast(xvec .<= te))
        !isnothing(tstart) && (idxt1 = findfirst(xvec .>= max(tstart, first(xvec))))
        !isnothing(tend) && (idxt2 = findlast(xvec .<= min(tend, last(xvec))))

        # Check for log
        logscale && (yvec .= max.(yvec, 1e-8))

        # edges
        localmin = min(localmin, minimum(yvec[idxt1:idxt2]))
        localmax = max(localmax, maximum(yvec[idxt1:idxt2]))

        # Define mov_avg window automatically if not given
        (mov_avg && mov_avg_window == 0) && (mov_avg_window = max(3, Int(size(yvec, 1) ÷ 50)))

        scatter && scatter!(ax, xvec, yvec, color=clrs[b], markersize=markersize, marker=marker, alpha = mov_avg ? 0.3 : 1.0)
        scatter!(ax, xvec[1], yvec[1], color=clrs[b], markersize=markersize, marker=marker, label=blck.metadata.name)
        line && lines!(ax, xvec, yvec, color=clrs[b], linewidth=linewidth, linestyle=linestyle, alpha = mov_avg ? 0.3 : 1.0)
        mov_avg && lines!(ax, xvec, EasyFit.movavg(yvec, mov_avg_window).x, color=clrs[b], linewidth=linewidth, linestyle=:solid)
    end

    isnothing(tstart) && (tstart = ts)
    isnothing(tend) && (tend     = te)
    Δ = 0.05(localmax - localmin)
    !logscale && (yreversed ? ylims!(ax, localmax+Δ, localmin-Δ) : ylims!(ax, localmin-Δ, localmax+Δ))
    xlims!(ax, tstart, tend)
    axislegend(ax; position = :lt, unique=true)
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
function field_vs_field(Dblock, fieldx::String, fieldy::String; fsize=(800, 600), xlabelsize=25, ylabelsize=25,
                            xticklabelsize=17, yticklabelsize=17, xticksize=10, yticksize=10, xlabelpadding=15, ylabelpadding=15,
                            mov_avg=false, mov_avg_window=0,subsample_size=1000, color=:blue, tstart=nothing, tend=nothing, logx=false, logy=false,
                            scatter=true,
                                markersize=10,
                                marker=:circle,
                            line=true,
                                linewidth=2.5,
                                linestyle=:solid,
                            fig=nothing, fpos=(1,1), disp=true, savein=""
                            )

    if Dblock isa Vector{DataBlock}
        field_vs_field_multiple(Dblock, fieldx, fieldy; fsize=fsize, xlabelsize=xlabelsize, ylabelsize=ylabelsize,
                            xticklabelsize=xticklabelsize, yticklabelsize=yticklabelsize, xticksize=xticksize, yticksize=yticksize, xlabelpadding=xlabelpadding, ylabelpadding=ylabelpadding,
                            mov_avg=mov_avg, mov_avg_window=mov_avg_window,subsample_size=subsample_size, color=color, tstart=tstart, tend=tend, logx=logx, logy=logy,
                            scatter=scatter,
                                markersize=markersize,
                                marker=marker,
                            line=line,
                                linewidth=linewidth,
                                linestyle=linestyle,
                            fig=fig, fpos=fpos, disp=disp, savein=savein
                            )
        return
    end

    # field 2 use
    fieldx2use = replace(fieldx, "_lm"=>"", "_tz"=>"", "_um"=>"", "_crust"=>"")
    fieldy2use = replace(fieldy, "_lm"=>"", "_tz"=>"", "_um"=>"", "_crust"=>"")

    # Field in plates.dat flag
    in_plates1 = ((fieldx2use=="Mob") || (fieldx2use=="Vsurf")) && !contains(fieldx, "_lm") && !contains(fieldx, "_tz") && !contains(fieldx, "_um") && !contains(fieldx, "_crust")
    in_plates2 = ((fieldy2use=="Mob") || (fieldy2use=="Vsurf")) && !contains(fieldy, "_lm") && !contains(fieldy, "_tz") && !contains(fieldy, "_um") && !contains(fieldy, "_crust")
    in_time1  = (fieldx2use in Dblock.timeheader) && !in_plates1 && !contains(fieldx, "_lm") && !contains(fieldx, "_tz") && !contains(fieldx, "_um") && !contains(fieldx, "_crust")
    in_time2  = (fieldy2use in Dblock.timeheader) && !in_plates2 && !contains(fieldy, "_lm") && !contains(fieldy, "_tz") && !contains(fieldy, "_um") && !contains(fieldy, "_crust")
    in_rprof1 = (fieldx2use in Dblock.rprofheader) && !in_plates1 && !in_time1
    in_rprof2 = (fieldy2use in Dblock.rprofheader) && !in_plates2 && !in_time2
    @assert in_plates1 || in_time1 || in_rprof1 "Field $fieldx2use not found in any header"
    @assert in_plates2 || in_time2 || in_rprof2 "Field $fieldy2use not found in any header"

    # Get encoding
    idxT, idxR, idxP = data_encoding(Dblock.timeheader, Dblock.rprofheader, Dblock.platesheader)

    # Find radial phase boundaries
    pv_ring, wad_ol, lith_ol = idx_ph_transitions(Dblock; timeidx=1)
    rprof_rangex = if contains(fieldx, "_lm")
        1:pv_ring
    elseif contains(fieldx, "_tz")
        pv_ring+1:wad_ol
    elseif contains(fieldx, "_um")
        wad_ol+1:lith_ol
    elseif contains(fieldx, "_crust")
        lith_ol+1:size(Dblock.rprofdata, 1)
    else
        1:size(Dblock.rprofdata, 1)
    end
    rprof_rangey = if contains(fieldy, "_lm")
        1:pv_ring
    elseif contains(fieldy, "_tz")
        pv_ring+1:wad_ol
    elseif contains(fieldy, "_um")
        wad_ol+1:lith_ol
    elseif contains(fieldy, "_crust")
        lith_ol+1:size(Dblock.rprofdata, 1)
    else
        1:size(Dblock.rprofdata, 1)
    end

    # Automatic subsample
    xvec = in_plates1 ? Dblock.platesdata[:, idxP[fieldx2use]] : in_rprof1 ? vec(mean(Dblock.rprofdata[rprof_rangex, :, idxR[fieldx2use]], dims=1)) : Dblock.timedata[:, idxT[fieldx2use]]
    yvec = in_plates2 ? Dblock.platesdata[:, idxP[fieldy2use]] : in_rprof2 ? vec(mean(Dblock.rprofdata[rprof_rangey, :, idxR[fieldy2use]], dims=1)) : Dblock.timedata[:, idxT[fieldy2use]]
    tvec1 = in_plates1 ? Dblock.platesdata[:, idxP["time"]] : in_rprof1 ? Dblock.rproftime : Dblock.timedata[:, idxT["time"]]
    tvec2 = in_plates2 ? Dblock.platesdata[:, idxP["time"]] : in_rprof2 ? Dblock.rproftime : Dblock.timedata[:, idxT["time"]]

    # OM
    (fieldx == "SurfOceanMass3D") && (xvec ./= om)
    (fieldy == "SurfOceanMass3D") && (yvec ./= om)

    # Time window
    tstart1, tstart2, tend1, tend2 = tstart, tstart, tend, tend
    if isnothing(tstart)
        tstart1 = first(tvec1)
        tstart2 = first(tvec2)
    end
    if isnothing(tend)
        tend1 = last(tvec1)
        tend2 = last(tvec2)
    end
    tstartidx1 = findfirst(tvec1 .>= tstart1); (isnothing(tstartidx1)) && error("tstart=$tstart1 beyond data range")
    tstartidx2 = findfirst(tvec2 .>= tstart2); (isnothing(tstartidx2)) && error("tstart=$tstart2 beyond data range")
    tendidx1   = findlast(tvec1 .<= tend1); (isnothing(tendidx1)) && error("tend=$tend1 beyond data range")
    tendidx2   = findlast(tvec2 .<= tend2); (isnothing(tendidx2)) && error("tend=$tend2 beyond data range")

    # Vectors
    itpx, itpy = interpolate((sort(tvec1),), xvec, Gridded(Linear())), interpolate((sort(tvec2),), yvec, Gridded(Linear()))
    subsample1 = max(1, Int(size(tvec1[tstartidx1:tendidx1], 1) ÷ min(subsample_size, length(tvec1), length(tvec2))))
    subsample2 = max(1, Int(size(tvec2[tstartidx2:tendidx2], 1) ÷ min(subsample_size, length(tvec1), length(tvec2))))
    tp1 = tvec1[tstartidx1:subsample1:tendidx1]
    tp2 = tvec2[tstartidx2:subsample2:tendidx2]
    tp = length(tp1) <= length(tp2) ? tp1 : tp2
    xp, yp = itpx(tp), itpy(tp)

    # Define mov_avg window automatically if not given
    (mov_avg && mov_avg_window == 0) && (mov_avg_window = max(3, Int(size(yp, 1) ÷ 50)))

    # Plot
    isnothing(fig) && (fig = Figure(size = fsize))
    ax = Axis(fig[fpos[1], fpos[2]], xlabel = in_rprof1 ? LRN[fieldx] : LTN[fieldx], ylabel = in_rprof2 ? LRN[fieldy] : LTN[fieldy], xlabelsize=xlabelsize, ylabelsize=ylabelsize,
               xticklabelsize=xticklabelsize, yticklabelsize=yticklabelsize, xticksize=xticksize, yticksize=yticksize, xscale=logx ? log10 : identity, yscale=logy ? log10 : identity,
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

function field_vs_field_multiple(Dblock, fieldx::String, fieldy::String; fsize=(800, 600), xlabelsize=25, ylabelsize=25,
                            xticklabelsize=17, yticklabelsize=17, xticksize=10, yticksize=10, xlabelpadding=15, ylabelpadding=15,
                            mov_avg=false, mov_avg_window=0,subsample_size=1000, color=:blue, tstart=nothing, tend=nothing, logx=false, logy=false,
                            scatter=true,
                                markersize=10,
                                marker=:circle,
                            line=true,
                                linewidth=2.5,
                                linestyle=:solid,
                            fig=nothing, fpos=(1,1), disp=true, savein=""
                            )

    # field 2 use
    fieldx2use = replace(fieldx, "_lm"=>"", "_tz"=>"", "_um"=>"", "_crust"=>"")
    fieldy2use = replace(fieldy, "_lm"=>"", "_tz"=>"", "_um"=>"", "_crust"=>"")

    # Field in plates.dat flag
    in_plates1 = ((fieldx2use=="Mob") || (fieldx2use=="Vsurf")) && !contains(fieldx, "_lm") && !contains(fieldx, "_tz") && !contains(fieldx, "_um") && !contains(fieldx, "_crust")
    in_plates2 = ((fieldy2use=="Mob") || (fieldy2use=="Vsurf")) && !contains(fieldy, "_lm") && !contains(fieldy, "_tz") && !contains(fieldy, "_um") && !contains(fieldy, "_crust")
    in_time1  = (fieldx2use in Dblock[1].timeheader) && !in_plates1 && !contains(fieldx, "_lm") && !contains(fieldx, "_tz") && !contains(fieldx, "_um") && !contains(fieldx, "_crust")
    in_time2  = (fieldy2use in Dblock[1].timeheader) && !in_plates2 && !contains(fieldy, "_lm") && !contains(fieldy, "_tz") && !contains(fieldy, "_um") && !contains(fieldy, "_crust")
    in_rprof1 = (fieldx2use in Dblock[1].rprofheader) && !in_plates1 && !in_time1
    in_rprof2 = (fieldy2use in Dblock[1].rprofheader) && !in_plates2 && !in_time2
    @assert in_plates1 || in_time1 || in_rprof1 "Field $fieldx2use not found in any header"
    @assert in_plates2 || in_time2 || in_rprof2 "Field $fieldy2use not found in any header"

    # Plot
    isnothing(fig) && (fig = Figure(size = fsize))
    color = cpalette(:vik100, length(Dblock))
    ax = Axis(fig[fpos[1], fpos[2]], xlabel = in_rprof1 ? LRN[fieldx] : LTN[fieldx], ylabel = in_rprof2 ? LRN[fieldy] : LTN[fieldy], xlabelsize=xlabelsize, ylabelsize=ylabelsize,
               xticklabelsize=xticklabelsize, yticklabelsize=yticklabelsize, xticksize=xticksize, yticksize=yticksize, xscale=logx ? log10 : identity, yscale=logy ? log10 : identity,
               xlabelpadding=xlabelpadding, ylabelpadding=ylabelpadding)

    for (b, blck) in enumerate(Dblock)

        # Get encoding
        idxT, idxR, idxP = data_encoding(blck.timeheader, blck.rprofheader, blck.platesheader)

        # Find radial phase boundaries
        pv_ring, wad_ol, lith_ol = idx_ph_transitions(blck; timeidx=1)
        rprof_rangex = if contains(fieldx, "_lm")
            1:pv_ring
        elseif contains(fieldx, "_tz")
            pv_ring+1:wad_ol
        elseif contains(fieldx, "_um")
            wad_ol+1:lith_ol
        elseif contains(fieldx, "_crust")
            lith_ol+1:size(blck.rprofdata, 1)
        else
            1:size(blck.rprofdata, 1)
        end
        rprof_rangey = if contains(fieldy, "_lm")
            1:pv_ring
        elseif contains(fieldy, "_tz")
            pv_ring+1:wad_ol
        elseif contains(fieldy, "_um")
            wad_ol+1:lith_ol
        elseif contains(fieldy, "_crust")
            lith_ol+1:size(blck.rprofdata, 1)
        else
            1:size(blck.rprofdata, 1)
        end

        # Automatic subsample
        xvec = in_plates1 ? blck.platesdata[:, idxP[fieldx2use]] : in_rprof1 ? vec(mean(blck.rprofdata[rprof_rangex, :, idxR[fieldx2use]], dims=1)) : blck.timedata[:, idxT[fieldx2use]]
        yvec = in_plates2 ? blck.platesdata[:, idxP[fieldy2use]] : in_rprof2 ? vec(mean(blck.rprofdata[rprof_rangey, :, idxR[fieldy2use]], dims=1)) : blck.timedata[:, idxT[fieldy2use]]
        tvec1 = in_plates1 ? blck.platesdata[:, idxP["time"]] : in_rprof1 ? blck.rproftime : blck.timedata[:, idxT["time"]]
        tvec2 = in_plates2 ? blck.platesdata[:, idxP["time"]] : in_rprof2 ? blck.rproftime : blck.timedata[:, idxT["time"]]

        # OM
        (fieldx == "SurfOceanMass3D") && (xvec ./= om)
        (fieldy == "SurfOceanMass3D") && (yvec ./= om)

        # Time window
        tstart1, tstart2, tend1, tend2 = tstart, tstart, tend, tend
        if isnothing(tstart)
            tstart1 = first(tvec1)
            tstart2 = first(tvec2)
        end
        if isnothing(tend)
            tend1 = last(tvec1)
            tend2 = last(tvec2)
        end
        tstartidx1 = findfirst(tvec1 .>= tstart1); (isnothing(tstartidx1)) && error("tstart=$tstart1 beyond data range")
        tstartidx2 = findfirst(tvec2 .>= tstart2); (isnothing(tstartidx2)) && error("tstart=$tstart2 beyond data range")
        tendidx1   = findlast(tvec1 .<= tend1); (isnothing(tendidx1)) && error("tend=$tend1 beyond data range")
        tendidx2   = findlast(tvec2 .<= tend2); (isnothing(tendidx2)) && error("tend=$tend2 beyond data range")

        # Vectors
        itpx, itpy = interpolate((sort(tvec1),), xvec, Gridded(Linear())), interpolate((sort(tvec2),), yvec, Gridded(Linear()))
        subsample1 = max(1, Int(size(tvec1[tstartidx1:tendidx1], 1) ÷ min(subsample_size, length(tvec1), length(tvec2))))
        subsample2 = max(1, Int(size(tvec2[tstartidx2:tendidx2], 1) ÷ min(subsample_size, length(tvec1), length(tvec2))))
        tp1 = tvec1[tstartidx1:subsample1:tendidx1]
        tp2 = tvec2[tstartidx2:subsample2:tendidx2]
        tp = length(tp1) <= length(tp2) ? tp1 : tp2
        xp, yp = itpx(tp), itpy(tp)

        # Define mov_avg window automatically if not given
        (mov_avg && mov_avg_window == 0) && (mov_avg_window = max(3, Int(size(yp, 1) ÷ 50)))

        # Plot
        scatter && scatter!(ax, xp, yp, color=color[b], markersize=markersize, marker=marker, alpha = mov_avg ? 0.3 : 1.0, label=blck.metadata.name)
        line && lines!(ax, xp, yp, color=color[b], linewidth=linewidth, linestyle=linestyle, alpha = mov_avg ? 0.3 : 1.0)
        mov_avg && lines!(ax, xp, EasyFit.movavg(yp, mov_avg_window).x, color=color[b], linewidth=linewidth, linestyle=:solid)

    end

    # Legend
    axislegend(ax; position = :ct, unique=true, orientation=:horizontal)

    # # Set zoom
    # localminx, localmaxx, localminy, localmaxy = 0.95minimum(xp), 1.05maximum(xp), 0.95minimum(yp), 1.05maximum(yp)
    # xlims!(ax, localminx, localmaxx); ylims!(ax, localminy, localmaxy)

    # Display and save
    disp && display(fig)
    (savein != "") && save(savein*".png", fig)
end

function ta_field_vs_field(Dblocks, fieldx::String, fieldy::String; fsize=(800, 600), xlabelsize=25, ylabelsize=25,
                            xticklabelsize=17, yticklabelsize=17, xticksize=10, yticksize=10, xlabelpadding=15, ylabelpadding=15,
                            mov_avg=false, mov_avg_window=0,subsample_size=1000, color=:blue, tstart=nothing, tend=nothing,
                            scatter=true,
                                markersize=10,
                                strokewidth=1.0,
                                marker=:circle,
                            line=true,
                                linewidth=2.5,
                                linestyle=:solid,
                            fig=nothing, fpos=(1,1), disp=true, savein="",
                            setlabs=nothing,
                            )

    if (Dblocks isa Vector{Vector{DataBlock}})
        color = (color isa Symbol) ? [:blue, :red, :green, :orange, :purple] : color
        ta_field_vs_field_multiple(Dblocks, fieldx, fieldy; fsize=fsize, xlabelsize=xlabelsize, ylabelsize=ylabelsize,
                            xticklabelsize=xticklabelsize, yticklabelsize=yticklabelsize, xticksize=xticksize, yticksize=yticksize,
                            xlabelpadding=xlabelpadding, ylabelpadding=ylabelpadding,
                            mov_avg=mov_avg, mov_avg_window=mov_avg_window, subsample_size=subsample_size, color=color,
                            tstart=tstart, tend=tend,
                            scatter=scatter,
                                markersize=markersize,
                                strokewidth=1.0,
                                marker=marker,
                            line=line,
                                linewidth=linewidth,
                                linestyle=linestyle,
                            fig=fig, fpos=fpos, disp=disp, savein=savein, setlabs=setlabs)
        return
    elseif Dblocks isa DataBlock
        error("if you want to time-average a single datablock, use t_avg() instead.")
    end

    # Field to use
    fieldx2use = replace(fieldx, "_lm"=>"", "_tz"=>"", "_um"=>"", "_crust"=>"", "_norm"=>"")
    fieldy2use = replace(fieldy, "_lm"=>"", "_tz"=>"", "_um"=>"", "_crust"=>"", "_norm"=>"")

    # Field in plates.dat flag
    in_plates1 = ((fieldx2use=="Mob") || (fieldx2use=="Vsurf")) && !contains(fieldx, "_lm") && !contains(fieldx, "_tz") && !contains(fieldx, "_um") && !contains(fieldx, "_crust")
    in_plates2 = ((fieldy2use=="Mob") || (fieldy2use=="Vsurf")) && !contains(fieldy, "_lm") && !contains(fieldy, "_tz") && !contains(fieldy, "_um") && !contains(fieldy, "_crust")
    in_time1  = (fieldx2use in Dblocks[1].timeheader) && !contains(fieldx, "_lm") && !contains(fieldx, "_tz") && !contains(fieldx, "_um") && !contains(fieldx, "_crust")
    in_time2  = (fieldy2use in Dblocks[1].timeheader) && !contains(fieldy, "_lm") && !contains(fieldy, "_tz") && !contains(fieldy, "_um") && !contains(fieldy, "_crust")
    in_rprof1 = (fieldx2use in Dblocks[1].rprofheader) && !in_time1 && !in_plates1
    in_rprof2 = (fieldy2use in Dblocks[1].rprofheader) && !in_time2 && !in_plates2
    @assert in_plates1 || in_time1 || in_rprof1 "Field $fieldx2use not found in any header"
    @assert in_plates2 || in_time2 || in_rprof2 "Field $fieldy2use not found in any header"

    tax = zeros(Float64, length(Dblocks))
    tay = zeros(Float64, length(Dblocks))

    # Set up normalization if wanted
    normx = contains(fieldx,"_norm")
    normy = contains(fieldy,"_norm")
    normx && (basex = 0.0)
    normy && (basey = 0.0)
    
    for (b, blck) in enumerate(Dblocks)

        # Get encoding
        idxT, idxR, idxP = data_encoding(blck.timeheader, blck.rprofheader, blck.platesheader)

        # Find radial phase boundaries
        pv_ring, wad_ol, lith_ol = idx_ph_transitions(blck; timeidx=1)
        rprof_rangex = if contains(fieldx, "_lm")
            1:pv_ring
        elseif contains(fieldx, "_tz")
            pv_ring+1:wad_ol
        elseif contains(fieldx, "_um")
            wad_ol+1:lith_ol
        elseif contains(fieldx, "_crust")
            lith_ol+1:size(blck.rprofdata, 1)
        else
            1:size(blck.rprofdata, 1)
        end
        rprof_rangey = if contains(fieldy, "_lm")
            1:pv_ring
        elseif contains(fieldy, "_tz")
            pv_ring+1:wad_ol
        elseif contains(fieldy, "_um")
            wad_ol+1:lith_ol
        elseif contains(fieldy, "_crust")
            lith_ol+1:size(blck.rprofdata, 1)
        else
            1:size(blck.rprofdata, 1)
        end

        # Time window
        time1 = in_plates1 ? blck.platesdata[:, idxP["time"]] : in_rprof1 ? blck.rproftime : blck.timedata[:, idxT["time"]]
        time2 = in_plates2 ? blck.platesdata[:, idxP["time"]] : in_rprof2 ? blck.rproftime : blck.timedata[:, idxT["time"]]
        isnothing(tstart) && (tstart = first(time1))
        isnothing(tend) && (tend = last(time1))
        tstartidx1 = findfirst(time1 .>= tstart); (isnothing(tstartidx1)) && error("tstart=$tstart beyond data range")
        tendidx1   = findlast(time1 .<= tend); (isnothing(tendidx1)) && error("tend=$tend beyond data range")
        tstartidx2 = findfirst(time2 .>= tstart); (isnothing(tstartidx2)) && error("tstart=$tstart beyond data range")
        tendidx2   = findlast(time2 .<= tend); (isnothing(tendidx2)) && error("tend=$tend beyond data range")

        # Vectors
        xvec = in_plates1 ? blck.platesdata[:, idxP[fieldx2use]] : in_rprof1 ? vec(mean(blck.rprofdata[rprof_rangex, :, idxR[fieldx2use]], dims=1)) : blck.timedata[:, idxT[fieldx2use]]
        yvec = in_plates2 ? blck.platesdata[:, idxP[fieldy2use]] : in_rprof2 ? vec(mean(blck.rprofdata[rprof_rangey, :, idxR[fieldy2use]], dims=1)) : blck.timedata[:, idxT[fieldy2use]]

        # OM
        (fieldx == "SurfOceanMass3D") && (xvec ./= om)
        (fieldy == "SurfOceanMass3D") && (yvec ./= om)

        # Time-average
        time1 = vcat(diff(time1), 0.0)[tstartidx1:tendidx1]
        time2 = vcat(diff(time2), 0.0)[tstartidx2:tendidx2]
        tax[b] = sum(xvec[tstartidx1:tendidx1].*time1 ./ sum(time1))
        tay[b] = sum(yvec[tstartidx2:tendidx2].*time2 ./ sum(time2))

        # Normalized
        normx && (b==1) && (basex = tax[b])
        normy && (b==1) && (basey = tay[b])
    end

    # Plot
    isnothing(fig) && (fig = Figure(size = fsize))
    ax = Axis(fig[fpos[1], fpos[2]], xlabel = in_rprof1 ? LRN[fieldx] : LTN[fieldx], ylabel = in_rprof2 ? LRN[fieldy] : LTN[fieldy], xlabelsize=xlabelsize, ylabelsize=ylabelsize,
        xticklabelsize=xticklabelsize, yticklabelsize=yticklabelsize, xticksize=xticksize, yticksize=yticksize,
        xlabelpadding=xlabelpadding, ylabelpadding=ylabelpadding)
    scatter && scatter!(ax, normx ? tax./basex : tax, normy ? tay./basey : tay, color=color, markersize=markersize, marker=marker, alpha = mov_avg ? 0.3 : 1.0, strokewidth=strokewidth)
    line && lines!(ax, normx ? tax./basex : tax, normy ? tay./basey : tay, color=color, linewidth=linewidth, linestyle=linestyle, alpha = mov_avg ? 0.3 : 1.0)

    # Display and save
    disp && display(fig)
    (savein != "") && save(savein*".png", fig)
end

function ta_field_vs_field_multiple(Dblocks, fieldx::String, fieldy::String; fsize=(800, 600), xlabelsize=25, ylabelsize=25,
                            xticklabelsize=17, yticklabelsize=17, xticksize=10, yticksize=10, xlabelpadding=15, ylabelpadding=15,
                            mov_avg=false, mov_avg_window=0,subsample_size=1000, color=[:blue, :red, :green, :orange, :purple], tstart=nothing, tend=nothing,
                            scatter=true,
                                markersize=10,
                                strokewidth=1.0,
                                marker=:circle,
                            line=true,
                                linewidth=2.5,
                                linestyle=:solid,
                            fig=nothing, fpos=(1,1), disp=true, savein="",
                            setlabs=nothing,
                            )

    # Field to use
    fieldx2use = replace(fieldx, "_lm"=>"", "_tz"=>"", "_um"=>"", "_crust"=>"", "_norm"=>"")
    fieldy2use = replace(fieldy, "_lm"=>"", "_tz"=>"", "_um"=>"", "_crust"=>"", "_norm"=>"")

    # Field in plates.dat flag
    in_plates1 = ((fieldx=="Mob") || (fieldx=="Vsurf")) && !contains(fieldx, "_lm") && !contains(fieldx, "_tz") && !contains(fieldx, "_um") && !contains(fieldx, "_crust")
    in_plates2 = ((fieldy=="Mob") || (fieldy=="Vsurf")) && !contains(fieldy, "_lm") && !contains(fieldy, "_tz") && !contains(fieldy, "_um") && !contains(fieldy, "_crust")
    in_time1  = (fieldx2use in Dblocks[1][1].timeheader) && !in_plates1 && !contains(fieldx, "_lm") && !contains(fieldx, "_tz") && !contains(fieldx, "_um") && !contains(fieldx, "_crust")
    in_time2  = (fieldy2use in Dblocks[1][1].timeheader) && !in_plates2 && !contains(fieldy, "_lm") && !contains(fieldy, "_tz") && !contains(fieldy, "_um") && !contains(fieldy, "_crust")
    in_rprof1 = (fieldx2use in Dblocks[1][1].rprofheader) && !in_time1 && !in_plates1
    in_rprof2 = (fieldy2use in Dblocks[1][1].rprofheader) && !in_time2 && !in_plates2
    @assert in_plates1 || in_time1 || in_rprof1 "Field $fieldx2use not found in any header"
    @assert in_plates2 || in_time2 || in_rprof2 "Field $fieldy2use not found in any header"

    maxl = maximum(length.(Dblocks))

    tax = fill(NaN, length(Dblocks), maxl)
    tay = fill(NaN, length(Dblocks), maxl)

    # Set up normalization if wanted
    normx = contains(fieldx,"_norm")
    normy = contains(fieldy,"_norm")
    normx && (basex = zeros(Float64, length(Dblocks)))
    normy && (basey = zeros(Float64, length(Dblocks)))
    
    for set in eachindex(Dblocks)
        for (b, blck) in enumerate(Dblocks[set])

            # Get encoding
            idxT, idxR, idxP = data_encoding(blck.timeheader, blck.rprofheader, blck.platesheader)

            # Find radial phase boundaries
            pv_ring, wad_ol, lith_ol = idx_ph_transitions(blck; timeidx=1)
            rprof_rangex = if contains(fieldx, "_lm")
                1:pv_ring
            elseif contains(fieldx, "_tz")
                pv_ring+1:wad_ol
            elseif contains(fieldx, "_um")
                wad_ol+1:lith_ol
            elseif contains(fieldx, "_crust")
                lith_ol+1:size(blck.rprofdata, 1)
            else
                1:size(blck.rprofdata, 1)
            end
            rprof_rangey = if contains(fieldy, "_lm")
                1:pv_ring
            elseif contains(fieldy, "_tz")
                pv_ring+1:wad_ol
            elseif contains(fieldy, "_um")
                wad_ol+1:lith_ol
            elseif contains(fieldy, "_crust")
                lith_ol+1:size(blck.rprofdata, 1)
            else
                1:size(blck.rprofdata, 1)
            end

            # Time window
            time1 = in_plates1 ? blck.platesdata[:, idxP["time"]] : in_rprof1 ? blck.rproftime : blck.timedata[:, idxT["time"]]
            time2 = in_plates2 ? blck.platesdata[:, idxP["time"]] : in_rprof2 ? blck.rproftime : blck.timedata[:, idxT["time"]]
            isnothing(tstart) && (tstart = first(time1))
            isnothing(tend) && (tend = last(time1))
            tstartidx1 = findfirst(time1 .>= tstart); (isnothing(tstartidx1)) && error("tstart=$tstart beyond data range")
            tendidx1   = findlast(time1 .<= tend); (isnothing(tendidx1)) && error("tend=$tend beyond data range")
            tstartidx2 = findfirst(time2 .>= tstart); (isnothing(tstartidx2)) && error("tstart=$tstart beyond data range")
            tendidx2   = findlast(time2 .<= tend); (isnothing(tendidx2)) && error("tend=$tend beyond data range")

            # Vectors
            xvec = in_plates1 ? blck.platesdata[:, idxP[fieldx2use]] : in_rprof1 ? vec(mean(blck.rprofdata[rprof_rangex, :, idxR[fieldx2use]], dims=1)) : blck.timedata[:, idxT[fieldx2use]]
            yvec = in_plates2 ? blck.platesdata[:, idxP[fieldy2use]] : in_rprof2 ? vec(mean(blck.rprofdata[rprof_rangey, :, idxR[fieldy2use]], dims=1)) : blck.timedata[:, idxT[fieldy2use]]

            # Time-average
            time1 = vcat(diff(time1), 0.0)[tstartidx1:tendidx1]
            time2 = vcat(diff(time2), 0.0)[tstartidx2:tendidx2]
            tax[set, b] = sum(xvec[tstartidx1:tendidx1].*time1 ./ sum(time1))
            tay[set, b] = sum(yvec[tstartidx2:tendidx2].*time2 ./ sum(time2))

            # Normalized
            normx && (b==1) && (basex[set] = tax[set, b])
            normy && (b==1) && (basey[set] = tay[set, b])

        end
    end

    # Plot
    isnothing(fig) && (fig = Figure(size = fsize))
    ax = Axis(fig[fpos[1], fpos[2]], xlabel = in_rprof1 ? LRN[fieldx] : LTN[fieldx], ylabel = in_rprof2 ? LRN[fieldy] : LTN[fieldy], xlabelsize=xlabelsize, ylabelsize=ylabelsize,
        xticklabelsize=xticklabelsize, yticklabelsize=yticklabelsize, xticksize=xticksize, yticksize=yticksize,
        xlabelpadding=xlabelpadding, ylabelpadding=ylabelpadding, xgridvisible=false, ygridvisible=false)
    color = cpalette(:managua, length(Dblocks))
    for set in eachindex(Dblocks)
        scatter && scatter!(ax, normx ? tax[set, :]./basex[set] : tax[set, :], normy ? tay[set, :]./basey[set] : tay[set, :], color=color[set], markersize=markersize, marker=marker, label=isnothing(setlabs) ? "Set $set" : setlabs[set], strokewidth=strokewidth)
        line && lines!(ax, normx ? tax[set, :]./basex[set] : tax[set, :], normy ? tay[set, :]./basey[set] : tay[set, :], color=color[set], linewidth=linewidth, linestyle=linestyle)
    end
    axislegend(ax; position = :rb, unique=true)

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
                            fig=nothing, fpos=(1,1), disp=true, savein="", interpolate=false, np=500, xlabelsize=25, ylabelsize=25, title="")

    !isnothing(fig) && (Paxis = false)

    # Get encoding
    idxT, idxR, idxP = data_encoding(Dblock.timeheader, Dblock.rprofheader, Dblock.platesheader)
    @assert haskey(idxR, field) "Field $field not found in rprof header"

    # Vectors
    yp = logscale ? log10.(replace(Dblock.rprofdata[:,:,idxR[field]], 0.0 => 1e-8)) : Dblock.rprofdata[:,:,idxR[field]]
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
    ax = Axis(fig[Paxis ? 2 : fpos[1], fpos[2]], xlabel = L"Time\;[Gyr]", ylabel = L"Radius\;[km]", xlabelsize=xlabelsize, ylabelsize=ylabelsize, xticklabelsize=12, yticklabelsize=12, xticksize=10, yticksize=10, xlabelpadding=15, ylabelpadding=15,
                title=title, titlesize=20, titlegap=14, xticklabelrotation=pi/3.5)
    logscale && (field="log"*field)
    cmap_reverse && (cmap = Reverse(cmap))
    hm = heatmap!(ax, interpolate ? tgrid : Dblock.rproftime, interpolate ? 1e-3rgrid : 1e-3Dblock.rprofdata[:,idxR["r"],1], yp', colormap=cmap, interpolate=interpolate, colorrange=colorrange)
    Colorbar(fig[fpos[1], Paxis ? 1 : fpos[2]+1], hm; label = LRN[field], labelsize=ylabelsize, ticklabelsize=15, vertical= Paxis ? false : true, labelpadding=13)
    Paxis && second_axis!(fig, fpos, Dblock.rproftime, 1e-3abs.(Dblock.rprofdata[:,idxR["r"],1] .- Dblock.rprofdata[end,idxR["r"],1]), "Pressure", 25, 17, 15, true, false, (tstart, tend), (nothing, nothing))
    xlims!(ax, tstart, tend)
    disp && display(fig)
    (savein != "") && save(savein*".png", fig)

end

# ================================================================================
# =======                    Statistical analysis                         ========
# ================================================================================

"""
    Plot autocorrelation of mantle water content in different mantle sectors over time. Also displays e-folding time and integratal timescale.
    
    \t Basic usage: \t H2O_memory_time(Dblock::DataBlock; kwargs...)

    Optional arguments (kwargs):

      • -- Canvas arguments --

        - tstart::Float64 \t-->\t Start time for analysis [default: 0.1 Gyr]
        - tend::Union{Nothing, Float64} \t-->\t End time for analysis [default: end of simulation]
        - fsize::Tuple{Int64,Int64} \t-->\t Figure size in pixels [default: (800, 500)]
        - fpos::Tuple{Int64,Int64} \t-->\t Figure position in a grid layout [default: (1, 1)]
        - fig::Union{Nothing, Figure} \t-->\t Pre-initialized figure to plot into [default: nothing]
        - disp::Bool \t\t\t-->\t Display the figure [default: true]
        - savein::String \t\t-->\t Path to save the figure (PNG) [default: ""]
"""
function H2O_memory_time(Dblock; tstart=0.1, tend=nothing, fig=nothing, fpos=(1,1), fsize=(800,500), disp=true, savein="", legend=true, xlabvis=true, ylabvis=true)

    if typeof(Dblock) == Vector{DataBlock}
        H2O_memory_time_multiple(Dblock; tstart=tstart, tend=tend, fig=fig, fpos=fpos, fsize=fsize, disp=disp, savein=savein)
        return
    end

    # Checks and indexing setup
    @assert Dblock.metadata.H₂O_tracked "Mantle water content tracking not enabled in simulation $(Dblock.metadata.Sname)"
    idxT, idxR, idxP = data_encoding(Dblock)

    # Initialize figure 
    isnothing(fig) && (fig = Figure(size = fsize))

    # Timing indexes
    isnothing(tend) && (tend = last(Dblock.rproftime))
    timeidx_s = findfirst(Dblock.rproftime .>= tstart)
    timeidx_e = findlast(Dblock.rproftime .<= tend)

    # Gather Data
    r = Dblock.rprofdata[:, timeidx_s:timeidx_e, idxR["r"]]
    ρ = Dblock.rprofdata[:, timeidx_s:timeidx_e, idxR["rhomean"]]
    water = Dblock.rprofdata[:, timeidx_s:timeidx_e, idxR["Water"]]
    s = Dblock.rprofdata[:, timeidx_s:timeidx_e, idxR["Wsol"]]
    ∂M = Dblock.rprofdata[:, timeidx_s:timeidx_e, idxR["dM"]]
    H2Okg = water./om.*∂M.*1e-2 # wt% to fraction
    sH2Okg = s./om.*∂M.*1e-2

    # Mass per sector
    sectH2O = zeros(Float64, 4, timeidx_e - timeidx_s + 1) # Rows: LM, TZ, UM, Crust
    for t in timeidx_s:timeidx_e
        # Find radial phase boundaries
        pv_ring, wad_ol, lith_ol = idx_ph_transitions(Dblock; timeidx=t)
        sectH2O[1, t - timeidx_s + 1] = sum(H2Okg[1:pv_ring, t - timeidx_s + 1])
        sectH2O[2, t - timeidx_s + 1] = sum(H2Okg[pv_ring+1:wad_ol, t - timeidx_s + 1])
        sectH2O[3, t - timeidx_s + 1] = sum(H2Okg[wad_ol+1:lith_ol, t - timeidx_s + 1])
        sectH2O[4, t - timeidx_s + 1] = sum(H2Okg[lith_ol+1:end, t - timeidx_s + 1])
    end

    # resample to uniform grid
    tvec = LinRange(Dblock.rproftime[timeidx_s], Dblock.rproftime[timeidx_e], 500)
    vvec = zeros(Float64, 4, length(tvec))
    for sector in 1:4
        itp = interpolate((Dblock.rproftime[timeidx_s:timeidx_e],), sectH2O[sector,:], Gridded(Linear()))
        vvec[sector, :] = itp.(tvec)
    end

    # Autocorrelation
    time_lags = LinRange(1, 1000, 50)./1e3 # Myr -> Gyr
    dt_lags = round.(Int, time_lags ./ step(tvec))
    cutidx = findfirst(dt_lags .> length(tvec))
    isnothing(cutidx) ? (cutidx = length(dt_lags)) : (cutidx = max(cutidx-2, 1))
    time_lags = time_lags[1:cutidx]
    dt_lags = dt_lags[1:cutidx]
    ac_sectH2O = zeros(Float64, 4, length(time_lags))
    for sector in 1:4
        ac_sectH2O[sector, :] = StatsBase.autocor(vvec[sector,:], dt_lags, demean=true)
    end
    efold = zeros(Float64, 4)

    # ACF
    ax = Axis(fig[fpos[1], fpos[2]], xlabel = L"Time\;lag\;[\mathrm{Gyr}]", ylabel = L"Sector\;H_2O\;ACF", xlabelsize=25, ylabelsize=25, xticklabelsize=17, yticklabelsize=17, xticksize=10, yticksize=10, xlabelpadding=15, ylabelpadding=15,
            xgridvisible=false, ygridvisible=false, titlesize=20, xlabelvisible=xlabvis, ylabelvisible=ylabvis)
    cmap = :berlin
    clrs = cpalette(cmap, 4)
    for sector in 1:4
        # e-folding time
        efold_idx = findfirst(ac_sectH2O[sector, :] .<= exp(-1))
        if !isnothing(efold_idx)
            t0, t1 = time_lags[efold_idx-1], time_lags[efold_idx]
            r0, r1 = ac_sectH2O[sector, efold_idx-1], ac_sectH2O[sector, efold_idx]
            efold[sector] = (r0==r1) ? t1 : t0 + (t1 - t0) * (r0 - exp(-1)) / (r0 - r1)
        end
        lines!(ax, [efold[sector], efold[sector]], [-0.1, 1.0], color=clrs[sector], linestyle=:dash, linewidth=1.5)
        textval = @sprintf("%.2f Gyr", efold[sector])
        text!(ax, efold[sector], -0.3, text=textval, color=clrs[sector], align = (:left, :bottom), fontsize=15, rotation=pi/2)
        # Integral timescale τ
        k0 = findfirst(ac_sectH2O[sector, 2:end] .<= 0.0)
        K = (k0 === nothing) ? length(ac_sectH2O[sector,:]) : k0
        τ = 0.0
        for i in 1:(K-1)
            dt = time_lags[i+1] - time_lags[i]
            τ += 0.5 * (ac_sectH2O[sector, i] + ac_sectH2O[sector, i+1]) * dt
        end
        τs = @sprintf("%.2f Gyr", τ)
        # lines!(ax, time_lags, ac_sectH2O[sector, :], color=clrs[sector], linewidth=2.0, label=sector==1 ? "Lower Mantle (τᵢₙₜ = "*τs*")" : sector==2 ? "Mantle Transition Zone (τᵢₙₜ = "*τs*")" : sector==3 ? "Upper Mantle (τᵢₙₜ = "*τs*")" : "Crust (τᵢₙₜ = "*τs*")")
        lines!(ax, time_lags, ac_sectH2O[sector, :], color=clrs[sector], linewidth=2.0, label=sector==1 ? "Lower Mantle" : sector==2 ? "Mantle Transition Zone" : sector==3 ? "Upper Mantle" : "Crust")
        scatter!(ax, time_lags, ac_sectH2O[sector, :], markersize=10, strokewidth=1.0, color=clrs[sector], marker=:rect)
    end
    lines!(ax, [-10, -9], [0.0, 0.1], color=:black, linestyle=:dash, linewidth=1.5, label="e-folding time")
    lines!(ax, [0.0, 10000], [0.0, 0.0], color=:grey, linewidth=1.5, alpha=0.3)
    xlims!(ax, 0.0, maximum(time_lags))
    legend && axislegend(ax, position=:rt, framevisible=true, fontsize=10, padding=10, rowgap=10, orientation=:horizontal)


    disp && display(fig)
    (savein != "") && save(savein*".png", fig)
end

function H2O_memory_time_multiple(Dblock; tstart=0.1, tend=nothing, fig=nothing, fpos=(1,1), fsize=(800,500), disp=true, savein="")

    # Checks and indexing setup
    for blck in Dblock
        @assert blck.metadata.H₂O_tracked "Mantle water content tracking not enabled in simulation $(blck.metadata.Sname)"
    end
    
    τe = zeros(Float64, 4, length(Dblock))
    ef = zeros(Float64, 4, length(Dblock))
    
    for (b, blck) in enumerate(Dblock)
        idxT, idxR, idxP = data_encoding(blck)

        # Timing indexes
        isnothing(tend) && (tend = last(blck.rproftime))
        timeidx_s = findfirst(blck.rproftime .>= tstart)
        timeidx_e = findlast(blck.rproftime .<= tend)

        # Gather Data
        r = blck.rprofdata[:, timeidx_s:timeidx_e, idxR["r"]]
        ρ = blck.rprofdata[:, timeidx_s:timeidx_e, idxR["rhomean"]]
        water = blck.rprofdata[:, timeidx_s:timeidx_e, idxR["Water"]]
        s = blck.rprofdata[:, timeidx_s:timeidx_e, idxR["Wsol"]]
        ∂M = blck.rprofdata[:, timeidx_s:timeidx_e, idxR["dM"]]
        H2Okg = water./om.*∂M.*1e-2 # wt% to fraction
        sH2Okg = s./om.*∂M.*1e-2

        # Mass per sector
        sectH2O = zeros(Float64, 4, timeidx_e - timeidx_s + 1) # Rows: LM, TZ, UM, Crust
        for t in timeidx_s:timeidx_e
            # Find radial phase boundaries
            pv_ring, wad_ol, lith_ol = idx_ph_transitions(blck; timeidx=t)
            sectH2O[1, t - timeidx_s + 1] = sum(H2Okg[1:pv_ring, t - timeidx_s + 1])
            sectH2O[2, t - timeidx_s + 1] = sum(H2Okg[pv_ring+1:wad_ol, t - timeidx_s + 1])
            sectH2O[3, t - timeidx_s + 1] = sum(H2Okg[wad_ol+1:lith_ol, t - timeidx_s + 1])
            sectH2O[4, t - timeidx_s + 1] = sum(H2Okg[lith_ol+1:end, t - timeidx_s + 1])
        end

        # resample to uniform grid
        tvec = LinRange(blck.rproftime[timeidx_s], blck.rproftime[timeidx_e], 500)
        vvec = zeros(Float64, 4, length(tvec))
        for sector in 1:4
            itp = interpolate((blck.rproftime[timeidx_s:timeidx_e],), sectH2O[sector,:], Gridded(Linear()))
            vvec[sector, :] = itp.(tvec)
        end

        # Autocorrelation
        time_lags = LinRange(1, 1000, 50)./1e3 # Myr -> Gyr
        dt_lags = round.(Int, time_lags ./ step(tvec))
        cutidx = findfirst(dt_lags .> length(tvec))
        isnothing(cutidx) ? (cutidx = length(dt_lags)) : (cutidx = max(cutidx-2, 1))
        time_lags = time_lags[1:cutidx]
        dt_lags = dt_lags[1:cutidx]
        ac_sectH2O = zeros(Float64, 4, length(time_lags))
        for sector in 1:4
            ac_sectH2O[sector, :] = StatsBase.autocor(vvec[sector,:], dt_lags, demean=true)
        end
        for sector in 1:4
            # e-folding time
            efold_idx = findfirst(ac_sectH2O[sector, :] .<= exp(-1))
            if !isnothing(efold_idx)
                t0, t1 = time_lags[efold_idx-1], time_lags[efold_idx]
                r0, r1 = ac_sectH2O[sector, efold_idx-1], ac_sectH2O[sector, efold_idx]
                ef[sector, b] = (r0==r1) ? t1 : t0 + (t1 - t0) * (r0 - exp(-1)) / (r0 - r1)
            end
            # Integral timescale τ
            k0 = findfirst(ac_sectH2O[sector, 2:end] .<= 0.0)
            K = (k0 === nothing) ? length(ac_sectH2O[sector,:]) : k0
            for i in 1:(K-1)
                dt = time_lags[i+1] - time_lags[i]
                τe[sector, b] += 0.5 * (ac_sectH2O[sector, i] + ac_sectH2O[sector, i+1]) * dt
            end
        end
    end

    # tick labels
    tlabs = [Dblock[i].metadata.name for i in 1:length(Dblock)]


    # Initialize figure 
    isnothing(fig) && (fig = Figure(size = fsize))

    # ACF
    ax = Axis(fig[fpos[1], fpos[2]], ylabel = L"\tau_{e-fold}\;[\mathrm{Gyr}]", xlabelsize=25, ylabelsize=25, xticklabelsize=17, yticklabelsize=17, xticksize=10, yticksize=10, xlabelpadding=15, ylabelpadding=15,
            xgridvisible=false, ygridvisible=false, titlesize=20, xticks=(1:length(Dblock), tlabs), xticklabelrotation=pi/3)
    ax2 = Axis(fig[fpos[1], fpos[2]+1], ylabel = L"\tau_{integral}\;[\mathrm{Gyr}]", xlabelsize=25, ylabelsize=25, xticklabelsize=17, yticklabelsize=17, xticksize=10, yticksize=10, xlabelpadding=15, ylabelpadding=15,
            xgridvisible=false, ygridvisible=false, titlesize=20, xticks=(1:length(Dblock), tlabs), xticklabelrotation=pi/3)
    cmap = :Set1_4
    clrs = cpalette(cmap, 4)
    for sector in 1:4
        # integrals
        lines!(ax, 1:length(Dblock), ef[sector, :], color=clrs[sector], linewidth=2.0)
        scatter!(ax, 1:length(Dblock), ef[sector, :], markersize=20, strokewidth=1.0, color=clrs[sector], marker=:rect, label=sector==1 ? "Lower Mantle" : sector==2 ? "Mantle Transition Zone" : sector==3 ? "Upper Mantle" : "Crust")
        # xHD, yHD = regression(1:length(Dblock), ef[sector,:])
        # lines!(ax, xHD, yHD, color=clrs[sector], linewidth=1.5)
        lines!(ax2, 1:length(Dblock), τe[sector, :], color=clrs[sector], linewidth=2.0)
        scatter!(ax2, 1:length(Dblock), τe[sector, :], markersize=20, strokewidth=1.0, color=clrs[sector], marker=:rect)
        # xHD, yHD = regression(1:length(Dblock), τe[sector,:])
        # lines!(ax2, xHD, yHD, color=clrs[sector], linewidth=1.5)
    end
    Legend(fig[fpos[1]+1, fpos[2]:fpos[2]+1], ax, position=:rt, framevisible=true, fontsize=15, padding=10, rowgap=10, orientation=:horizontal)

    disp && display(fig)
    (savein != "") && save(savein*".png", fig)

end

function H2O_PCA(Dblocks, type, clrby; disp=true, fig=nothing, fpos=(1,1), fsize=(800,600), savein="", loadings_plot=false)

    if type == "twindow"

        # PC1 → Strongest variance. Loadings are ambiguous however. Most likely secular evolution
        # PC2 → Upper-mantle concentration relative to the rest
        # PC3 → How water is partitioned between TZ and deep mantle

        # Initialise 2-D matrix
        window_size = 10 * 1e-3 # Myr
        lM = 0; for blck in Dblocks
            lM += Int(ceil(last(blck.rproftime)/window_size))
            @assert blck.metadata.H₂O_tracked "Mantle water content tracking not enabled in simulation $(blck.metadata.name)"
        end
        M = zeros(Float64, lM, 4) # Columns: LM, TZ, UM, Crust
        tM, yM, oM, mM = zeros(Float64, lM), zeros(Float64, lM), zeros(Float64, lM), zeros(Float64, lM)
        tmax = 0.0

        # Extract H₂O contents
        iM = 1
        for blck in Dblocks
            # Checks and indexing setup
            idxT, idxR, idxP = data_encoding(blck)

            # Time vector
            tvec = 0.0:window_size:last(blck.rproftime)

            for t in tvec
                timeidx = findfirst(blck.rproftime .>= t)
                # Gather Data
                water = blck.rprofdata[:, timeidx, idxR["Water"]]

                # Find radial phase boundaries
                pv_ring, wad_ol, lith_ol = idx_ph_transitions(blck; timeidx=timeidx)

                # Divide by sectors
                Wlm = mean(water[1:pv_ring])
                Wtz = mean(water[pv_ring+1:wad_ol])
                Wum = mean(water[wad_ol+1:lith_ol])
                Wcrust = mean(water[lith_ol+1:end])

                # Store in matrix
                M[iM, 1] = Wlm
                M[iM, 2] = Wtz
                M[iM, 3] = Wum
                M[iM, 4] = Wcrust
                # Colors
                tM[iM] = t
                yM[iM] = blck.metadata.ystress
                oM[iM] = blck.metadata.totH₂O/om
                mM[iM] = blck.platesdata[findfirst(blck.platesdata[:,2].>=t), idxP["Mob"]]
                iM += 1

                # Track max time
                (t>tmax) && (tmax = t)
            end
        end

        # Subtract mean + divide by σ
        Mmean = mean(M, dims=1)
        Mstd = std(M, dims=1)
        Mstd[Mstd .== 0.0] .= 1.0 # Prevent division by zero
        Mnorm = (M .- Mmean) ./ Mstd

        # pca call
        pca = fit(PCA, Mnorm'; maxoutdim=3)
        scores = transform(pca, Mnorm')'
        loadings = projection(pca)

        # Figure
        GLMakie.activate!()
        nclrs = 100
        cmap = :vik100
        isnothing(fig) && (fig = Figure(size = fsize))
        if loadings_plot
            ax = Axis3(fig[fpos[1], fpos[2]], xlabel = L"Lower\;Mantle\;c^{H_2O}", ylabel = L"Mantle\;Transition\;Zone\;c^{H_2O}", zlabel = L"Upper\;Mantle\;c^{H_2O}", xlabelsize=20, ylabelsize=20, zlabelsize=20, xticklabelsize=15, yticklabelsize=15, zticklabelsize=15, xticksize=8, yticksize=8, zticksize=8, 
                        limits=((-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0)), xgridvisible=false, ygridvisible=false, zgridvisible=false, viewmode=:free)
            # PCA vectors
            lines!(ax, [0.0, loadings[1,1]], [0.0, loadings[2,1]], [0.0, loadings[3,1]], color=:red, linewidth=2.3)
            lines!(ax, [0.0, loadings[1,2]], [0.0, loadings[2,2]], [0.0, loadings[3,2]], color=:green, linewidth=2.3)
            lines!(ax, [0.0, loadings[1,3]], [0.0, loadings[2,3]], [0.0, loadings[3,3]], color=:blue, linewidth=2.3)
            scatter!(ax, [loadings[1,1]], [loadings[2,1]], [loadings[3,1]], markersize=25, color=:red, marker=:circle, label="PC-1", strokewidth=1.5)
            scatter!(ax, [loadings[1,2]], [loadings[2,2]], [loadings[3,2]], markersize=25, color=:green, marker=:circle, label="PC-2", strokewidth=1.5)
            scatter!(ax, [loadings[1,3]], [loadings[2,3]], [loadings[3,3]], markersize=25, color=:blue, marker=:circle, label="PC-3", strokewidth=1.5)
            # Dashed coordinates
            lines!(ax, [loadings[1,1], loadings[1,1]], [loadings[2,1], loadings[2,1]], [loadings[3,1], -1], color=:red, linewidth=1.7, linestyle=:dash, alpha=0.5)
            lines!(ax, [loadings[1,1], loadings[1,1]], [loadings[2,1], 1.0], [loadings[3,1], loadings[3,1]], color=:red, linewidth=1.7, linestyle=:dash, alpha=0.5)
            lines!(ax, [loadings[1,1], 1.0], [loadings[2,1], loadings[2,1]], [loadings[3,1], loadings[3,1]], color=:red, linewidth=1.27, linestyle=:dash, alpha=0.5)
            lines!(ax, [loadings[1,2], loadings[1,2]], [loadings[2,2], loadings[2,2]], [loadings[3,2], -1], color=:green, linewidth=1.7, linestyle=:dash, alpha=0.5)
            lines!(ax, [loadings[1,2], loadings[1,2]], [loadings[2,2], 1.0], [loadings[3,2], loadings[3,2]], color=:green, linewidth=1.7, linestyle=:dash, alpha=0.5)
            lines!(ax, [loadings[1,2], 1.0], [loadings[2,2], loadings[2,2]], [loadings[3,2], loadings[3,2]], color=:green, linewidth=1.7, linestyle=:dash, alpha=0.5)
            lines!(ax, [loadings[1,3], loadings[1,3]], [loadings[2,3], loadings[2,3]], [loadings[3,3], -1], color=:blue, linewidth=1.7, linestyle=:dash, alpha=0.5)
            lines!(ax, [loadings[1,3], loadings[1,3]], [loadings[2,3], 1.0], [loadings[3,3], loadings[3,3]], color=:blue, linewidth=1.7, linestyle=:dash, alpha=0.5)
            lines!(ax, [loadings[1,3], 1.0], [loadings[2,3], loadings[2,3]], [loadings[3,3], loadings[3,3]], color=:blue, linewidth=1.7, linestyle=:dash, alpha=0.5)
            # border scatters
            scatter!(ax, loadings[1,1], loadings[2,1], -1, color=:red, alpha=0.5)
            scatter!(ax, loadings[1,1], 1.0, loadings[3,1], color=:red, alpha=0.5)
            scatter!(ax, 1.0, loadings[2,1], loadings[3,1], color=:red, alpha=0.5)
            scatter!(ax, loadings[1,2], loadings[2,2], -1, color=:green, alpha=0.5)
            scatter!(ax, loadings[1,2], 1.0, loadings[3,2], color=:green, alpha=0.5)
            scatter!(ax, 1.0, loadings[2,2], loadings[3,2], color=:green, alpha=0.5)
            scatter!(ax, loadings[1,3], loadings[2,3], -1, color=:blue, alpha=0.5)
            scatter!(ax, loadings[1,3], 1.0, loadings[3,3], color=:blue, alpha=0.5)
            scatter!(ax, 1.0, loadings[2,3], loadings[3,3], color=:blue, alpha=0.5)
            # zero lines
            lines!(ax, [-1, 1], [1, 1], [0, 0], color=:grey, linewidth=1.3, alpha=0.5)
            lines!(ax, [1, 1], [-1, 1], [0, 0], color=:grey, linewidth=1.3, alpha=0.5)
            lines!(ax, [0, 0], [1, 1], [-1, 1], color=:grey, linewidth=1.3, alpha=0.5)
            lines!(ax, [0, 0], [-1, 1], [-1, -1], color=:grey, linewidth=1.3, alpha=0.5)
            lines!(ax, [-1, 1], [0, 0], [-1, -1], color=:grey, linewidth=1.3, alpha=0.5)
            lines!(ax, [1, 1], [0, 0], [-1, 1], color=:grey, linewidth=1.3, alpha=0.5)

            axislegend(ax, position=:ct, framevisible=false, fontsize=15, padding=10, rowgap=10, orientation=:horizontal)



            disp && display(fig)
            (savein != "") && save(savein*".png", fig)
        else
            (clrby=="t") && (cM = tM)
            (clrby=="ys") && (cM = yM)
            (clrby=="om") && (cM = oM)
            (clrby=="mob") && (cM = mM)
            (clrsmap = hcat(LinRange(minimum(cM), maximum(cM), nclrs), cpalette(cmap, nclrs)))
            clrs = [clrsmap[findfirst(clrsmap[:,1] .>= cM[i]), 2] for i in 1:length(cM)]
            ax = Axis3(fig[fpos[1], fpos[2]], xlabel = L"PC-1\;(Inverse\;of\;UM\;c_{H_2O}\;dominance)", ylabel = L"PC-2\;(Secular\;evolution)", zlabel = L"PC-3\;(Inverse\;of\;LM\;c_{H_2O}\;dominance)", xlabelsize=20, ylabelsize=20, zlabelsize=20, xticklabelsize=15, yticklabelsize=15, zticklabelsize=15, xticksize=8, yticksize=8, zticksize=8, viewmode=:free)
            scatter!(ax, scores[:,1], scores[:,2], scores[:,3], markersize=15, color=clrs, marker=:circle)
            # Max and mix for legend
            idxmax, idxmin = argmax(cM[:,1]), argmin(cM[:,1])
            scatter!(ax, [scores[idxmax,1]], [scores[idxmax,2]], [scores[idxmax,3]], markersize=10, color=clrs[idxmax], marker=:circle, label=clrby*" = "*@sprintf("%.2f", maximum(cM[:,1])))
            scatter!(ax, [scores[idxmin,1]], [scores[idxmin,2]], [scores[idxmin,3]], markersize=10, color=clrs[idxmin], marker=:circle, label=clrby*" = "*@sprintf("%.2f", minimum(cM[:,1])))
            axislegend(ax, position=:ct, framevisible=true, fontsize=15, padding=10, rowgap=10, orientation=:horizontal)
            disp && display(fig)
            (savein != "") && save(savein*".png", fig)
        end

        # Return
        return loadings

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
        - Xox::Array{String,1} \t-->\t List of oxides in composition [default: ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "Cr2O3", "O", "H2O"]]
        - XB::Array{Float64,1} \t-->\t StagYY's enriched endmember [default: [49.33, 15.31, 10.82, 7.41, 10.33, 0.19, 2.53, 1.46, 0.0, 0.0, 100.0]]
        - XH::Array{Float64,1} \t-->\t Bulk composition for HARZBURGITE (wt%) [default: [45.5, 2.59, 4.05, 35.22, 7.26, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0]]
        - ccomp::Array{Float64,1} \t-->\t Custom bulk composition (wt%) [default: zeros(Float64, length(Xox))]
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
function minmap(sector, em::String; nP=70, nT=70,
                Xox=["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "Cr2O3", "O", "H2O"],
                XB=[49.33, 15.31, 10.82, 7.41, 10.33, 0.19, 2.53, 1.46, 0.0, 0.0, 100.0], # From stxirtude & Bertelloni 2024
                XH=[45.5, 2.59, 4.05, 35.22, 7.26, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0], # From stxirtude & Bertelloni 2024
                ccomp=zeros(Float64, length(Xox)),DBswitchP=7.0, interp=false, cmap=:BuPu,
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
        out = multi_point_minimization(10Pv, Tv.-273.15, data, X=Xv, Xoxides=Xox, B=ones(length(Pv)), sys_in="wt", name_solvus=true, rm_list=rm_list)
    else
        data = Initialize_MAGEMin(sector=="tz" ? "sb24" : "sb24", verbose=false);
        if em=="XH"
            rm_list = remove_phases(["st"], "sb24")
            out = multi_point_minimization(10Pv, Tv.-273.15, data, X=Xv, Xoxides=Xox, name_solvus=true, rm_list=rm_list)
        else
            out = multi_point_minimization(10Pv, Tv.-273.15, data, X=Xv, Xoxides=Xox, name_solvus=true)
        end
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
    fig = Figure(size=(1500,800)); i=1; j=1
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
    Solve a single P(GPa)-T(K) point for given composition and return the minimization output.
    
    \t Basic usage: \t solve_point(P::Float64, T::Float64, em::String; kwargs...)

    Optional arguments (kwargs):

        - Xox::Array{String,1} \t-->\t List of oxides in the composition [default: ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "Cr2O3", "O", "H2O"]]
        - XB::Array{Float64,1} \t\t-->\t StagYY's default enriched endmember [default: [49.33, 15.31, 10.82, 7.41, 10.33, 0.19, 2.53, 1.46, 0.0, 0.0, 100.0]]
        - XH::Array{Float64,1} \t\t-->\t StagYY's default depleted endmember [default: [45.5, 2.59, 4.05, 35.22, 7.26, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0]]
        - ccomp::Array{Float64,1} \t-->\t Custom composition in wt% [default: zeros(Float64, length(Xox))]
        - DBswitchP::Float64 \t\t-->\t Pressure (GPa) at which to switch from upper mantle to transition zone database [default: 7.0 GPa]
        - phase_out::Array{String,1} \t-->\t List of phases to remove from the minimization [default: ["chl"]]
        - H2Osat::Bool \t\t-->\t Considers H2O saturation in upper mantle calculations [default: true]
"""
function solve_point(P, T, em;
                    Xox=["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "Cr2O3", "O", "H2O"],
                    XB=[49.33, 15.31, 10.82, 7.41, 10.33, 0.19, 2.53, 1.46, 0.0, 0.0, 100.0], # From stxirtude & Bertelloni 2024
                    XH=[45.5, 2.59, 4.05, 35.22, 7.26, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0], # From stxirtude & Bertelloni 2024
                    ccomp=zeros(Float64, length(Xox)),DBswitchP=7.0,phase_out=["chl"],H2Osat=true,R=0.02
                    )

    # Checks
    @assert length(XB) == length(Xox) "Length of Xox and XB must match."
    @assert length(XH) == length(Xox) "Length of Xox and XH must match."
    @assert em in ["XH", "XB", "Custom"] "Endmember must be either: XB, XH or Custom"
    @assert (em!="Custom" || (em=="Custom" && sum(ccomp)!=0.0)) "Please provide a custom composition if chosen"

    # Transform int inputs to real if needed
    P = Float64(P); T = Float64(T)

    # Minimize
    X = em=="XB" ? XB : em=="XH" ? XH : ccomp
    H2Osat && (X[end] = 100.0) # Set H2O to 100 for saturation calculations
    if P <= DBswitchP
        data = H2Osat ? Initialize_MAGEMin("um", verbose=false, buffer="aH2O") : Initialize_MAGEMin("um", verbose=false);
        rm_list = remove_phases(phase_out, "um")
        out = single_point_minimization(10P, T-273.15, data, X=X, Xoxides=Xox, B=1.0, name_solvus=true, rm_list=rm_list)
    else
        data = Initialize_MAGEMin(P<=25. ? "sb24" : "sb24", verbose=false);
        Xo, Xlist = oxidize_bulk(X, Xox, R; wt=false, onlyvals=false)
        if em=="XH"
            rm_list = remove_phases(["st"], "sb24")
            out = single_point_minimization(10P, T-273.15, data, X=Xo, Xoxides=Xlist, name_solvus=true, rm_list=rm_list)
        else
            out = single_point_minimization(10P, T-273.15, data, X=Xo, Xoxides=Xlist, name_solvus=true)
        end
    end
    Finalize_MAGEMin(data);

    return out
end