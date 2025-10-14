# time.dat field --> String
LTN = Dict(
    "time"  => L"Evolution\;time\;[Gyr]",
    "F_top" => L"Surface\;heat\;flux\;[mW/m^2]",
    "F_bot" => L"CMB\;heat\;flux\;[mW/m^2]",
    "Tmean" => L"Mean\;mantle\;temperature\;[K]",
    "Vrms"  => L"RMS\;velocity\;[cm/yr]",
    "eta_amean" => L"Mean\;viscosity\;[Pa\cdot s]",
    "eta_mean"  => L"Mean\;viscosity\;[Pa\cdot s]",
    "erupta" => L"Erupted\;material\;[kg]",
    "erupt_rate" => L"Eruption\;rate\;[kg/s]",
    "Tcmb" => L"CMB\;temperature\;[K]",
    "Tsurf" => L"Surface\;temperature\;[K]",
    "Tpotl" => L"Potential\;temperature\;[K]",
    "cH2O_mean" => L"Mean\;mantle\;H_2O\;content\;[wt%]",
    "IngassedH2O" => L"Total\;ingassed\;H_2O\;[kg]",
    "SaturationOutgassH2O" => L"Saturation\;outgassed\;H_2O\;[kg]",
    "EruptedH2O" => L"Erupted\;H_2O\;[kg]",
    "SurfOceanMass3D" => L"Surface\;ocean\;mass\;[kg]",
    "TotH2OMassOceanPlusMantle" => L"Total\;system\;H_2O\;[kg]",
    "Mob" => L"Lithosphere\;mobility\;[v_{surf}/v_{rms}]",
    "Vsurf" => L"Surface\;velocity\;[cm/yr]",
)

# rprof.dat field --> String
LRN = Dict(
    "Tmean" => L"Mean\;mantle\;temperature\;[K]", "logTmean" => L"Mean\;mantle\;temperature\;[log_{10}K]",
    "Vrms"  => L"RMS\;velocity\;[cm/yr]", "logvrms" => L"RMS\;velocity\;[log_{10}(cm/yr)]",
    "eta_amean" => L"Mean\;viscosity\;[Pa\cdot s]", "logeta_amean" => L"Mean\;viscosity\;[log_{10}(Pa\cdot s)]",
    "eta_mean"  => L"Mean\;viscosity\;[Pa\cdot s]", "logeta_mean"  => L"Mean\;viscosity\;[log_{10}(Pa\cdot s)]",
    "Water" => L"H_2O\;content\;[wt%]", "logWater" => L"H_2O\;content\;[log_{10}wt%]",
    "Wsol" => L"H_2O\;storage\;capacity\;[wt%]", "logWsol" => L"H_2O\;storage\;capacity\;[log_{10}wt%]",
    "fO2" => L"fO_2\;[log_{10}FMQ]", "logfO2" => L"fO_2\;[log_{10}FMQ]",
    "H2Ofree" => L"Free\;H_2O\;[wt%]", "logH2Ofree" => L"Free\;H_2O\;[log_{10}wt%]",
    "rhomean" => L"Mean\;density\;[kg/m^3]", "logrhomean" => L"Mean\;density\;[log_{10}(kg/m^3)]",
    "bsmean" => L"Mean\;basalt\;fraction", "logbsmean" => L"Mean\;basalt\;fraction\;[log_{10}]",
)

# =============
# ==== API ====
# =============

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
        - tstart::Union{Float64, Nothing} \t-->\t Start time [default: nothing (first time)]
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
    fig = Figure(size = fsize)
    ax = Axis(fig[1, 1], xlabel = L"Time\;[Gyr]", ylabel = LTN[field], xlabelsize=xlabelsize, ylabelsize=ylabelsize, xgridvisible=xgrid, ygridvisible=ygrid, yreversed=yreversed,
               xticklabelsize=xticklabelsize, yticklabelsize=yticklabelsize, xticksize=xticksize, yticksize=yticksize, xlabelpadding=xlabelpadding, ylabelpadding=ylabelpadding)
    scatter && scatter!(ax, xvec, yvec, color=color, markersize=markersize, marker=marker, alpha = mov_avg ? 0.3 : 1.0)
    line && lines!(ax, xvec, yvec, color=color, linewidth=linewidth, linestyle=linestyle, alpha = mov_avg ? 0.3 : 1.0)
    mov_avg && lines!(ax, xvec, EasyFit.movavg(yvec, mov_avg_window).x, color=color, linewidth=linewidth, linestyle=:solid)
    tstartidx = findfirst(xvec .>= tstart); (isnothing(tstartidx)) && error("tstart=$tstart beyond data range")
    tendidx   = findlast(xvec .<= tend); (isnothing(tendidx)) && error("tend=$tend beyond data range")
    localmin, localmax = 0.9minimum(yvec[tstartidx:tendidx]), 1.1maximum(yvec[tstartidx:tendidx])
    second_axis!(fig, xvec, yvec, field, ylabelsize, yticklabelsize, ylabelpadding, yreversed, true, (tstart, tend), (localmin, localmax))
    xlims!(ax, tstart, tend)
    (tstart!=0. || tend!=4.5) && (yreversed ? ylims!(ax, localmax, localmin) : ylims!(ax, localmin, localmax))
    display(fig)
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
        - log::Bool \t\t\t-->\t Plots log₁₀ of the field [default: false]
        - tstart::Union{Float64, Nothing} \t-->\t Start time [default: nothing (first time)]
        - tend::Union{Float64, Nothing} \t-->\t End time [default: nothing (last time)]

"""
function rprof_vs_field(Dblock::DataBlock, field::String; fsize=(900, 600), cmap=:vik100, log=false, cmap_reverse=false, Paxis=true, tstart=nothing, tend=nothing)

    # Get encoding
    idxT, idxR, idxP = data_encoding(Dblock.timeheader, Dblock.rprofheader, Dblock.platesheader)
    @assert haskey(idxR, field) "Field $field not found in rprof header"

    # Vectors
    yp = log ? log10.(Dblock.rprofdata[:,:,idxR[field]]) : Dblock.rprofdata[:,:,idxR[field]]
    isnothing(tstart) && (tstart = first(Dblock.rproftime))
    isnothing(tend) && (tend = last(Dblock.rproftime))

    # Plot
    fig = Figure(size = fsize)
    ax = Axis(fig[Paxis ? 2 : 1, 1], xlabel = L"Time\;[Gyr]", ylabel = L"Radius\;[km]", xlabelsize=25, ylabelsize=25, xticklabelsize=17, yticklabelsize=17, xticksize=10, yticksize=10, xlabelpadding=15, ylabelpadding=15)
    log && (field="log"*field)
    cmap_reverse && (cmap = Reverse(cmap))
    hm = heatmap!(ax, Dblock.rproftime, 1e-3Dblock.rprofdata[:,idxR["r"],1], yp', colormap=cmap)
    Colorbar(fig[1, Paxis ? 1 : 2], hm; label = LRN[field], labelsize=25, ticklabelsize=15, vertical= Paxis ? false : true, labelpadding=13)
    Paxis && second_axis!(fig, Dblock.rproftime, 1e-3abs.(Dblock.rprofdata[:,idxR["r"],1] .- Dblock.rprofdata[end,idxR["r"],1]), "Pressure", 25, 17, 15, true, false, (tstart, tend), (nothing, nothing))
    xlims!(ax, tstart, tend)
    display(fig)

end

# ======================
# ==== Auxilliaries ====
# ======================

function second_axis!(fig, x, y, field, ylabelsize, yticklabelsize, ylabelpadding, yreversed, timeplot, trange, locals)

    function transform_data(y, field)
        yp = [first(y), last(y)]
        (field=="SurfOceanMass3D") && (return yp./om, L"Ocean\;masses") # 10^21 kg
        (field=="Pressure") && (return km2GPa.(yp), L"Pressure\;[GPa]") # GPa
        return nothing, nothing
    end
    yp, ylab = transform_data(y, field); (isnothing(yp)) && return
    ax2 = Axis(fig[timeplot ? 1 : 2, 1], yticklabelcolor = :black, yaxisposition = :right, ylabel = ylab, ylabelsize=ylabelsize, yticklabelsize=yticklabelsize, ylabelpadding=ylabelpadding,
                xgridvisible=false, ygridvisible=false, yreversed=yreversed)
    scatter!(ax2, first(x)*ones(2), yp, alpha=0.0)
    yreversed && ylims!(ax2, maximum(yp), minimum(yp))
    (trange[1]!=0. || trange[2]!=4.5) && timeplot && (yreverse ? ylims!(ax2, reverse(locals)) : (ylims!(ax2, locals)))
    xlims!(ax2, trange)
    hidespines!(ax2); hidexdecorations!(ax2)
end

# Convert depth in km to pressure in GPa using PREM
function km2GPa(depth)

    r = 6371 - max(depth,0.0)
    r < 0.0 ? error("Depth beyond Earth's radius") : nothing
    #  Inner core
    if (r >= 0.0 && r <= 1221.5)
    pressure = -575.170316 + 0.29421732*depth -2.30448727e-05*depth*depth
    
    #  Outer core
    elseif (r >= 1221.5 && r < 3480.0)
    pressure = -293.862745 +0.183037581*depth  -1.202957e-05*depth*depth
    
    #  Lower mantle
    elseif (r >= 3480.0 && r <= 3630.0)
    pressure = -33.49617 + 0.05853999*depth
    elseif (r >= 3630.0 && r <= 5600.0)
    pressure = -3.56294332+0.0390269653*depth+3.11829162e-06*depth*depth
    elseif (r >= 5600.0 && r <= 5701.0)
    pressure =  -5.758207 + 0.04415793*depth
    
    # Transition zone
    elseif (r >= 5701.0 && r <= 5771.0)
    pressure =  -2.875294 + 0.03985716*depth
    elseif (r >= 5771.0 && r <= 5971.0)
    pressure =  -2.056346 + 0.03845555*depth
    elseif (r >= 5971.0 && r <= 6151.0)
    pressure = -0.5319232 + 0.03467901*depth
    
    # LVZ & LID
    elseif (r >= 6151.0 && r <= 6346.6)
    pressure =   -0.213082 + 0.0332782*depth
    elseif (r >= 6346.6 && r <= 6371.0)
    pressure = 2.6e3 * 9.8 * (6371.0 - r ) *1000/1.0e9 # rho=2.6g/cm3, in GPa
    end
    
    return pressure
end

# Convert pressure in GPa to depth in km using PREM
function GPa2km(pressure)
    depth = LinRange(0., 6371., 800)
    P = sort(km2GPa.(depth))
    itp = interpolate((P,), depth, Gridded(Linear()))
    return itp(pressure)
end