LTN = Dict(
    "time"  => L"Evolution\;time\;[Gyr]",
    "F_top" => L"Surface\;heat\;flux\;[mW/m^2]",
    "F_bot" => L"CMB\;heat\;flux\;[mW/m^2]",
    "Tmean" => L"Mean\;mantle\;temperature\;[K]",
    "Vrms"  => L"RMS\;velocity\;[m/s]",
    "eta_amean" => L"Mean\;viscosity\;[Pa\cdot s]",
    "eta_mean"  => L"Mean\;viscosity\;[Pa\cdot s]",
    "erupta" => L"Erupted\;material\;[kg]",
    "erupt_rate" => L"Eruption\;rate\;[kg/s]",
    "T_cmb" => L"CMB\;temperature\;[K]",
    "T_surf" => L"Surface\;temperature\;[K]",
    "Tpotl" => L"Potential\;temperature\;[K]",
    "cH2O_mean" => L"Mean\;mantle\;H_2O\;content\;[wt%]",
    "IngassedH2O" => L"Total\;ingassed\;H_2O\;[kg]",
    "SaturationOutgassedH2O" => L"Saturation\;outgassed\;H_2O\;[kg]",
    "EruptedH2O" => L"Erupted\;H_2O\;[kg]",
    "SurfOceanMass3D" => L"Surface\;ocean\;mass\;[kg]",
    "TotH2OMassOceanPlusMantle" => L"Total\;system\;H_2O\;[kg]",
)

# =============
# ==== API ====
# =============

function time_vs_field(Dblock::DataBlock, field::String; subsample=0, color=:blue,
                            scatter=false,
                                markersize=10, marker=:circle,
                            line=true,
                                linewidth=2, linestyle=:solid
                            )
    # Get encoding
    idxT, idxR = data_encoding(Dblock.timeheader, Dblock.rprofheader)
    @assert haskey(idxT, field) "Field $field not found in time header"

    # Automatic subsample
    (subsample==0) && (subsample = max(1, Int(size(Dblock.timedata, 1) ÷ 500)))

    # Plot
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], xlabel = L"Time\;[Gyr]", ylabel = LTN[field])
    scatter && scatter!(ax, Dblock.timedata[1:subsample:end,2], Dblock.timedata[1:subsample:end, idxT[field]], color = color, markersize=markersize, marker=marker)
    line && lines!(ax, Dblock.timedata[1:subsample:end,2], Dblock.timedata[1:subsample:end, idxT[field]], color = color, linewidth=linewidth, linestyle=linestyle)
    display(fig)
end