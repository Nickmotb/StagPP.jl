module StagPP

    # Dependencies
    using Mmap, Parsers, Interpolations, ForwardDiff, LinearAlgebra, StaticArrays, LightXML
    using Statistics, EasyFit, LsqFit, StatsBase, MultivariateStats, MAGEMin_C, ColorSchemes
    using CairoMakie, GLMakie, DelimitedFiles, Printf, NumericalIntegration

    # Constants
    const R = 8.31446261815324 # J mol⁻¹ K⁻¹
    const H₂O_mm = 18.01528 # g/mol
    const sec2Gyr = 3.1536e-17 # Gyr s⁻¹
    const om = 1.31e21 # Ocean mass in kg
    const m_s2cm_yr = 3.1536e9 # (cm yr⁻¹) / (m s⁻¹)
    const ocD2toD3 = 244.45 # Conversion factor of surface budgets from 2D to 3D
    const rootdir = @__DIR__
    const savedir = joinpath(rootdir, "../", "+op"); !isdir(savedir) && mkdir(savedir)
    # Molar masses (g/mol)
    const mm = Dict(
            "SiO2" => 60.08, "Al2O3" => 101.96, "CaO" => 56.08, "MgO" => 40.30, "FeO" => 71.85, "Fe2O3" => 159.69,
            "K2O" => 94.2, "Na2O" => 61.98, "TiO2" => 79.88, "O" => 16.0, "Cr2O3" => 151.99, "MnO" => 70.937,
            "H2O" => 18.015, "CO2" => 44.01, "S" => 32.06, "P2O5" => 141.9445, "Fe" => 55.845
            )

    # === Structures
        include("struct.jl")
    # === Backend
        include("backend.jl")
    # === Checks
        include("checks.jl")
    # === sᴴ²ᴼ
        include("sf_H2O.jl")
    # === Misc
        include("misc.jl")
    # === Public API
        include("api.jl")

    # Structure export
    export DataBlock, sᴴ²ᴼ, min_sᴴ²ᴼ
    # Constant export
    export sec2Gyr, om, m_s2cm_yr
    # Function export
    export load_sim, load_local, data_encoding, solve_sH2O_fO2, min_sᴴ²ᴼ_assembler, solve_point
    # Misc export
    export SM_fO2_Oex_solution_space, Hirschmann_fO2_to_R
    # Plot exports
    export time_vs_field, rprof_vs_field, field_vs_field, ta_field_vs_field, ta_series_vs_field, 
            mantle_water, mantle_water_at_t, plot_sf, minmap, 
            snapshot, IOplot, omplot, H2O_memory_time, H2O_memory_time_multiple, H2O_PCA,
    # Auxilliaries
            idx_ph_transitions, t_avg, oxidize_bulk, c2k, k2c, km2GPa, GPa2km
end
