module StagPP

    # Dependencies
    using Mmap, Parsers, Interpolations, ForwardDiff, LinearAlgebra, StaticArrays, LightXML
    using Statistics, EasyFit, LsqFit, StatsBase, MultivariateStats, MAGEMin_C, ColorSchemes
    using GLMakie, DelimitedFiles, Printf, NumericalIntegration, LaTeXStrings

    # === Structures
    include("struct.jl")

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
    const mm = Cbulk(60.08, 101.96, 71.85, 55.85, 56.08, 101.96, 61.98, 151.99, 159.69, 16.0, 79.88, 94.2, 18.015, 44.009, 12.011)
    const Cbulk_iterator = String.(fieldnames(Cbulk))
    

    # === Backend
        include("backend.jl")
    # === Checks
        include("checks.jl")
    # === Public API
        include("api.jl")
    # === sᴴ²ᴼ
        include("sf_H2O.jl")
    # === Misc
        include("misc.jl")

    # Structure export
    export DataBlock, sᴴ²ᴼ, min_sᴴ²ᴼ
    # Constant export
    export sec2Gyr, om, m_s2cm_yr
    # Function export
    export load_sim, load_local, data_encoding, solve_sH2O_fO2, min_sᴴ²ᴼ_assembler, solve_point
    # Misc export
    export partition_Oₑₓ, P_T_ϕ_TOₑₓ_TC_topoplogy
    # Plot exports
    export time_vs_field, rprof_vs_field, field_vs_field, ta_field_vs_field, ta_series_vs_field, 
            mantle_water, mantle_water_at_t, plot_sf, minmap, time_vs_ratio,
            snapshot, IOplot, H2O_sector_ACF, H2O_PCA,
    # Auxilliaries
            idx_ph_transitions, t_avg, oxidize_bulk, c2k, k2c, km2GPa, GPa2km
end
