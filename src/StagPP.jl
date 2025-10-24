module StagPP

    # Dependencies
    using Mmap, Parsers, Interpolations, ForwardDiff
    using Statistics, EasyFit, LsqFit, MAGEMin_C, ImageFiltering
    using CairoMakie, DelimitedFiles, Printf, BenchmarkTools

    # Constants
    const R = 8.31446261815324 # J mol⁻¹ K⁻¹
    const H₂O_mm = 18.01528 # g/mol
    const sec2Gyr = 3.1536e-17 # Gyr s⁻¹
    const om = 1.31e21 # Ocean mass in kg
    const m_s2cm_yr = 3.1536e9 # (cm yr⁻¹) / (m s⁻¹)

    # === Structures
        include("struct.jl")
    # === Backend
        include("backend.jl")
    # === Checks
        include("checks.jl")
    # === sᴴ²ᴼ
        include("s_H2O.jl")
    # === fO₂
        include("fO2.jl")
    # === Public API
        include("api.jl")

    # Structure export
    export DataBlock
    # Constant export
    export sec2Gyr, om, m_s2cm_yr
    # Function export
    export load_sim, data_encoding
    # Plot exports
    export time_vs_field, rprof_vs_field, field_vs_field, mantle_water
end
