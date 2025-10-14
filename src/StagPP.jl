module StagPP

    # Dependencies
    using CairoMakie, DelimitedFiles, Printf, BenchmarkTools, Mmap, Parsers, Interpolations
    using Statistics, EasyFit

    # Constants
    const sec2Gyr = 3.1536e-17 # Gyr s⁻¹
    const om = 1.31e21 # Ocean mass in kg
    const m_s2cm_yr = 3.1536e9 # (cm yr⁻¹) / (m s⁻¹)

    # === Structures
        include("struct.jl")
    # === Backend
        include("backend.jl")
    # === Checks
        include("checks.jl")
    # === Public API
        include("api.jl")

    # Structure export
    export DataBlock
    # Constant export
    export sec2Gyr, om, m_s2cm_yr
    # Function export
    export load_sim, time_vs_field, rprof_vs_field

end
