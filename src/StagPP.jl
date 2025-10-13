module StagPP

    # Dependencies
    using CairoMakie, DelimitedFiles, Printf, BenchmarkTools, Mmap, Parsers, Interpolations

    # Constants
    const sec2Gyr = 3.1536e-17 # Gyr s⁻¹
    const om = 1.37e21 # Ocean mass in kg
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
    export StagData, DataBlock
    # Function export
    export aggregate_StagData, time_vs_field, rprof_vs_field

end
