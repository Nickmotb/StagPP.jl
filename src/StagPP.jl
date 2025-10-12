module StagPP

    # Dependencies
    using CairoMakie, DelimitedFiles, Printf, BenchmarkTools, Mmap, Parsers

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
    export aggregate_StagData, time_vs_field

end
