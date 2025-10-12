module StagPP

    # Dependencies
    using CairoMakie, DelimitedFiles, Printf, BenchmarkTools, Mmap, Parsers

    # === Structures
        include("struct.jl")
    # === Backend
        include("backend.jl")
    # === Checks
        include("checks.jl")

    export read_StagYY_timefile, aggregate_StagData, StagData

end
