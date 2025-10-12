module StagPP

    # Dependencies
    using CairoMakie, DelimitedFiles, Printf, BenchmarkTools, Mmap, Parsers

    # === Structures
        include("struct.jl")
    # === Backend
        include("backend.jl")
    # === Checks
        include("checks.jl")

    export aggregate_StagData, StagData, DataBlock

end
