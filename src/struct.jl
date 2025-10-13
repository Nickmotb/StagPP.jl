# ====================
# ==== Structures ====
# ====================

struct DataBlock
    name        :: String               # Name of the data block
    # Headers
    timeheader  :: Union{Array{String,1}, Nothing}      # Time header
    rprofheader :: Union{Array{String,1}, Nothing}      # Azimuthal-averaged profile header
    platesheader:: Union{Array{String,1}, Nothing}      # Plates header (can be nothing)
    # Data
    timedata    :: Union{Array{Float64,2}, Nothing}     # Time data (2D array)
    rprofdata   :: Union{Array{Float64,3}, Nothing}     # Azimuthal-averaged profile block (nfields x 2D arrays)
    platesdata  :: Union{Array{Float64,3}, Nothing} # Plates block (nfields x 2D arrays) (can be nothing)
    # Extra
    rproftime   :: Union{Array{Float64,1}, Nothing}     # Time vector for rprof/plates data (Sampled)
end

struct StagData
    # Metadata
    name        :: String                   # Name of the data set
    # Numerics
    shape       :: String                   # Annulus, sperical, Yin-Yang or Cartesian
    rkm         :: Float64                  # Radial length of the domain in km
    rcmb        :: Float64                  # Core-mantle boundary radius in km
    nz          :: Int64                    # Number of vertical levels
    tend        :: Union{Float64, Nothing}  # End time of the simulation
    ndts        :: Int64                    # Number of time steps
    # Simulation Parameters
    T_tracked   :: Bool                     # Whether temperature is tracked
    H₂O_tracked :: Bool                     # Whether water is tracked
    Crb_tracked :: Bool                     # Whether carbon is tracked
end