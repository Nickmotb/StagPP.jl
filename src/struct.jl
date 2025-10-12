# ====================
# ==== Structures ====
# ====================

struct DataBlock
    name        :: String               # Name of the data block
    timeheader  :: Array{String,1}      # Time header
    rprofheader :: Array{String,1}      # Azimuthal-averaged profile header
    timedata    :: Array{Float64,2}     # Time data (2D array)
    rprofdata   :: Array{Float64,3}     # Azimuthal-averaged profile block (nfields x 2D arrays)
end

struct StagData
    # Metadata
    name        :: String               # Name of the data set
    # Numerics
    shape       :: String               # Annulus, sperical, Yin-Yang or Cartesian
    rkm         :: Float64              # Radial length of the domain in km
    rcmb        :: Float64              # Core-mantle boundary radius in km
    nz          :: Int64                # Number of vertical levels
    tend        :: Float64              # End time of the simulation
    ndts        :: Int64                # Number of time steps
    # Simulation Parameters
    T_tracked   :: Bool                 # Whether temperature is tracked
    H₂O_tracked :: Bool                 # Whether water is tracked
    Crb_tracked :: Bool                 # Whether carbon is tracked
end