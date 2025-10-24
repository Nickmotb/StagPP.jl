# ========================
# ==== API Structures ====
# ========================

# Metadata
struct StagData
    # Metadata
    name        :: String                   # Name of the data set
    # Numerics
    shape       :: String                   # Annulus, sperical, Yin-Yang or Cartesian
    rkm         :: Float64                  # Radial length of the domain in km
    rcmb        :: Float64                  # Core-mantle boundary radius in km
    nz          :: Int64                    # Number of vertical levels
    tend        :: Union{Float64, Nothing}  # End time of the simulation
    ndts        :: Union{Int64, Nothing}    # Number of time steps
    totH₂O        :: Float64                # Total water content in kg
    # Simulation Parameters
    T_tracked   :: Bool                     # Whether temperature is tracked
    H₂O_tracked :: Bool                     # Whether water is tracked
    Crb_tracked :: Bool                     # Whether carbon is tracked
end

# Simulation data block
struct DataBlock
    name        :: String               # Name of the data block
    # Headers
    timeheader  :: Union{Array{String,1}, Nothing}      # Time header
    rprofheader :: Union{Array{String,1}, Nothing}      # Azimuthal-averaged profile header
    platesheader:: Union{Array{String,1}, Nothing}      # Plates header (can be nothing)
    # Data
    timedata    :: Union{Array{Float64,2}, Nothing}     # Time data (2D array)
    rprofdata   :: Union{Array{Float64,3}, Nothing}     # Azimuthal-averaged profile block (nfields x 2D arrays)
    platesdata  :: Union{Array{Float64,2}, Nothing} # Plates block (nfields x 2D arrays) (can be nothing)
    # Extra
    rproftime   :: Union{Array{Float64,1}, Nothing}     # Time vector for rprof/plates data (Sampled)
    rprofdV     :: Union{Array{Float64,1}, Nothing}     # Radial volume differentials for rprof data
    # Metadata
    metadata    :: StagData     # Metadata structure
end

# =============================
# ==== sᴴ²ᴼ/fO₂ Structures ====
# =============================

struct sᴴ²ᴼ
    # maps
    um       :: Array{Float64,3}    # Upper mantle sᴴ²ᴼ grid (nP x nT)
    tz       :: Array{Float64,3}    # Transition zone sᴴ²ᴼ grid (nP x nT)
    lm       :: Array{Float64,3}    # Lower mantle sᴴ²ᴼ grid (nP x nT)
    # Vectors
    Pum      :: Array{Float64,1}    # Upper mantle pressure vector (kbar)
    Tum      :: Array{Float64,1}    # Upper mantle temperature vector (K
    Ptz      :: Array{Float64,1}    # Transition zone pressure vector (kbar)
    Ttz      :: Array{Float64,1}    # Transition zone temperature vector (
    Plm      :: Array{Float64,1}    # Lower mantle pressure vector (kbar)
    Tlm      :: Array{Float64,1}    # Lower mantle temperature vector (    
end

struct min_sᴴ²ᴼ{T}
    # Interpolator fields
    ol          :: T
    wad         :: T
    rw          :: T
    grt         :: T
    opx         :: T
    opx_al      :: T
    cpx_lp_di   :: T
    cpx_lp_jd   :: T
    cpx_hp      :: T
    st          :: T
    CaCl₂_st    :: T
    α_PbO₂_st   :: T
    coe         :: T
    Dppv_al     :: T
    Dppv_noal   :: T
    # Constant fields
    pv          :: Float64
    cpv         :: Float64
    cf          :: Float64
    crst        :: Float64
    fp          :: Float64
    D_rw_aki    :: Float64
    cor         :: Float64
end