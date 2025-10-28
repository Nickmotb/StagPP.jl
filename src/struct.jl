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
    # DHMS
    PhA         :: Float64
    PhE         :: Float64
    shB         :: Float64
    PhD         :: Float64
    PhH         :: Float64
end

struct DHMS_boundaries{T}
    # See Komabayashi 2005
    # --- Entries
        r1 :: T     # 5 Atg                     →   14 PhA + 142 En + 113 H₂O   (2 on paper)
        r2 :: T     # 5790 phA + 21831 En       →   18757 hy-Wad + 11630 phE    (8 on paper)
        r3 :: T     # 2 phA + 11 Fo             →   3 shB + 6 En                (7 on paper)
        r4 :: T     # 19265 phE + 62762 hy-Wad  →   16288 shB + 38172 St        (55 on paper)
        r5 :: T     # 1017 shB + 7839 St        →   1800 phD + 8046 Pv          (83 on paper)
        r6 :: T     # 1 PhD                     →   1 PhH + 1 St                (decomposition)
    # --- Exit
        e :: T
end

struct PTpath{T}
    T           :: T                   # Temperature (K) at requrest P (GPa). Interpolator object
    Bᵢ          :: Float64             # Initial antigorite budget of path
    rh          :: Array{String, 1}    # Reaction history of path, entries as string numbers
    en_seed     :: Float64             # Max seen enstatite to seed for chaining outside stable range
    wad_seed    :: Float64             # Max seen wadsleyite to seed for chaining outside stable range
    fo_seed     :: Float64             # Max seen forsterite to seed for chaining outside stable range
    st_seed     :: Float64             # Max seen stishovite to seed for chaining outside stable range
end