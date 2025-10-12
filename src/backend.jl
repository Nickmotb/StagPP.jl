# ====================
# ==== Converters ====
# ====================

# Seconds to Gyr converter
sec2Gyr(s) = s / (3.1536e16)

# =========================
# ==== Reading Backend ====
# =========================

# StagData constructor (parameters.dat)
function aggregate_StagData(Sname::String, filename::String, timedata::Array{Float64,2})
    # Read strings
    lines = split.(readlines(filename), "=")
    newlines = fill("", length(lines), 2)
    # Cleanup
    for (l, line) in enumerate(lines)
        if length(line)==1
            newlines[l, 1] = replace.(line[1], "," => "", r"\s+" => "")
        else
            [occursin.(raw"\"", line[i]) ? (line[i] = replace.(line[i], raw"\"" => "", "," => "", r"\s+" => "")) : (line[i] = replace.(line[i], "," => "", r"\s+" => "")) for i in eachindex(line)]
            newlines[l, :] = line
        end
    end
    # Pipeline
    function pass_field(f::String)
        # NOTE: Currently works for scalar Integers, Floats, Bools and Strings
        # Take value after "="
        val = newlines[findall(f.==@view newlines[:,1]), 2][1]
        # Try parse to Bool
        (length(val)==1) && (val=="T" ? (return true) : (return false))
        # Try parse to float
        field = tryparse(Int64, val); !isnothing(field) && (return field)
        field = tryparse(Float64, val); !isnothing(field) && (return field)
        # Else return string
        return val
    end
    name    = Sname
    # Numerics
    shape   = pass_field("SHAPE")
    rkm     = 1e-3pass_field("D_DIMENSIONAL")
    rcmb    = 1e-3pass_field("R_CMB")
    nz      = pass_field("NZTOT")
    tend = sec2Gyr(timedata[end,2])
    ndts = timedata[end,1]
    # Simulation Parameters
    T_tracked   = pass_field("TRACERS_TEMPERATURE")
    H₂O_tracked = pass_field("TRACEELEMENT_WATER")
    Crb_tracked = pass_field("TRACEELEMENT_CARBON")

    return StagData(name, shape, rkm, rcmb, nz, tend, ndts, T_tracked, H₂O_tracked, Crb_tracked)
end

# Backend for reading StagYY "time.dat"s
function read_StagYY_timefile(filename::String, rm_Cmass_error::Bool=true)
    # Memory map setup
        map = Mmap.mmap(filename)
    # --- Compute header + fixed line bytes
        lf = findnext(==(0x0A), map, 1)
        lf2 = findnext(==(0x0A), map, lf+1)
        @assert !isnothing(lf) "$filename appears to be empty"
        @assert !isnothing(lf2) "$filename appears to have no data"
        line_bytes = lf2 - lf
    # --- Extract Header content
        header = split(String(@view map[1:(map[lf-1]==0x0D ? lf-2 : lf-1)]))
        Cmass_idx = findfirst(header.=="Cmass_error")
    # --- Compute number of rows and Character range for columns
        nrows = (length(map) - lf) ÷ line_bytes

    # Preallocate data array
        data = Array{Float64}(undef, nrows, length(header))

    # Chunk configuration
        runners = Threads.nthreads()
        div, rem = divrem(nrows, runners)
        startrow(r) = (r-1)*div + min(r-1, rem) + 1
        nrows_on_thread(r) = div + (r <= rem ? 1 : 0)
        endrow(r) = startrow(r) + nrows_on_thread(r) - 1

    # Parallel-read
        Threads.@threads for r in 1:runners
            # Thread range
            sr, er = startrow(r), endrow(r)
            # Destination block
            dest = @view data[sr:er, :]
            # Chunk stream
            @inbounds for i in sr:er
                # Define row in bytes
                srow = lf + (i-1)*line_bytes + 1
                erow = srow + line_bytes - 1
                # Remove endline characters (\n = 0x0A) for both Unix and Windows + extra windows return if present (\r = 0x0D).
                e = erow
                (map[e]==0x0A) && (e -= 1)
                (e>=srow && map[e]==0x0D) && (e -= 1)
                # Extract row string
                rowstr = String(@view map[srow:e])
                colranges = collect(findall(r"\S+", rowstr))
                # Column parsing
                j = 1
                @inbounds for rcol in colranges
                    if rm_Cmass_error && !isnothing(Cmass_idx) && j==Cmass_idx
                        dest[i-sr+1, j] = 0.0
                    else
                        fld = @view rowstr[rcol]
                        dest[i-sr+1, j] = Parsers.parse(Float64, fld)
                    end
                    j+=1
                end
            end
        end
    
    return data, header
end

# Backend for reading StagYY "rprof.dat"s
function read_StagYY_rproffile(filename::String, Stag::StagData)
    # Memory map setup
        map = Mmap.mmap(filename)
    # --- Compute header + slice separator + fixed line bytes
        lf = findnext(==(0x0A), map, 1)
        lf2 = findnext(==(0x0A), map, lf+1)
        lf3 = findnext(==(0x0A), map, lf2+1)
        @assert !isnothing(lf) "$filename appears to be empty"
        @assert !isnothing(lf2) "$filename appears to have no time separator"
        @assert !isnothing(lf3) "$filename appears to have no profile data"
        line_bytes, separator_bytes = lf3 - lf2, lf2 - lf
    # --- Extract Header content
        header = split(String(@view map[1:(map[lf-1]==0x0D ? lf-2 : lf-1)]))
    # --- Compute number of rows and byte range for separator lines (Export rate of rprof is different than time.dat)
        nblocks = (length(map) - lf) ÷ (Stag.nz*line_bytes + separator_bytes) # Number of data blocks

    # Preallocate data array
        data = Array{Float64}(undef, Stag.nz, nblocks, length(header))
        time = Array{Float64}(undef, nblocks)

    # Chunk configuration
        runners = Threads.nthreads()
        div, rem = divrem(nblocks, runners)
    # --- Raw file coordinates
        startblock(r) = (r-1)*div + min(r-1, rem) + 1
        nblocks_on_thread(r) = div + (r <= rem ? 1 : 0)
        endblock(r) = startblock(r) + nblocks_on_thread(r) - 1
    # --- Block coordinates
        startbyte(b) = lf + (b-1)*(Stag.nz*line_bytes+separator_bytes) + 1
        function endbyte(b) 
            e = (startbyte(b) + (Stag.nz*line_bytes+separator_bytes) - 1)
            (map[e]==0x0A) && (e -= 1)
            (e>=startbyte(b) && map[e]==0x0D) && (e -= 1)
            return e
        end
        block_to_byte_range(b) = startbyte(b) : endbyte(b)
        function row_to_byte_range(row, block_start_byte)
            Δ = separator_bytes + (row-1)*line_bytes
            return (block_start_byte + Δ) : (block_start_byte + Δ + line_bytes - 2)
        end

    # Parallel-read
    Threads.@threads for r in 1:runners
        # Thread range
        sb, eb = startblock(r), endblock(r)
        # Destination view array
        dest = @view data[:, sb:eb, :]
        # Block iterator
        @inbounds for block in sb:eb
            # Block range
            block_range = block_to_byte_range(block)
            # Time entry
            time_range = first(block_range):first(block_range)+separator_bytes
            t = tryparse(Float64, split(String(map[time_range]))[6])
            @assert !isnothing(t) "Failed to parse time entry in rprof.dat"
            time[block] = t
            # Row iterator
            @inbounds for row in 1:Stag.nz
                row_range = row_to_byte_range(row, first(block_range))
                fld = split(String(@view map[row_range]))
                dest[row, block-sb+1, :] .= Parsers.parse.(Float64, fld)
            end
        end
    end
    return data, header, time

end

# High-level output reader
function aggregate_StagData(sroot::String, Sname::String)
    Tfname = sroot * "_time.dat"
    Rfname = sroot * "_rprof.dat"
    time_data, time_header = read_StagYY_timefile(Tfname);
    Stag = aggregate_StagData(Sname, sroot * "_parameters.dat", time_data);
    rprof_data, rprof_header, rprof_time = read_StagYY_rproffile(Rfname, Stag);
    time_data[:,2] .= sec2Gyr.(time_data[:,2]); # Convert time to Gyr
    rprof_time .= sec2Gyr.(rprof_time); # Convert time to Gyr
    return Stag, DataBlock(Sname, time_header, rprof_header, time_data, rprof_data, rprof_time)
end

# ==================
# ==== Encoding ====
# ==================

# Create dictionary encoding for time and rprof blocks
function data_encoding(time_header, rprof_header)
    # Time header encoding
        time_encoding = Dict{String,Int64}()
        [time_encoding[name] = i for (i, name) in enumerate(time_header)]
    # Rprof header encoding
        rprof_encoding = Dict{String,Int64}()
        [rprof_encoding[name] = i for (i, name) in enumerate(rprof_header)]
    return time_encoding, rprof_encoding
end

# ======================
# ==== Auxilliaries ====
# ======================
