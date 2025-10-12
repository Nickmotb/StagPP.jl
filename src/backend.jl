# Seconds to Gyr
sec2Gyr(s) = s / (3.1536e16)

# StagData constructor (hardcoded)
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
        nrows = nblocks * (Stag.nz+1) + 1 # Total number of rows including separators
        separator_ranges = Vector{UnitRange{Int}}(undef, nblocks) # Byte ranges for separators
        sepline_idx_at_tstep(t) = (lf+1) + (t-1)*(separator_bytes + Stag.nz*line_bytes)
        for block in 1:nblocks
            s, e = sepline_idx_at_tstep(block), sepline_idx_at_tstep(block)+separator_bytes-2
            test = String(@view map[s:(map[e]==0x0D ? e-1 : e)])
            separator_ranges[block] = s:e
        end

    # Preallocate data array
        data = Array{Float64}(undef, Stag.nz, nblocks, length(header))

    # Chunk configuration
        runners = Threads.nthreads()
        div, rem = divrem(nrows, runners)
        # Raw file coordinates
        startrow(r) = (r-1)*div + min(r-1, rem) + 1
        nrows_on_thread(r) = div + (r <= rem ? 1 : 0)
        endrow(r) = startrow(r) + nrows_on_thread(r) - 1
        # Data coordinates
        startdatarow(r) = startrow(r)

    # Parallel-read
    for r in 1:runners
        # Thread range
        sr, er = startrow(r), endrow(r)
        firstsep = findfirst(sr .<= first.(separator_ranges))
        start_R, start_t = 
        # Destination block
        dest = @view data[sr:er, :]
    end

end