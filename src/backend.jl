# =========================
# ==== Reading Backend ====
# =========================

# StagData constructor (parameters.dat)
function aggregate_StagData(Sname::String, filename::String, timevec::Union{Array{Float64,1}, Nothing})
    # Parameter file
    filename = filename * "_parameters.dat"
    @assert isfile(filename) "$filename not found"
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
    tend    = isnothing(timevec) ? nothing : sec2Gyr*timevec[end]
    ndts    = isnothing(timevec) ? nothing : length(timevec)
    totH₂O  = pass_field("INITIALOCEANMASS")
    # Simulation Parameters
    T_tracked   = pass_field("TRACERS_TEMPERATURE")
    H₂O_tracked = pass_field("TRACEELEMENT_WATER")
    Crb_tracked = pass_field("TRACEELEMENT_CARBON")

    return StagData(name, shape, rkm, rcmb, nz, tend, ndts, totH₂O,
                        T_tracked, H₂O_tracked, Crb_tracked)
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
        isnothing(Cmass_idx) && (rm_Cmass_error = false) # Disable if not present
    # --- Compute number of rows and Character range for columns
        nrows = (length(map) - lf) ÷ line_bytes

    # Preallocate data array
        data = Array{Float64}(undef, nrows, length(header))

    # Chunk configuration
        runners = Threads.nthreads()
        div, rem = divrem(nrows, runners)

    # Parallel-read
        Threads.@threads :static for r in 1:runners
            # Thread range
            sr, er = startrow(r, div, rem), endrow(r, div, rem)
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
                    if rm_Cmass_error && j==Cmass_idx
                        dest[i-sr+1, j] = 0.0
                    else
                        fld = @view rowstr[rcol]
                        dest[i-sr+1, j] = Parsers.parse(Float64, fld)
                    end
                    j+=1
                end
            end
        end
    
    return data, String.(header)
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

    # Parallel-read
    Threads.@threads :static for r in 1:runners
        # Thread range
        sb, eb = startrow(r, div, rem), endrow(r, div, rem)
        # Destination view array
        dest = @view data[:, sb:eb, :]
        # Block iterator
        @inbounds for block in sb:eb
            # Block range
            block_range = block_to_byte_range(map, block, lf, Stag.nz, line_bytes, separator_bytes)
            # Time entry
            time_range = first(block_range):first(block_range)+separator_bytes
            t = tryparse(Float64, split(String(map[time_range]))[6])
            @assert !isnothing(t) "Failed to parse time entry in rprof.dat"
            time[block] = t
            # Row iterator
            @inbounds for row in 1:Stag.nz
                row_range = row_to_byte_range(row, first(block_range), line_bytes, separator_bytes)
                fld = split(String(@view map[row_range]))
                dest[row, block-sb+1, :] .= Parsers.parse.(Float64, fld)
            end
        end
    end
    return data, String.(header), time

end

# Backend for reading StagYY "plates_analyse.dat"s
function read_StagYY_platesfile(filename::String)
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
        header = replace.(header, "#" => "")
    # --- Compute number of rows and Character range for columns
        nrows = (length(map) - lf) ÷ line_bytes

    # Preallocate data array
        data = Array{Float64}(undef, nrows, length(header))

    # Chunk configuration
        runners = Threads.nthreads()
        div, rem = divrem(nrows, runners)

    # Parallel-read
        Threads.@threads :static for r in 1:runners
            # Thread range
            sr, er = startrow(r, div, rem), endrow(r, div, rem)
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
                    fld = @view rowstr[rcol]
                    dest[i-sr+1, j] = Parsers.parse(Float64, fld)
                    j+=1
                end
            end
        end
    return data, String.(header)
end

# High-level output reader
function load_sim(sroot::String, Sname::String; time::Bool=true, rprof::Bool=true, plates::Bool=true)
    @assert time || rprof || plates "No data block requested."
    # Filenames
    Tfname = sroot * "_time.dat"
    Rfname = sroot * "_rprof.dat"
    Pfname = sroot * "_plates_analyse.dat"
    # Check files + construction
    if time; @assert isfile(Tfname) "$Tfname not found"; time_data, time_header = read_StagYY_timefile(Tfname); end
    if plates; @assert isfile(Pfname) "$Pfname not found"; plates_data, plates_header = read_StagYY_platesfile(Pfname); end
    tvec = time ? time_data[:,findfirst(time_header.=="time")] : plates ? plates_data[:,findfirst(plates_header.=="time")] : nothing
    Stag = aggregate_StagData(Sname, sroot, tvec)
    if rprof; @assert isfile(Rfname) "$Rfname not found"; rprof_data, rprof_header, rprof_time = read_StagYY_rproffile(Rfname, Stag); end
    # Encoding
    idxT, idxR, idxP = data_encoding(time ? time_header : nothing, rprof ? rprof_header : nothing, plates ? plates_header : nothing)
    # Unit conversions
    if time
        time_data[:,idxT["time"]] .= sec2Gyr*time_data[:,idxT["time"]]; # Convert time to Gyr
        time_data[:,idxT["Vrms"]] .= m_s2cm_yr*time_data[:,idxT["Vrms"]]; # Convert time to Gyr
        # Outgassing
        if Stag.H₂O_tracked
            time_header = vcat(time_header, "OutgassedH2O")
            time_data = hcat(time_data, (time_data[:,idxT["OutputtedNotEruptedH2O"]].+time_data[:,idxT["EruptedH2O"]].+
                                (("SaturationOutgasH2O" in time_header) ? time_data[:,idxT["SaturationOutgasH2O"]] : time_data[:,idxT["SaturationOutgassH2O"]])))
        end
    end
    if rprof
        rprof_data[:,:,idxR["vrms"]] .= m_s2cm_yr*rprof_data[:,:,idxR["vrms"]]
        rprof_header[idxR["vrms"]] = "Vrms"
        if haskey(idxR, "rhomean")
            rprof_header = vcat(rprof_header, "dM")
            rprof_data = cat(rprof_data, ∂M(rprof_data[:, :, idxR["r"]], Stag.rcmb, rprof_data[:, :, idxR["rhomean"]]), dims=3)
        end
    end
    if plates
        plates_header[idxP["Vsurf_rms"]] = "Vsurf"
        plates_data[:,idxP["time"]] .= sec2Gyr*plates_data[:,idxP["time"]]
        plates_data[:,idxP["Vrms"]] .= m_s2cm_yr*plates_data[:,idxP["Vrms"]]
        plates_data[:,idxP["Vsurf_rms"]] .= m_s2cm_yr*plates_data[:,idxP["Vsurf_rms"]]
    end

    if time && rprof && haskey(idxR, "rhomean")
        # Compensates a current bug for ocean mass not being reduced after interior hydration (remove after fix)
        Δ = ((time_data[1,idxT["SurfOceanMass3D"]] + sum(rprof_data[:, 1, idxR["Water"]].*rprof_data[:, 1, end].*1e-2)) - Stag.totH₂O)
        time_data[:,idxT["SurfOceanMass3D"]] .-= Δ
    end
    return DataBlock(Sname, 
                    time ? time_header : nothing, rprof ? rprof_header : nothing, plates ? plates_header : nothing,
                    time ? time_data : nothing, rprof ? rprof_data : nothing, plates ? plates_data : nothing, 
                    plates ? plates_data[:,2] : rprof ? sec2Gyr*rprof_time : nothing,
                    haskey(idxR, "r") ? ∂V(rprof_data[:, 1, idxR["r"]], Stag.rcmb) : nothing,
                    Stag)
end

# ==================
# ==== Encoding ====
# ==================

# Create dictionary encoding for time and rprof blocks
function data_encoding(time_header, rprof_header, plates_header)
    # Initialize
    time_encoding, rprof_encoding, plates_encoding = nothing, nothing, nothing
    # Time header encoding
    if !isnothing(time_header)
        time_encoding = Dict{String,Int64}()
        [time_encoding[name] = i for (i, name) in enumerate(time_header)]
    end
    # Rprof header encoding
    if !isnothing(rprof_header)
        rprof_encoding = Dict{String,Int64}()
        [rprof_encoding[name] = i for (i, name) in enumerate(rprof_header)]
    end
    # Plates header encoding
    if !isnothing(plates_header)
        plates_encoding = Dict{String,Int64}()
        [plates_encoding[name] = i for (i, name) in enumerate(plates_header)]
    end
    return time_encoding, rprof_encoding, plates_encoding
end

"""
  Assess DataBlock header indexing for time, rprof and plates data.
    
    idxT, idxR, idxP = data_encoding(D::DataBlock)
"""
data_encoding(D::DataBlock) = data_encoding(D.timeheader, D.rprofheader, D.platesheader)

# ===========================
# ========== VTK ============
# ===========================

function readVTK(fname::String)

    # LightXML is not able to read the XML file coming out of StagYY utilities. Thus the workaround is:
    #   StructuredGrid format append all binary data at the bottom, and just specify the byte offsets in the XML tags.
    #   Therefore we read the offsets using LightXML, then read the binary part manually.

    # Helper function. Finds the beginning of each tag, errors if not found.
    function first_tag(node, tag)
        els = LightXML.get_elements_by_tagname(node, tag)
        length(els) == 0 && error("tag $tag not found")
        return els[1]
    end

    # 1. read the file
    raw = Mmap.mmap(fname)

    # 2. Find the data start
    head = String(raw[1:10_000]) # Restrict exploration to first 10k bytes, as XML part is always at the top.
    rng = findfirst("<AppendedData", head) # Find <AppendedData ...> tag
    rng === nothing && error("no <AppendedData> found") # Check
    app_tag_start = first(rng) # findfirst gives unitrange, we need the very first byte

    # 3. LightXML crashes if the binary part is integrated. Thus we truncate the XML string before the binary part and read the fake file.
    xml_str = String(raw[1:app_tag_start-1]) * "</VTKFile>"
    doc  = parse_string(xml_str)
    rootnode = LightXML.root(doc)

    # === 4. Find CellData and Points tags ===
    sgrid    = first_tag(rootnode, "StructuredGrid")
    piece    = first_tag(sgrid, "Piece")
    celldata = first_tag(piece, "CellData")  # your file clearly has CellData
    points_tag = first_tag(piece, "Points")  # Point information

    # === 5. find where the REAL binary starts (the '_' after <AppendedData ...>) ===
    sub = raw[app_tag_start:end]
    u = findfirst(x -> x == UInt8('_'), sub)
    u === nothing && error("no '_' after <AppendedData>")
    bin_start = app_tag_start + u             # first byte of actual appended data

    # === 5.2 Extract Point data ===
    da_pts = LightXML.get_elements_by_tagname(points_tag, "DataArray")[1]
    off  = parse(Int, LightXML.attribute(da_pts, "offset"))
    pos  = bin_start + off
    nbytes = reinterpret(UInt32, raw[pos:pos+3])[1]
    pos += 4
    buf = raw[pos:pos+nbytes-1]
    points_flat = reinterpret(Float32, buf)
    points = Float32.(reshape(points_flat, 3, :))

    # === 6. read each CellData array (UInt32 len + Float32 data) ===
    cells = Dict{String,Array{Float32}}()
    for da in LightXML.get_elements_by_tagname(celldata, "DataArray")
        name = LightXML.attribute(da, "Name")
        off  = parse(Int, LightXML.attribute(da, "offset"))

        pos = bin_start + off

        # 4-byte little-endian length prefix (your file uses this)
        nbytes = reinterpret(UInt32, raw[pos:pos+3])[1]
        pos += 4

        buf  = raw[pos:pos+nbytes-1]
        vals = reinterpret(Float32, buf)

        cells[name] = vals
    end
    # Add points to the dictionary
    cells["Points"] = points

    return cells
end


# ======================
# ==== Auxilliaries ====
# ======================
    # Sequential row/block <-> byte-finder for chunking
        startrow(r, div, rem) = (r-1)*div + min(r-1, rem) + 1   # Return first row on time.dat or plates_analyse.dat; block on rprof.dat
        nrows_on_thread(r, div, rem) = div + (r <= rem ? 1 : 0) #
        endrow(r, div, rem) = startrow(r, div, rem) + nrows_on_thread(r, div, rem) - 1 # Return last row on time.dat or plates_analyse.dat; block on rprof.dat
    # Irregular block <-> byte range finder (rprof.dat)
        startbyte(b, header_bytes, nz, line_bytes, separator_bytes) = header_bytes + (b-1)*(nz*line_bytes+separator_bytes) + 1
        function endbyte(map, b, header_bytes, nz, line_bytes, separator_bytes) 
            e = (startbyte(b, header_bytes, nz, line_bytes, separator_bytes) + (nz*line_bytes+separator_bytes) - 1)
            (map[e]==0x0A) && (e -= 1)
            (e>=startbyte(b, header_bytes, nz, line_bytes, separator_bytes) && map[e]==0x0D) && (e -= 1)
            return e
        end
        block_to_byte_range(map, b, header_bytes, nz, line_bytes, separator_bytes) = startbyte(b, header_bytes, nz, line_bytes, separator_bytes) : endbyte(map, b, header_bytes, nz, line_bytes, separator_bytes)
        function row_to_byte_range(row, block_start_byte, line_bytes, separator_bytes)
            Δ = separator_bytes + (row-1)*line_bytes
            return (block_start_byte + Δ) : (block_start_byte + Δ + line_bytes - 2)
        end

    # Shell volume/mass differentials
       ∂V(r, rcmb) = 4π/3 * diff(vcat(1e3rcmb, r.+1e3rcmb).^3)
       ∂M(r, rcmb, ρ) = 4π/3 * diff(cat(1e3rcmb, r.+1e3rcmb, dims=1).^3, dims=1) .* ρ

    # Get phase transitions radial coordinates
        function idx_ph_transitions(Dblock::DataBlock; timeidx::Int=1)
            idxR = data_encoding(Dblock)[2]
            r = Dblock.rprofdata[:, timeidx, idxR["r"]]
            midpoint_r = 0.5(r[1:end-1] .+ r[2:end])
            ρ = Dblock.rprofdata[:, timeidx, idxR["rhomean"]]
            ∂ρ_∂r = diff(ρ)./diff(r)
            ∂²ρ_∂r² = diff(∂ρ_∂r)./diff(midpoint_r)
            # Find radial phase boundaries
            baseline = maximum(∂²ρ_∂r²[1:findfirst(r.>=250000.0)]) # Corresponds to the decaying ppv -> pv transition from bottom to top
            shifted_∂²ρ_∂r² = ∂²ρ_∂r² .- baseline
            ph_bounds, prev, current_ph = zeros(Int, 3), -1.0, 1
            for i in eachindex(shifted_∂²ρ_∂r²)
                # Entering peak
                ((shifted_∂²ρ_∂r²[i] <= 0.0) && (prev < 0.0)) && continue
                # Forward step if next value is higher
                if (shifted_∂²ρ_∂r²[i] > prev); (prev = shifted_∂²ρ_∂r²[i]); continue; end
                # Descend until closes to zero, take first if end of array
                prev = ∂²ρ_∂r²[i];
                if (abs(∂²ρ_∂r²[i]) < prev); (prev = abs(∂²ρ_∂r²[i])); continue; end
                # Record boundary index
                ph_bounds[current_ph] = i-1
                current_ph += 1
                prev = -1.0
            end
            (ph_bounds[3]==0) && (ph_bounds[3] = length(shifted_∂²ρ_∂r²)) # If 3rd boundary not found, set to last index
            return ph_bounds
        end

    # Vectorized to reshaped indexing
        @inline function get_ip_it_nrows(nP, i)
            iT = Int(ceil(i/nP)); iP = i - (iT-1)*nP;
            return iP, iT
        end

        @inline function get_ip_it_ncols(nT, i)
            iP = Int(ceil(i/nT)); iT = i - (iP-1)*nT;
            return iP, iT
        end

    # Right axis for plots
        function second_axis!(fig, fpos, x, y, field, ylabelsize, yticklabelsize, ylabelpadding, yreversed, timeplot, trange, locals)

            function transform_data(y, field)
                yp = [first(y), last(y)]
                (field=="SurfOceanMass3D") && (return [locals[1], locals[2]]./om, L"Ocean\;masses") # 10^21 kg
                (field=="Pressure") && (return reverse(km2GPa.(yp)), L"Pressure\;[GPa]") # GPa
                return nothing, nothing
            end
            yp, ylab = transform_data(y, field); (isnothing(yp)) && return
            ax2 = Axis(fig[timeplot ? fpos[1] : fpos[1]+1, fpos[2]], yticklabelcolor = :black, yaxisposition = :right, ylabel = ylab, ylabelsize=ylabelsize, yticklabelsize=yticklabelsize, ylabelpadding=ylabelpadding,
                        xgridvisible=false, ygridvisible=false, yreversed=yreversed)
            scatter!(ax2, first(x)*ones(2), yp, alpha=0.0)
            # yreversed && ylims!(ax2, maximum(yp), minimum(yp))
            (timeplot || field=="Pressure") && (yreversed ? ylims!(ax2, reverse(yp)) : (ylims!(ax2, yp)))
            xlims!(ax2, trange)
            hidespines!(ax2); hidexdecorations!(ax2)
        end

    # Convert depth in km to pressure in GPa using PREM
        function km2GPa(depth)

            r = 6371 - max(depth,0.0)
            r < 0.0 ? error("Depth beyond Earth's radius") : nothing
            #  Inner core
            if (r >= 0.0 && r <= 1221.5)
            pressure = -575.170316 + 0.29421732*depth -2.30448727e-05*depth*depth
            
            #  Outer core
            elseif (r >= 1221.5 && r < 3480.0)
            pressure = -293.862745 +0.183037581*depth  -1.202957e-05*depth*depth
            
            #  Lower mantle
            elseif (r >= 3480.0 && r <= 3630.0)
            pressure = -33.49617 + 0.05853999*depth
            elseif (r >= 3630.0 && r <= 5600.0)
            pressure = -3.56294332+0.0390269653*depth+3.11829162e-06*depth*depth
            elseif (r >= 5600.0 && r <= 5701.0)
            pressure =  -5.758207 + 0.04415793*depth
            
            # Transition zone
            elseif (r >= 5701.0 && r <= 5771.0)
            pressure =  -2.875294 + 0.03985716*depth
            elseif (r >= 5771.0 && r <= 5971.0)
            pressure =  -2.056346 + 0.03845555*depth
            elseif (r >= 5971.0 && r <= 6151.0)
            pressure = -0.5319232 + 0.03467901*depth
            
            # LVZ & LID
            elseif (r >= 6151.0 && r <= 6346.6)
            pressure =   -0.213082 + 0.0332782*depth
            elseif (r >= 6346.6 && r <= 6371.0)
            pressure = 2.6e3 * 9.8 * (6371.0 - r ) *1000/1.0e9 # rho=2.6g/cm3, in GPa
            end
            
            return pressure
        end

    # Convert pressure in GPa to depth in km using PREM
        function GPa2km(pressure)
            depth = LinRange(0., 6371., 800)
            P = sort(km2GPa.(depth))
            itp = interpolate((P,), depth, Gridded(Linear()))
            return itp(pressure)
        end