using StagPP, MAGEMin_C, GLMakie, DelimitedFiles

# Auxilliaries
@inline function get_ip_it_nrows(nP, i)
    iT = Int(ceil(i/nP)); iP = i - (iT-1)*nP;
    return iP, iT
end

function vectorize!(v1, v2, v3, XB, XH, X50, n1, n2, n3, Xv)
    if nH > 1
        # P and T
        vv1, vv2 = repeat(repeat(v1, outer=n2), n3), repeat(repeat(v2, inner=n1), n3)
        # Compositions
        n1n2    = n1*n2
        n1n2n3  = n1*n2*n3
        # Construct possible comp vectors
        Xp      = zeros(3nH, length(XB))
        for i in 1:nH
            Xp[i,:]      .= copy(XH); Xp[i,end] = v3[i]
            Xp[i+nH,:]   .= copy(XB); Xp[i+nH,end] = v3[i]
            Xp[i+2nH,:]  .= copy(X50); Xp[i+2nH,end] = v3[i]
        end
        # Fill first H page
        for i in 1:n1n2n3
            # Key
            key = i<n1n2n3 ? Int(floor(i/n1n2)+1) : nH
            # Assign
            Xv[i]          = Xp[key, :];
            Xv[i+n1n2n3]   = Xp[nH+key, :];
            Xv[i+2n1n2n3]  = Xp[2nH+key, :];
        end
    else
        # P and T
        vv1, vv2  = repeat(v1, outer=n2), repeat(v2, inner=n1)
        # Compositions
        n1n2      = n1*n2
        # Construct possible comp vectors
        Xp        = zeros(2, length(XB))
        Xp[1,:]  .= copy(XH)
        Xp[2,:]  .= copy(XB)
        # Fill first H page
        for i in 1:n1n2
            # Assign
            Xv[i]         = Xp[1, :];
            Xv[i+n1n2]  = Xp[2, :];
        end

    end
    return vcat(vv1, vv1, vv1), vcat(vv2, vv2, vv2)
end

function lunarwater_phases(nP, nT, nH, mantle_comp, H2Orange, phase_out, minimize, display_plot, ncols, Hslice, cmap, interp, db, savein, outHB)
    
    # Vectorized resolution
    nPnT      = nP*nT
    nPnTnH    = nP*nT*nH
    
    # Checks
    (!minimize && isnothing(outHB)) && (minimize=true; println("Turning ON minimize as no output is currently stored"))
    (!minimize && !isnothing(outHB)) && (@assert length(outHB)==3nPnTnH "Current stored output differs from requested resolution.")
    
    # Define P-T conditions
    P       = LinRange(1e-3, 5.0, nP)
    T       = LinRange(680., 2000., nT)
    H       = nH>1 ? LinRange(1e-4H2Orange[1], 1e-4H2Orange[2], nH) : nothing
    
    # Define compositions
    Xox     = ["SiO2", "MgO", "FeO", "CaO", "Al2O3", "Na2O", "Cr2O3", "TiO2", "MnO", "H2O"]
    XB      = [50.47, 7.58, 10.43, 11.39, 14.70, 2.79, 0.0, 0.0, 1.68, 0.184]
    XH      = mantle_comp == "Taylor" ? [43.4, 32., 13., 4.5, 6., 0.09, 0.61, 0.3, 0.15, 0.00] : [46.8, 36., 9.2, 3.06, 3.87, 0.05, 0.44, 0.18, 0.132, 0.0]
    X50     = 0.5XB .+ 0.5XH

    if minimize
        # Vectorize
        Xv      = Vector{Vector{Float64}}(undef, 3nP*nT*nH)
        vP, vT  = vectorize!(P, T, H, XB, XH, X50, nP, nT, nH, Xv)

        # Initialize minimizer
        println("Calculating upper mantle sᴴ²ᴼ...")
        data    = Initialize_MAGEMin(db, verbose=false)
        if isempty(phase_out)
            outHB   = multi_point_minimization(10vP, vT.-273.15, data, X=Xv, Xoxides=Xox, name_solvus=true, sys_in="wt") # kbar and K
        else
            rm_list = remove_phases(phase_out, db)
            outHB   = multi_point_minimization(10vP, vT.-273.15, data, X=Xv, Xoxides=Xox, name_solvus=true, sys_in="wt", rm_list=rm_list) # kbar and K
        end
        Finalize_MAGEMin(data);
    end

    if display_plot=="p"|| display_plot=="P"
        # Figure
        fig = Figure(size=(1500,800));
        @assert (em=="H" || em=="B" || em=="X50")
        st = em=="H" ? 1 : em=="X50" ? 2nPnTnH+1 : nPnTnH+1
        if nH>1
            # Assess closes Hslice
            idx = argmin(abs.(H.-1e-4Hslice))
            st+=(idx-1)*nPnT
        end
        # Phase iterator
        observed_phases = Dict{String, Array{Float64, 2}}()
        for n in st:nPnT+(st-1)
            # Vectorized to mesh indexing
            li, lj = get_ip_it_nrows(nP, n-(st-1))
            # Phase iterator
            for (idx, phase) in enumerate(outHB[n].ph)
                # Initialise phase if not present
                !haskey(observed_phases, phase) && (observed_phases[phase] = zeros(Float64, nP, nT))
                # Assign solution phase value
                observed_phases[phase][li, lj] = outHB[n].ph_frac_wt[idx]
            end
        end
        # plot
        xlabsz, ylabsz, titlesz, xticklabsz, yticklabsz, xticksz, yticksz = 20, 20, 22, 16, 16, 12, 12
        i=1; j=1
        for (ph, map) in observed_phases
            ax = Axis(fig[i, j], title=ph, ylabel=L"Pressure\;[GPa]", xlabel=L"Temperature\;[K]", xlabelsize=xlabsz, ylabelsize=ylabsz,
                        xticklabelsize=xticklabsz, yticklabelsize=yticklabsz, xticksize=xticksz, yticksize=yticksz, titlesize=titlesz, yreversed=true)
            hm = heatmap!(ax, T, P, map', colormap=cmap, interpolate=interp)
            Colorbar(fig[i, j+1], hm; labelsize=20, ticklabelsize=15, vertical=true, labelpadding=13)
            j+=2; (j==2ncols+1) && (i+=1; j=1)
        end
        savein!="" && GLMakie.save(savein, fig)
        display(fig)
    end
    return outHB, true, P, T, nH>1 ? H : nothing
end

function sH2O_integration(output, nP, nT, nH, proceed)

    # Check
    @assert proceed "Integration not authorized."

    # Utilities
    nPnT    = nP*nT
    nPnTnH  = nP*nT*nH

    # NAMs solubility assemble
    min_s = min_sᴴ²ᴼ_assembler(Int64(ceil(1.1max(nP, nT))))

    # Vectorized maps
    Hmap, Xmap, Bmap = zeros(Float64, nPnTnH), zeros(Float64, nPnTnH), zeros(Float64, nPnTnH)

    # Integrate hydrous + anhydrous phase contribution
    Threads.@threads for i in 1:nPnTnH
        # ------- Hydrous --------
            # Depleted (Harzburgite) or custom
            H₂O = 0.0; for j in eachindex(output[i].SS_vec)
                ("H2O" in output[i].SS_vec[j].emNames) && continue
                H₂O += 1e2output[i].SS_vec[j].Comp_wt[end-1]*output[i].ph_frac_wt[j]
            end; Hmap[i]=H₂O
            # Enriched (Basalt)
            H₂O = 0.0; for j in eachindex(output[i+nPnTnH].SS_vec)
                ("H2O" in output[i+nPnTnH].SS_vec[j].emNames) && continue
                H₂O += 1e2output[i+nPnTnH].SS_vec[j].Comp_wt[end-1]*output[i+nPnTnH].ph_frac_wt[j]
            end; Bmap[i] = H₂O
            # Half-enriched (Basalt)
            H₂O = 0.0; for j in eachindex(output[i+2nPnTnH].SS_vec)
                ("H2O" in output[i+2nPnTnH].SS_vec[j].emNames) && continue
                H₂O += 1e2output[i+2nPnTnH].SS_vec[j].Comp_wt[end-1]*output[i+2nPnTnH].ph_frac_wt[j]
            end; Xmap[i] = H₂O

        # ------- Anhydrous --------
            idx = 1
            lt = output[i].T_C + 273.15 # K
            lp = 0.1output[i].P_kbar # GPa
            # Filler Harz
            ("ol"   in output[i].ph) && (idx=findfirst(output[i].ph.=="ol");        Hmap[i]+=output[i].ph_frac[idx]*min_s.ol(lp, lt))
            ("opx"  in output[i].ph) && (idx=findfirst(output[i].ph.=="opx");       Hmap[i]+=output[i].ph_frac[idx]*min_s.opx(lp, lt) + min(output[i].SS_vec[idx].Comp[2]/0.1005, 1.0)*min_s.opx_al(lp, lt))
            ("g"    in output[i].ph) && (idx=findfirst(output[i].ph.=="g");         Hmap[i]+=output[i].ph_frac[idx]*min_s.grt(lp, lt))
            ("coe"  in output[i].ph) && (idx=findfirst(output[i].ph.=="coe");       Hmap[i]+=output[i].ph_frac[idx]*min_s.coe(lp,lt))
            ("crst" in output[i].ph) && (idx=findfirst(output[i].ph.=="crst");      Hmap[i]+=output[i].ph_frac[idx]*min_s.crst)
            # Filler MORB
            ("ol"   in output[i+nPnTnH].ph) && (idx=findfirst(output[i+nPnTnH].ph.=="ol");  Bmap[i]+=output[i+nPnTnH].ph_frac[idx]*min_s.ol(lp, lt))
            ("opx"  in output[i+nPnTnH].ph) && (idx=findfirst(output[i+nPnTnH].ph.=="opx"); Bmap[i]+=output[i+nPnTnH].ph_frac[idx]*min_s.opx(lp, lt) + min(output[i+nPnTnH].SS_vec[idx].Comp[2]/0.1005, 1.0)*min_s.opx_al(lp, lt))
            ("g"    in output[i+nPnTnH].ph) && (idx=findfirst(output[i+nPnTnH].ph.=="g");   Bmap[i]+=output[i+nPnTnH].ph_frac[idx]*min_s.grt(lp, lt))
            ("coe"  in output[i+nPnTnH].ph) && (idx=findfirst(output[i+nPnTnH].ph.=="coe"); Bmap[i]+=output[i+nPnTnH].ph_frac[idx]*min_s.coe(lp,lt))
            ("crst" in output[i+nPnTnH].ph) && (idx=findfirst(output[i+nPnTnH].ph.=="crst");Bmap[i]+=output[i+nPnTnH].ph_frac[idx]*min_s.crst)
            # Filler X50
            ("ol"   in output[i+2nPnTnH].ph) && (idx=findfirst(output[i+2nPnTnH].ph.=="ol");  Xmap[i]+=output[i+2nPnTnH].ph_frac[idx]*min_s.ol(lp, lt))
            ("opx"  in output[i+2nPnTnH].ph) && (idx=findfirst(output[i+2nPnTnH].ph.=="opx"); Xmap[i]+=output[i+2nPnTnH].ph_frac[idx]*min_s.opx(lp, lt) + min(output[i+2nPnTnH].SS_vec[idx].Comp[2]/0.1005, 1.0)*min_s.opx_al(lp, lt))
            ("g"    in output[i+2nPnTnH].ph) && (idx=findfirst(output[i+2nPnTnH].ph.=="g");   Xmap[i]+=output[i+2nPnTnH].ph_frac[idx]*min_s.grt(lp, lt))
            ("coe"  in output[i+2nPnTnH].ph) && (idx=findfirst(output[i+2nPnTnH].ph.=="coe"); Xmap[i]+=output[i+2nPnTnH].ph_frac[idx]*min_s.coe(lp,lt))
            ("crst" in output[i+2nPnTnH].ph) && (idx=findfirst(output[i+2nPnTnH].ph.=="crst");Xmap[i]+=output[i+2nPnTnH].ph_frac[idx]*min_s.crst)
    end

    # Reshape
    return reshape(Hmap, nP, nT, nH), reshape(Xmap, nP, nT, nH), reshape(Bmap, nP, nT, nH)
    
end

function plotit(Hm, Xm, Bm, P, T, oH, display_plot, cmap, interp, logscale, savein, wt)
    
    # Inputs
    sz=(1000,600);
    # Correct H
    H = wt ? oH : 1e4oH
    # Check
    !(display_plot=="s" || display_plot=="S") && return
    # Initialise figure
    fig = Figure(size=sz)
    # Limits
    if wt
        Hmap, Xmap, Bmap = Hm, Xm, Bm
    else
        Hmap, Xmap, Bmap = 1e4.*Hm, 1e4.*Xm, 1e4.*Bm
    end
    if logscale
        Hmap = log10.(Hmap); Xmap = log10.(Xmap); Bmap = log10.(Bmap)
    end
    mx = max(maximum(Hmap), maximum(Bmap), maximum(Xmap))
    mn = min(minimum(Hmap), minimum(Bmap), minimum(Xmap))

    if length(H)==1
        # Axis
        ax1 = Axis(fig[2,1], xlabel=L"Temperature\;\mathrm{[K]}", ylabel=L"Pressure\;\mathrm{[GPa]}", title="Mantle", yreversed=true, titlesize=20, titlefont=:bold_italic)
        ax2 = Axis(fig[2,2], xlabel=L"Temperature\;\mathrm{[K]}", ylabel=L"Pressure\;\mathrm{[GPa]}", title="50%", yreversed=true, titlesize=20, titlefont=:bold_italic)
        ax3 = Axis(fig[2,3], xlabel=L"Temperature\;\mathrm{[K]}", ylabel=L"Pressure\;\mathrm{[GPa]}", title="Basalt", yreversed=true, titlesize=20, titlefont=:bold_italic)
        # Maps
        hmH = heatmap!(ax1, T, P, Hmap[:,:,1]', colormap=cmap, interpolate=interp, colorrange=(mn, mx))
        hmX = heatmap!(ax2, T, P, Xmap[:,:,1]', colormap=cmap, interpolate=interp, colorrange=(mn, mx))
        hmB = heatmap!(ax3, T, P, Bmap[:,:,1]', colormap=cmap, interpolate=interp, colorrange=(mn, mx))
    else
        nP, nT, nH = length(P), length(T), length(H)
        # Axis
        ax1 = Axis3(fig[2,1], xlabel=L"Temperature\;\mathrm{[K]}", ylabel=L"Pressure\;\mathrm{[GPa]}", title="Mantle", zlabel=L"Water\;content\;\mathrm{[μg/g]}", yreversed=true, titlesize=20, titlefont=:bold_italic)
        ax2 = Axis3(fig[2,2], xlabel=L"Temperature\;\mathrm{[K]}", ylabel=L"Pressure\;\mathrm{[GPa]}", title="50%", zlabel=L"Water\;content\;\mathrm{[μg/g]}", yreversed=true, titlesize=20, titlefont=:bold_italic)
        ax3 = Axis3(fig[2,3], xlabel=L"Temperature\;\mathrm{[K]}", ylabel=L"Pressure\;\mathrm{[GPa]}", title="Basalt", zlabel=L"Water\;content\;\mathrm{[μg/g]}", yreversed=true, titlesize=20, titlefont=:bold_italic)
        # Maps
        for h in 1:nH
            hmH = surface!(ax1, T, P, H[h]*ones(nP, nT)', color=Hmap[:,:,h]', colormap=cmap, interpolate=interp, colorrange=(mn, mx))
            hmX = surface!(ax2, T, P, H[h]*ones(nP, nT)', color=Xmap[:,:,h]', colormap=cmap, interpolate=interp, colorrange=(mn, mx))
            hmB = surface!(ax3, T, P, H[h]*ones(nP, nT)', color=Bmap[:,:,h]', colormap=cmap, interpolate=interp, colorrange=(mn, mx))
        end
        # text H
        text!(ax1, 1400., 0.9, H[end], text="Ol + Opx",       color=:white, font=:bold_italic)
        text!(ax1, 1160., 2.7, H[end], text="Ol + Grt",    color=:white, font=:bold_italic)
        text!(ax1, 680.,  0.8, H[end], text="Ol + Opx",  color=:white, font=:bold_italic)
        text!(ax1, 700.,  2.5, H[end], text="Ol + Grt",  color=:white, font=:bold_italic)
        # text B
        text!(ax3, 1600., 0.7, H[end], text="Opx",       color=:white, font=:bold_italic)
        text!(ax3, 680.,  1.5, H[end], text="Anth + ta", color=:white, font=:bold_italic)
        text!(ax3, 1300., 2.1, H[end], text="Grt",    color=:grey40, font=:bold_italic)
        text!(ax3, 1000., 3.3, H[end], text="Coe",       color=:grey40, font=:bold_italic)

    end
    clab = L"Water\;content"
    clab = logscale ? L"%$clab\;\mathrm{[log_{10}" : L"%$clab\;\mathrm{["
    clab = wt ? L"%$clab\;wt\%]}" : L"%$clab\;μg/g]}"

    Colorbar(fig[1,1:3], hmH, vertical=false, label=clab, tellwidth=false, width=Relative(1/2))
    
    # Save
    savein!="" && GLMakie.save(savein, fig)

    # Display
    display(fig)

end

function exportit(Hm, Xm, Bm, P, T, H, savein)

    savein=="" && return

    # Use wt%
    Hmap, Xmap, Bmap, Hw = Hm, Xm, Bm, H

    # Array dimensions
    nP, nT, nR = length(P), length(T), length(H)
    pmap = zeros(Float64, nP, nT)
    ivec = [2, 3, 1, 2, 3, 1, 2, 3, 1]
    # Open file
    open(savein, "w") do io
        # Header
        println(io, nP, " ", nT, " ", nH, "\n");
        println(io, P[1], " ", P[end], " ", T[1], " ", T[end])
        println(io, 2P[end], " ", 2P[end]+nP, " ", 2T[end], " ", 2T[end]+nT)
        println(io, 4P[end], " ", 4P[end]+nP, " ", 4T[end], " ", 4T[end]+nT)
        println(io, Hw[1], " ", Hw[end], "\n")
    
        # UM Data
        for em in 1:3 # MORB -> Harz
            for h in 1:nH # H2O space
                pmap .= em==1 ? Hmap[:,:,h] : em==2 ? Bmap[:,:,h] : Xmap[:,:,h]
                for i in 1:nP
                    for j in 1:nT
                        print(io, pmap[i,j], " ")
                    end; println(io, "")
                end; println(io, "")
            end
        end

        # Rest
        for em in 1:6 # MORB -> Harz
            pmap .= 0.0
            for i in 1:nP
                for j in 1:nT
                    print(io, pmap[i,j], " ")
                end; println(io, "")
            end; println(io, "")
        end
    end
end

# ------------- Inputs --------------
mantle_comp  = "Taylor" # "Warren" | "Taylor"
phase_out    = []
minimize     = true
display_plot = "s"; em = "B"; Hslice=1000
nP, nT, nH   = 40, 40, 5
cmap, interp = :Blues, false
H2Orange     = (10, 2000) # in μg/g
db           = "um"
logscale     = false
wt           = true
savein_plot  = ""#"./scripts/Lunar Water/phases.png"
savein_file  = "" #"./scripts/Lunar Water/StagH2O_moon.dat"
# -----------------------------------

ncols=4
(!@isdefined output) && (output=nothing) # Initialize output variable
output, proceed, P, T, H = lunarwater_phases(nP, nT, nH, mantle_comp, H2Orange, phase_out, minimize, display_plot, ncols, Hslice, cmap, interp, db, savein_plot, output) # Minimization
Hmap, Xmap, Bmap = sH2O_integration(output, nP, nT, nH, proceed) # Integrate
plotit(Hmap, Xmap, Bmap, P, T, H, display_plot, cmap, interp, logscale, savein_plot, wt) # Plot
exportit(Hmap, Xmap, Bmap, P, T, H, savein_file)

# b = solve_sH2O_fO2(100, 100, fO2=false)
# plot_sf(b)

# format:
# full UM basalt maps
# full UM harz maps
# TZ basalt map
# TZ harz map
# LM basalt map
# LM harz map