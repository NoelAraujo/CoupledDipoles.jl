function apply_laser_over_atoms(laser::Laser{Gaussian3D}, atoms)
    @debug "start: apply_laser_over_atoms - Laser{Gaussian3D}"
    
    w₀, s, Δ = laser.pump.w₀, laser.s, laser.Δ

    E = zeros(Complex{eltype(atoms.r)}, atoms.N)
    E₀ = estimate_E₀(laser)
    v_atoms = view(Float64.(atoms.r), :, :)

    Threads.@threads for n in 1:(atoms.N)
        oneAtom = v_atoms[:,n]
        @inbounds E[n] = _core_Gaussian3D(oneAtom, E₀, w₀)
    end

    @debug "end  : apply_laser_over_atoms - Laser{Gaussian3D}"
    return E
end

function apply_laser_over_oneAtom(laser::Laser{Gaussian3D}, position::AbstractArray)
    w₀ = laser.pump.w₀

    E₀ = estimate_E₀(laser)
    E = _core_Gaussian3D(position, E₀, w₀)
    
    return E
end

function _core_Gaussian3D(oneAtom, E₀, w₀)
    x, y, z = oneAtom[1], oneAtom[2], oneAtom[3]

    ## This formula is stable for z==0
    # Ref: Eq 3.11 from "CHAPTER 3. PROPAGATION AND FOCUSING OF OPTICAL FIELDS"
    denominator_factor = 1 .+ 2im .* z / (k₀ * w₀^2)
    Eᵢ = E₀ * cis(+k₀ * z)
    Eᵢ = Eᵢ .* exp.(-(x .^ 2 + y .^ 2) ./ (denominator_factor .* w₀^2))
    Eᵢ = Eᵢ ./ denominator_factor
    
    return Eᵢ
end