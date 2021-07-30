function apply_laser_over_atoms(laser::Laser{PlaneWave3D}, atoms)
    @debug "start: apply_laser_over_atoms - Laser{PlaneWave3D}"
    
    s, Δ = laser.s, laser.Δ
    
    v_atoms = view(Float64.(atoms.r), :, :)
    direction = view(laser.pump.direction, :)
    E = zeros(Complex{eltype(atoms.r)}, atoms.N)
    E₀ = estimate_E₀(laser)
    
    Threads.@threads for n in 1:(atoms.N)
        oneAtom = v_atoms[:,n]
        @inbounds E[n] = _core_PlaneWave3D(oneAtom, E₀, direction)
    end

    @debug "end  : apply_laser_over_atoms - Laser{PlaneWave3D}"
    return E
end

function apply_laser_over_oneAtom(laser::Laser{PlaneWave3D}, oneAtom::AbstractArray)
    E₀ = estimate_E₀(laser)
    direction = laser.pump.direction
    E = _core_PlaneWave3D(oneAtom, E₀, direction)

    return E
end

"""
    E₀ * exp(+im * k₀ *dot(oneAtom, direction))
"""
function _core_PlaneWave3D(oneAtom, E₀, direction)
    E₀ * cis(+k₀ * (oneAtom[1]*direction[1] + oneAtom[2]*direction[2] + oneAtom[3]*direction[3]))
end


