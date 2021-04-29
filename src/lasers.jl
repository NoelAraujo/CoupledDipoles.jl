function laser_over_sensors(laser::PlaneWave, sensors)
    direction, s, Δ = laser.direction, laser.s, laser.Δ

    E₀ = estimate_E₀(laser)
    atomic_component = select_sensor_axes(sensors, direction)

    Threads.@threads for n in 1:Nsensors #
        @inbounds E[n] = E₀ * cis(k₀ * atomic_component[n])
    end
    return E
end

function laser_over_sensors(laser::Gaussian3D, sensors)
    ω₀, s, Δ = laser.ω₀, laser.s, laser.Δ
    Nsensors = length(sensors)

    E = zeros(ComplexF64, Nsensors)
    E₀ = estimate_E₀(laser)

    Threads.@threads for n in 1:Nsensors #
        oneSensor = get_one_sensor(sensors, n)
        @inbounds E[n] = _core_Gaussian(oneSensor, E₀, ω₀)
    end
    return E
end

### --------------- PLANE WAVES ---------------
laser_over_atoms(laser::PlaneWave, atoms) = apply_plane_wave_field(laser, atoms)

function apply_plane_wave_field(laser, atoms)
    direction, s, Δ = laser.direction, laser.s, laser.Δ

    E₀ = estimate_E₀(laser)
    atomic_component = select_atoms_axes(atoms, direction)
    return E₀ .* exp.(+im * k₀ * atomic_component)
end

estimate_E₀(laser) = √(laser.s * (1 + 4(laser.Δ / Γ)^2) / 2)
### --------------- GAUSSIANS ---------------
laser_over_atoms(laser::Gaussian3D, atoms) = apply_gaussian3D_field(laser, atoms)

"""
    Gaussian Beam along z direction
"""
function apply_gaussian3D_field(laser, atoms)
    ω₀, s, Δ = laser.ω₀, laser.s, laser.Δ

    E = zeros(ComplexF64, atoms.N)
    E₀ = estimate_E₀(laser)

    Threads.@threads for n in 1:(atoms.N) # 
        oneAtom = get_one_atom(atoms, n)
        @inbounds E[n] = _core_Gaussian(oneAtom, E₀, ω₀)
    end
    return E
end

function _core_Gaussian(oneAtom, E₀, ω₀)
    x, y, z = oneAtom[1], oneAtom[2], oneAtom[3]

    ## This formula is stable for z==0
    # Ref: Eq 3.11 from "CHAPTER 3. PROPAGATION AND FOCUSING OF OPTICAL FIELDS"
    denominator_factor = 1 .+ 2im .* z / (k₀ * ω₀^2)
    Eᵢ = E₀ .* exp.(+im * k₀ * z)
    Eᵢ = Eᵢ .* exp.(-(x .^ 2 + y .^ 2) ./ (denominator_factor .* ω₀^2))
    Eᵢ = Eᵢ ./ denominator_factor
    return Eᵢ
end

function estimate_waist(atoms::Cube)
    return size(atoms) / 4
end
function estimate_waist(atoms::Sphere)
    return size(atoms) / 4
end
