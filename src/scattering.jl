function scattering_fuction(distance::Symbol, dimension::Symbol)
    if dimension == :ThreeD
        if distance == :nearField
            sf = default_farField3D_Field
        elseif distance == :farField
            sf = default_farField3D_Field
        else
            @error "Only options are `:nearField` and `:farField`"
        end
    else
        @error "Only :ThreeD is supported."
    end
    return sf
end

"""
    default_nearField3D_Field(atoms::AbstractMatrix, β::AbstractArray, sensor::AbstractArray)

- atoms: each column is an atom
- atomic_states: β for Scalar Model, ou [β,z] for Mean Field Model
- sensor: measurement position

Returns a Complex Value of the Eletric Field: +(Γ/2) * ∑ⱼ exp(-i*k₀* n̂⋅R⃗ⱼ)/(k₀* sensor⋅R⃗ⱼ)
- R⃗ⱼ : atom j
"""
function default_nearField3D_Field(atoms::AbstractMatrix, β::AbstractArray, sensor::AbstractArray)
    E_scatt = zero(eltype(β))

    j = 1
    @inbounds for atom in eachcol(atoms)
        d_SensorAtom = sqrt((sensor[1] - atom[1])^2 + (sensor[2] - atom[2])^2 + (sensor[3] - atom[3])^2)
        E_scatt += cis(k₀ * d_SensorAtom) * (β[j] / d_SensorAtom)
        j += 1
    end
    E_scatt = +(Γ / 2) * im * E_scatt
    return E_scatt
end

"""
    default_farField3D_Field(atoms::AbstractMatrix, β::AbstractArray, sensor::AbstractArray)

- atoms: each column is an atom
- atomic_states: β for Scalar Model, ou [β,z] for Mean Field Model
- sensor: measurement position

Returns a Complex Value of the Eletric Field: -(Γ/2) * (exp(ikr) / ikr) * ∑ⱼ exp(-i*k₀* n̂⋅R⃗ⱼ)

- r : distance sensor to origin (r = norm(sensor))
- n̂ : norm of sensor ( n̂ = sensor / norm(sensor) )
- R⃗ⱼ : atom j
"""
function default_farField3D_Field(atoms::AbstractMatrix, β::AbstractArray, sensor::AbstractArray)
    E_scatt = zero(eltype(β))
    n̂ = sensor / norm(sensor)

    dot_n_r = zero(eltype(β))
    j = 1
    @inbounds for atom in eachcol(atoms)
        dot_n_r = n̂[1] * atom[1] + n̂[2] * atom[2] + n̂[3] * atom[3]
        dot_n_r = cis(-k₀ * dot_n_r)
        E_scatt += dot_n_r * β[j]
        j += 1
    end

    ikr = im * k₀ * norm(sensor)
    E_scatt = -(Γ / 2) * E_scatt * exp(ikr) / ikr
    return E_scatt
end

"""
    scattering_intensity(problem, atomic_states, measurement_positions, scattering_func::Function)

    Computes the Total Electric Field (Laser Pump + Scattering), then returns the Intensity
"""
@views function scattering_intensity(problem, atomic_states, measurement_positions, scattering_func::Function; useThreads=true)
    _laser = problem.laser
    _r = problem.atoms.r
    _physics = problem.physic

    _sensors = measurement_positions
    _states = atomic_states
    _func = scattering_func

    n_sensors = _get_number_elements(_sensors)

    if n_sensors == 1
        return _OnePoint_Intensity(_physics, _laser, _r, _sensors, _states, _func)
    else
        scat_int = Array{Float64}(undef, n_sensors)
        if useThreads
            Threads.@threads for i in 1:n_sensors
                scat_int[i] = _OnePoint_Intensity(_physics, _laser, _r, view(_sensors, :, i), _states, _func)
            end
        else
            for i in 1:n_sensors
                scat_int[i] = _OnePoint_Intensity(_physics, _laser, _r, view(_sensors, :, i), _states, _func)
            end
        end
        return scat_int
    end
end


