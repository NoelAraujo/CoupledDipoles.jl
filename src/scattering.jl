### --------------- SCALAR ---------------
function get_scattered_intensity(problem::T, atoms_states, sensors) where {T<:Scalar}
    ## Pump Field
    E_L = im * laser_over_sensors(problem.laser, sensors)
    n_sensors = get_number_sensors(sensors)
    intensities = zeros(n_sensors)

    Threads.@threads for n in 1:n_sensors #
        sensor_position = get_one_sensor(sensors, n)
        intensities[n] = get_scattered_light_scalar(
            problem, atoms_states, sensor_position, E_L[n]
        )
    end

    return intensities
end

function get_scattered_light_scalar(problem::T, β, sensor_position, E_L) where {T<:Scalar}
    E_scatt = zero(ComplexF64)
    sensor_versor = sensor_position / norm(sensor_position)
    for n in 1:(problem.atoms.N)
        E_scatt += get_E_scatterd(problem.atoms, sensor_versor, sensor_position, β, n)
    end
    return abs2(E_L + E_scatt)
end

function get_E_scatterd(atoms::T, sensor_versor, sensor_position, β, n) where {T<:ThreeD}
    r_atom = atoms.r[n] # == get_one_atom(atoms, n)

    dot_n_r = cis(-dot(sensor_versor, r_atom))
    E_scatt = dot_n_r * β[n]

    ikr = im * k₀ * norm(sensor_position)
    E_scattFinal = E_scatt * cis(k₀ * norm(sensor_position)) / ikr
    return E_scattFinal
end


function get_scattered_intensity(problem::SimulationScalar, atoms_states, θ::Number; Gₙₘ=nothing)
    if isnothing(Gₙₘ)
        Gₙₘ = get_geometric_factor(problem.atoms, θ)
    end

    β = atoms_states
    
    βₙₘ = transpose(β*β')
    βₙₘ[diagind(βₙₘ)] .= abs2.(β)

    intensity = real(sum(βₙₘ.*Gₙₘ))
end




### --------------- MEAN FIELD ---------------
## TO DO
function get_scattered_intensity(problem::T, atoms_states, sensors::AbstractArray) where {T<:MeanFieldProblem}   
    n_sensors = get_number_sensors(sensors)
    intensities = zeros(n_sensors)

    Threads.@threads for n in 1:n_sensors #
        sensor_position = get_one_sensor(sensors, n)
        intensities[n] = get_intensity_over_point_in_space_MF(
            problem.atoms.r, atoms_states, sensor_position
        )
    end

    return intensities
end

function get_intensity_over_point_in_space_MF(atoms_positions, atoms_states, sensor_position)
    N = length(atoms_states)÷2
    β = view(atoms_states, 1:N)
    z = view(atoms_states, (N+1):2N)

    intensity = zero(eltype(β))
    for m=1:N
        r_m = atoms_positions[m]
        for n = (m+1):N
            r_nm = r_m - atoms_positions[n]
            dot_n_r = cis(dot(sensor_position, r_nm))
            intensity += conj(β[n])*β[m]*dot_n_r
        end
    end   
    for n=1:N
        intensity += (1 + z[n])/4
    end
    return 2real(intensity)
end

function get_scattered_intensity(problem::T, atoms_states, θ::Number; Gₙₘ=nothing) where {T<:MeanFieldProblem} 
    if isnothing(Gₙₘ)
        Gₙₘ = get_geometric_factor(problem.atoms, θ)
    end

    β = atoms_states.β
    z = atoms_states.z

    βₙₘ = transpose(β*β')
    βₙₘ[diagind(βₙₘ)] .= (1 .+ z)./2


    intensity = real(sum(βₙₘ.*Gₙₘ))
    return intensity
end

function get_geometric_factor(atoms, Θ)
    N = atoms.N
    r = atoms.r

    Gₙₘ = zeros(ComplexF64, N,N)
    for n=1:N
        for m=n:N 
            rx = r[n][1] - r[m][1]
            ry = r[n][2] - r[m][2]
            rz = r[n][3] - r[m][3]
            
            argument_bessel = sqrt( Complex( -k₀^2*(sin(Θ)^2)*(rx^2 +ry^2)) )
            Gₙₘ[n,m] = cis(k₀*rz*cos(Θ))*besseli(0, argument_bessel )
        end
    end
    return Hermitian(Gₙₘ)
end
