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
"""
    get_scattered_intensity(problem::T, atoms_states, sensors::AbstractArray) where {T<:MeanFieldProblem}

atoms_states: vcat(β, z)
sensor: [ [sensor1], [sensor2], ... ,[sensorN]  ]

**important**: this function does not handles single sensor.
If needed, create a dummy position or duplicate the sensor - then, ignore the one that you don't need.
sensor = [ [sensor1], [sensor_dummy]  ] or [ [sensor1], [sensor1]  ]

"""
function get_scattered_intensity(problem::T, atoms_states, sensors::AbstractArray) where {T<:MeanFieldProblem}
    N = problem.atoms.N
    β = view(atoms_states, 1:N)
    z = view(atoms_states, (N+1):2N)    
    r = problem.atoms.r
    number_configurations = ((N^2)÷2 - N÷2)

    βₙₘ = Array{ComplexF64}(undef, number_configurations)
    cont = 1
    for n=1:N
        for m=(n+1):N            
            βₙₘ[cont] = conj(β[n])*β[m]
            cont += 1
        end
    end
    vβₙₘ = view(βₙₘ, :)
    
    rₙₘ = Array{ComplexF64}(undef, 3, number_configurations)
    cont = 1
    for n=1:N
        r_n = r[n]
        for m=(n+1):N
            rₙₘ[:,cont] = r_n - r[m]
            cont += 1
        end
    end
    vrₙₘ = view(rₙₘ,:, :)
    
    #=
        We don't need to compute each sensor in parallel, because the Folds.mapreduce
        is already doing an excelent job with multi-threading.
    =#
    n_sensors = CoupledDipole.get_number_sensors(sensors)
    intensities = Float64[]
    for i=1:n_sensors
        intensity = ComplexF64(0)
        n_hat = sensors[i]./norm(sensors[i])
        intensity = Folds.mapreduce(+, 1:number_configurations) do k
            (
                begin 
                dot_n_r = n_hat[1]*vrₙₘ[1, k] + n_hat[2]*vrₙₘ[2, k] + n_hat[3]*vrₙₘ[3, k];
                vβₙₘ[k]*cis( k₀*dot_n_r) 
                end 
            )
        end
        intensity +=  sum( (1 .+ z)./4 )
        push!(intensities, 2real(intensity)   )
    end

    return intensities
end

function get_scattered_intensity(problem::T, atoms_states, θ::Number; Gₙₘ=nothing) where {T<:MeanFieldProblem} 
    if isnothing(Gₙₘ)
        Gₙₘ = view(get_geometric_factor(problem.atoms, θ),:,:)
    end

    β = view(atoms_states.β,:)
    z = view(atoms_states.z,:)

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
