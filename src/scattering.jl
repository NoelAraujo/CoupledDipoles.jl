function get_intensities_over_sensors(problem::LinearOptics, β::AbstractArray, all_sensors::AbstractMatrix)
	@debug "start : get intensities over sensors - LinearOptics"

    n_sensors = size(all_sensors,2)
    intensities = zeros(n_sensors)
	
	v_r = view(problem.atoms.r, :, :)
	
	Threads.@threads for i = 1:n_sensors
		one_sensors = view(all_sensors, :, i)
		intensities[i] = _get_intensity_over_sensor(problem.atoms.shape, problem.laser, v_r, one_sensors, β)
    end

	@debug "end  : get intensities over sensors - LinearOptics"
    return intensities
end
#= WORKS ONLY FOR 3D CLOUD DISTRIBUTIONS =#
function _get_intensity_over_sensor(shape::T, laser::Laser, atoms::AbstractMatrix, sensor::AbstractArray,  β::AbstractArray) where T <: ThreeD
	## Laser Pump
	E_L = (im/Γ)*apply_laser_over_oneAtom(laser, sensor)
    
    # ## Scattered
	E_scatt = zero(ComplexF64)
	n̂ = sensor/norm(sensor)
	
    dot_n_r = zero(ComplexF64)
    j = 1
	for atom in eachcol(atoms)
        dot_n_r = n̂[1]*atom[1] + n̂[2]*atom[2] + n̂[3]*atom[3]
        dot_n_r = cis(-k₀*dot_n_r)
		E_scatt += dot_n_r*β[j]
        j += 1
	end
    
	ikr = im*k₀*norm(sensor)
	E_scatt = E_scatt*exp(ikr)/ikr
	return abs2(E_L + E_scatt)
end


# function get_intensity_over_angle(problem::SimulationScalar, atoms_states, θ::Number; Gₙₘ=nothing)
#     if isnothing(Gₙₘ)
#         Gₙₘ = get_geometric_factor(problem.atoms, θ)
#     end

#     β = atoms_states
    
#     βₙₘ = transpose(β*β')
#     βₙₘ[diagind(βₙₘ)] .= abs2.(β)

#     intensity = real(sum(βₙₘ.*Gₙₘ))
# end




# ### --------------- MEAN FIELD ---------------
"""
	get_intensities_over_sensors(problem::NonLinearOptics{MeanField}, atoms_states::AbstractArray, sensors::AbstractArray

atoms_states: vcat(β, z)
sensor: [ [sensor1], [sensor2], ... ,[sensorN]  ]

**important**: check benchmark folders to understand how this code was produced
"""
function get_intensities_over_sensors(problem::NonLinearOptics{MeanField}, atoms_states::AbstractArray, sensors::AbstractArray)
	@debug "start : get intensities over sensors - NonLinearOptics"

    N, r = problem.atoms.N, problem.atoms.r
    β, z = view(atoms_states, 1:N), view(atoms_states, (N+1):2N)    

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
        r_n = @view(r[:,n])
        for m=(n+1):N
            rₙₘ[1,cont] = r_n[1] - r[1,m]
            rₙₘ[2,cont] = r_n[2] - r[2,m]
            rₙₘ[3,cont] = r_n[3] - r[3,m]
            cont += 1
        end
    end
    vrₙₘ = view(rₙₘ,:, :)
    
    
    n_sensors = size(sensors, 2)
    intensities = Float64[]
    
    if n_sensors==1
        n_hat = sensors./norm(sensors)
        intensity = _oneSensor_MeanField_scattering(n_hat, vβₙₘ, vrₙₘ, number_configurations)
        intensity +=  sum( (1 .+ z)./4 )

        push!(intensities, 2real(intensity)   )
    else
        const_sum = sum( (1 .+ z)./4 )
        for oneSensor in eachcol(sensors)
            n_hat =  oneSensor./norm(oneSensor)
            
            intensity = _manySensors_MeanField_scattering(n_hat, vβₙₘ, vrₙₘ, number_configurations)
            intensity +=  const_sum

            push!(intensities, 2real(intensity)   )
        end
    end

	@debug "end  : get intensities over sensors - NonLinearOptics"
    return intensities
end
function _oneSensor_MeanField_scattering(n_hat, vβₙₘ, vrₙₘ, number_configurations; k₀=1)
    intensity = ComplexF64(0)
    for cont = 1:number_configurations
        dot_n_r = n_hat[1]*vrₙₘ[1, cont] + n_hat[2]*vrₙₘ[2, cont] + n_hat[3]*vrₙₘ[3, cont]
        intensity += vβₙₘ[cont]*cis( k₀*dot_n_r  )
    end    
    return intensity
end

function _manySensors_MeanField_scattering(n_hat, vβₙₘ, vrₙₘ, number_configurations; k₀=1)
    intensity = ComplexF64(0)
    intensity = Folds.mapreduce(+, 1:number_configurations) do k
        (
            begin 
                @inbounds vβₙₘ[k]*cis( k₀*(n_hat[1]*vrₙₘ[1, k] + n_hat[2]*vrₙₘ[2, k] + n_hat[3]*vrₙₘ[3, k])) 
            end 
        )
    end
    return intensity
end

# function get_scattered_intensity(problem::T, atoms_states, θ::Number; Gₙₘ=nothing) where {T<:MeanFieldProblem} 
#     if isnothing(Gₙₘ)
#         Gₙₘ = view(get_geometric_factor(problem.atoms, θ),:,:)
#     end

#     N = problem.atoms.N
#     β = view(atoms_states, 1:N)
#     z = view(atoms_states, (N+1):2N)

#     βₙₘ = transpose(β*β')
#     βₙₘ[diagind(βₙₘ)] .= (1 .+ z)./2

#     intensity = real(sum(βₙₘ.*Gₙₘ))
#     return intensity
# end

# function get_geometric_factor(atoms, Θ)
#     N = atoms.N
#     r = atoms.r

#     Gₙₘ = zeros(ComplexF64, N,N)
#     for n=1:N
#         for m=n:N 
#             rx = r[n][1] - r[m][1]
#             ry = r[n][2] - r[m][2]
#             rz = r[n][3] - r[m][3]
            
#             argument_bessel = sqrt( Complex( -k₀^2*(sin(Θ)^2)*(rx^2 +ry^2)) )
#             Gₙₘ[n,m] = cis(k₀*rz*cos(Θ))*besseli(0, argument_bessel )
#         end
#     end
#     return Hermitian(Gₙₘ)
# end
