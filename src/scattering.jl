function get_intensities_over_sensors(problem, β::AbstractArray, all_sensors::AbstractMatrix)
	@debug "start : get intensities over sensors"

    n_sensors = size(all_sensors,2)
    intensities = zeros(n_sensors)
	
	v_r = view(problem.atoms.r, :, :)
	
	Threads.@threads for i = 1:n_sensors
		one_sensors = view(all_sensors, :, i)
		intensities[i] = _get_intensity_over_sensor(problem.atoms.shape, problem.laser, v_r, one_sensors, β)
    end

	@debug "end  : get intensities over sensors"
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

### --------------- SCALAR ---------------
# function get_scattered_intensity(problem::T, atoms_states, sensors) where {T<:Scalar}
#     ## Pump Field
#     E_L = im * laser_over_sensors(problem.laser, sensors)
#     n_sensors = get_number_sensors(sensors)
#     intensities = Float64[]
    
#     ## Probably the code below will change in the future
#     ## It is only in its simplest form to later tests
#     for i in 1:n_sensors #
#         sensor_position = sensors[i]
#         r =  norm(sensor_position)
#         n_hat = sensor_position / r
        

#         E_scatt = Folds.mapreduce(+, 1:problem.atoms.N) do j
#             (
#                 begin 
#                     r_atom = problem.atoms.r[j]
#                     cis(-dot(n_hat, r_atom))*atoms_states[j]
#                 end 
#             )
#         end

    
#         ikr = im * k₀ * r
        
#         intensity = abs2(E_L[i] + (exp(ikr)/ikr)*E_scatt)
#         push!(intensities, intensity)
#     end

#     return intensities
# end

# function get_scattered_light_scalar(problem::T, β, sensor_position, E_L) where {T<:Scalar}
#     E_scatt = zero(ComplexF64)
#     sensor_versor = sensor_position / norm(sensor_position)
#     for n in 1:(problem.atoms.N)
#         E_scatt += get_E_scatterd(problem.atoms, sensor_versor, sensor_position, β, n)
#     end
#     return abs2(E_L + E_scatt)
# end

# function get_E_scatterd(atoms::T, sensor_versor, sensor_position, β, n) where {T<:ThreeD}
#     r_atom = atoms.r[n] # == get_one_atom(atoms, n)

#     dot_n_r = cis(-dot(sensor_versor, r_atom))
#     E_scatt = dot_n_r * β[n]

#     ikr = im * k₀ * norm(sensor_position)
#     E_scattFinal = E_scatt * cis(k₀ * norm(sensor_position)) / ikr
#     return E_scattFinal
# end


# function get_scattered_intensity(problem::SimulationScalar, atoms_states, θ::Number; Gₙₘ=nothing)
#     if isnothing(Gₙₘ)
#         Gₙₘ = get_geometric_factor(problem.atoms, θ)
#     end

#     β = atoms_states
    
#     βₙₘ = transpose(β*β')
#     βₙₘ[diagind(βₙₘ)] .= abs2.(β)

#     intensity = real(sum(βₙₘ.*Gₙₘ))
# end




# ### --------------- MEAN FIELD ---------------
# """
#     get_scattered_intensity(problem::T, atoms_states, sensors::AbstractArray) where {T<:MeanFieldProblem}

# atoms_states: vcat(β, z)
# sensor: [ [sensor1], [sensor2], ... ,[sensorN]  ]

# **important**: this function does not handles single sensor.
# If needed, create a dummy position or duplicate the sensor - then, ignore the one that you don't need.
# sensor = [ [sensor1], [sensor_dummy]  ] or [ [sensor1], [sensor1]  ]

# """
# function get_scattered_intensity(problem::T, atoms_states, sensors::AbstractArray) where {T<:MeanFieldProblem}
#     N = problem.atoms.N
#     β = view(atoms_states, 1:N)
#     z = view(atoms_states, (N+1):2N)    
#     r = get_atoms_matrix(problem.atoms)
#     number_configurations = ((N^2)÷2 - N÷2)
    
#     βₙₘ = Array{ComplexF64}(undef, number_configurations)
#     cont_i = 0
#     cont_f = 0
#     for n=1:N-1
#         mm = (n+1):N

#         cont_i = cont_i + 1
#         cont_f = cont_f + length(mm)

#         βₙₘ[cont_i:cont_f] = conj(β[n])*β[mm]
    
#         cont_i = cont_f
#     end
#     vβₙₘ = view(βₙₘ, :)
    
#     cont_i = 0
#     cont_f = 0
#     rₙₘ = Array{ComplexF64}(undef, number_configurations, 3)
#     for n=1:N-1
#         m = (n+1):N

#         cont_i = cont_i + 1
#         cont_f = cont_f + length(m)

#         rₙₘ[cont_i:cont_f,:] = repeat( r[n,:]', length(m), 1) - @view r[m,:] # 

#         cont_i = cont_f
#     end
#     vrₙₘ = view(transpose(rₙₘ),:, :) # transpose necessary only on CPU mode, but I remove it with CUDA
#     #=
#         We don't need to compute each sensor in parallel, because the Folds.mapreduce
#         is already doing an excelent job with multi-threading.
#     =#
#     n_sensors = CoupledDipole.get_number_sensors(sensors)
#     intensities = Float64[]
#     for i=1:n_sensors
#         intensity = ComplexF64(0)
#         n_hat = sensors[i]./norm(sensors[i])
#         intensity = Folds.mapreduce(+, 1:number_configurations) do k
#             (
#                 begin 
#                 dot_n_r = n_hat[1]*vrₙₘ[1, k] + n_hat[2]*vrₙₘ[2, k] + n_hat[3]*vrₙₘ[3, k];
#                 vβₙₘ[k]*cis( k₀*dot_n_r) 
#                 end 
#             )
#         end
#         intensity +=  sum( (1 .+ z)./4 )
#         push!(intensities, 2real(intensity)   )
#     end

#     return intensities
# end

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
