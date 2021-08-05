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

function get_intensity_over_an_angle(problem::LinearOptics{Scalar}, atoms_states::AbstractVector, θ::Number)  
    N = problem.atoms.N
    β = view(atoms_states, 1:N)

    ϕ_range = range(0,2π,length=30)
    vr = view(problem.atoms.r,:,:)
    
    complex_intensity = zeros(ComplexF64, N)
    total_intensity = 0.0
    
    for k= 1:length(ϕ_range)
        ϕ = ϕ_range[k]
        Threads.@threads for j = 1:N
            complex_intensity[j] = cis(-k₀*(vr[1,j]*sin(θ)*cos(ϕ) + vr[2,j]*sin(θ)*sin(ϕ) + vr[3,j]*cos(θ)) )*β[j] 
        end                 
        total_intensity += abs2( sum(complex_intensity) )
    end
    return total_intensity/length(ϕ_range)
end

function get_intensity_over_an_angle(problem::LinearOptics{Scalar}, atoms_states::Matrix, θ::Number)  
    timeSteps = size(atoms_states, 2)
    intensities = zeros(timeSteps)
    for i = 1:timeSteps
        oneState = view(atoms_states,:, i)
        intensities[i] = get_intensity_over_an_angle(problem, oneState, θ)
    end

    return intensities
end


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



function get_intensity_over_an_angle(problem::NonLinearOptics{MeanField}, atoms_states::Vector, θ::Float64)
    @debug "start : get intensity over an angle - NonLinearOptics{MeanField}"
    
    if is_integration_const_term_available(problem)
        Gₙₘ = problem.data[:Gₙₘ]
    else
        Gₙₘ = _get_meanfield_constant_term(problem.atoms, θ)
        problem.data[:Gₙₘ] = Gₙₘ
    end

    N = problem.atoms.N
    β = view(atoms_states, 1:N)
    z = view(atoms_states, (N+1):2N)

    βₙₘ = transpose(β*β') # I have to do "transpose" and NOT "adjoint = complex+tranpose"
    βₙₘ[diagind(βₙₘ)] .= (1 .+ z)./2
    # IMPORTANT: one does ELEMENT WISE multiplication
    # and one can use lesse memory with inplace multiplication
    βₙₘ .*=Gₙₘ 
    intensity = real(sum(βₙₘ)) 

    @debug "end  : get intensity over an angle - NonLinearOptics{MeanField}"
    return intensity
end

function is_integration_const_term_available(problem)
    if haskey(problem.data, :Gₙₘ)
        return true
    else
        return false
    end
end

function _get_meanfield_constant_term(atoms, Θ)
    N, r = atoms.N, atoms.r

    xₙₘ, yₙₘ, zₙₘ = get_xyz_distances(r)
    k₀sinΘ = k₀*sin(Θ)
    cos_Θ = cos(Θ)

    Gₙₘ_shared = SharedArray{ComplexF64,2}(N, N)    
    @sync for n=1:N
        Threads.@spawn for m=1:N # we had to compute all terms, and not the upper part
            @inbounds Gₙₘ_shared[n,m] = _constant_term_core_computation(xₙₘ,yₙₘ,zₙₘ, n,m, cos_Θ, k₀sinΘ)
        end
    end
    Gₙₘ = Array(Gₙₘ_shared) 
    #=
        IF I want to compute only the upper part,
        I have to multiply all terms by π/2:   Gₙₘ = (π/2)Array(Gₙₘ_shared)

        I decided to don't make this, to don't appear with factors
        not mentioned on theory.
    =#
    #= 
        Before returning, we HAVE to do some memory cleaning,      
        EVEN to get more performance. 

        Without this cleaning, the garbage collector gets lost 
        outside this function when many simulation occurs at the same time.
    =#
    Gₙₘ_shared = xₙₘ =  yₙₘ = zₙₘ = 1; GC.gc()  # DO NOT DELETE
    return Gₙₘ
end
function _constant_term_core_computation(xₙₘ,yₙₘ,zₙₘ, n,m, cos_Θ, k₀sinΘ)
    a = zero(ComplexF64)
    a = cis(k₀*zₙₘ[n,m]*cos_Θ)*besselj(0, k₀sinΘ*sqrt(xₙₘ[n,m]^2 +yₙₘ[n,m]^2) )
    return a
end
function get_xyz_distances(r) # @memoize  --> creating warning, let's ignore it right now
    dimensions = size(r, 1)
    N = size(r, 2)
    
    xₙₘ = SharedArray{Float64,2}(N, N)
    yₙₘ = SharedArray{Float64,2}(N, N)
    zₙₘ = SharedArray{Float64,2}(N, N)

    r_shared = SharedArray{Float64,2}(dimensions, N)
    r_shared .= r
    @sync for n=1:N
        r_n = view(r_shared,:,n)
        Threads.@spawn for m=1:N
            xₙₘ[n,m] = r_n[1] - r_shared[1,m]
            yₙₘ[n,m] = r_n[2] - r_shared[2,m]
            zₙₘ[n,m] = r_n[3] - r_shared[3,m]
        end
    end
    
    return xₙₘ, yₙₘ, zₙₘ
end
