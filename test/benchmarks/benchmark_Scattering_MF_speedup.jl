using CoupledDipoles, Revise, LinearAlgebra, Folds
using Random


#=
This file contain steps to compute improve the computation of
    Light scatterd on the far field approximation:

I(θ, φ) = ∑ₙ∑ₘ conj(βₙ)βₘ exp( im*k₀*n⋅(rₙ - rₘ) ) + ∑ₙ (1 + zₙ)/2

∑ₙ∑ₘ = n=1:N; m=1:N; n≠m

βⱼ = atom expected value of σ⁻ (j ∈ [n,m])
zₙ = atom expected value of σᶻ
k₀ = unit of length of the system (k₀ = 1)
n = versor of the sensor ( n = r_sensor./norm(r_sensor) )
rⱼ = atom position (j ∈ [n,m])

=#

"""
    Following the definition of I(θ, φ) without any optimization
"""
function scattering_version_naive(β,z, n_hat, r; k₀=1)
    N = length(β)
    intensity = ComplexF64(0)

    for n=1:N
        for m=1:N
            if n≠m
                intensity += conj(β[n])*β[m]*exp( im*k₀*dot(n_hat, r[:, n]-r[:, m])   )
            end
        end
    end
    for n=1:N
        intensity += ( 1 + z[n] )/2
    end
    return intensity
end


"""
   - I compute `conj(β[n])*β[m]` once
   - avoid accessing position rₙ inside loop m (r_n = @view(r[:, n]))
   - put @view to access all collumns
"""
function scattering_optimization_1(β, z, n_hat, r; k₀=1)
    N = length(β)
    
    βₙₘ = zeros(ComplexF64, N, N)
    for n=1:N
        for m=1:N
            if n≠m
                βₙₘ[n,m] = conj(β[n])*β[m]
            end
        end
    end

    intensity = ComplexF64(0)
    for n=1:N
        r_n = @view(r[:, n])
        for m=1:N
            if n≠m
                intensity += βₙₘ[n,m]*exp( im*k₀*dot(n_hat, r_n - @view(r[:,m])  )   )
            end
        end
    end
    for n=1:N
        intensity += (1+z[n])/2
    end
    return intensity
end




"""
   - `zeros(N,N)` --> `Array{ComplexF64}(undef, N,N)`
   - expand the dot product
   - change `exp(im*x)` by `cis(x)`
   - removed the loop for `∑ₙ (1 + zₙ)/2`, and used a `sum` function
"""
function scattering_optimization_2(β, z, n_hat, r; k₀=1)
    N = length(β)
    
    βₙₘ = Array{ComplexF64}(undef, N,N)
    for n=1:N
        for m=1:N
            if n≠m
                βₙₘ[n,m] = conj(β[n])*β[m]
            end
        end
    end

    intensity = ComplexF64(0)
    r_nm = zeros(3)
    for n=1:N
        r_n = @view(r[:, n])
        for m=1:N
            if n≠m
                r_nm .= r_n - @view(r[:,m])
                dot_n_r = n_hat[1]*r_nm[1] + n_hat[2]*r_nm[2] + n_hat[3]*r_nm[3]
                intensity += βₙₘ[n,m]*cis( k₀*dot_n_r  )
            end
        end
    end
    
    intensity += sum( (1 .+ z)./2 )
    
    return intensity
end


"""
   - stored in a matriz the vectors of `rₙ - rₘ`
"""
function scattering_optimization_3(β, z, n_hat, r; k₀=1)
    N = length(β)
    

    βₙₘ = Array{ComplexF64}(undef, N,N)
    for n=1:N
        for m=1:N            
            βₙₘ[n,m] = conj(β[n])*β[m]
        end
    end
    

    rₙₘ = Array{ComplexF64}(undef, 3, N^2)
    cont = 1
    for n=1:N
        r_n = @view(r[:,n])
        for m=1:N
            rₙₘ[:,cont] = r_n - @view(r[:,m])
            cont += 1
        end
    end

    intensity = ComplexF64(0)
    cont = 1    
    for n=1:N
        for m=1:N
            if n≠m
                dot_n_r = n_hat[1]*rₙₘ[1, cont] + n_hat[2]*rₙₘ[2, cont] + n_hat[3]*rₙₘ[3, cont]
                intensity += βₙₘ[n,m]*cis( k₀*dot_n_r  )
            end
            cont += 1
        end
    end
    
    intensity += sum( (1 .+ z)./2 )
    
    return intensity
end


"""
    - using complex number relation: 2real(z) = z + conj(z)
        - Then I return: 2real(intensity)
    - I don't need to compute `conj(z)` terms, which are the lower diagonal
        of βₙₘ and rₙₘ. 
    - Also, to avoid the `if n≠m`, I begin with upper diagonal  `m = (n+1):N`
    for n=1:N
        for m = (n+1):N  # <--- here is the change
            ...
        end
    end
    - I cannot forget to divide the diagonal by 2:
        `sum( (1 .+ z)./2 )/2` --> `sum( (1 .+ z)./4 )`
    - Isolate the main computation inside the function `aux4`
"""
function scattering_optimization_4(β, z, n_hat, r; k₀=1)
    N = length(β)

    number_configurations = ((N^2)÷2 - N÷2)

    βₙₘ = Array{ComplexF64}(undef, number_configurations)
    cont = 1
    for n=1:N
        for m=(n+1):N            
            βₙₘ[cont] = conj(β[n])*β[m]
            cont += 1
        end
    end    

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
    

    intensity = ComplexF64(0)
    cont = 1
    for n=1:N
        for m=(n+1):N
            intensity += aux4(n_hat, rₙₘ, βₙₘ, cont)
            cont += 1
        end
    end
    
    intensity +=  sum( (1 .+ z)./4 )
    return 2real(intensity)
end
function aux4(n_hat, rₙₘ, βₙₘ, cont; k₀=1)
    dot_n_r = n_hat[1]*rₙₘ[1, cont] + n_hat[2]*rₙₘ[2, cont] + n_hat[3]*rₙₘ[3, cont]
    return βₙₘ[cont]*cis( k₀*dot_n_r  )
end

"""
    - create views of βₙₘ and rₙₘ. 
    - do not call for external function to compute the intensity
"""
function scattering_optimization_5(β, z, n_hat, r; k₀=1)
    N = length(β)

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

    intensity = ComplexF64(0)
    for cont = 1:number_configurations
        dot_n_r = n_hat[1]*vrₙₘ[1, cont] + n_hat[2]*vrₙₘ[2, cont] + n_hat[3]*vrₙₘ[3, cont]
        intensity += vβₙₘ[cont]*cis( k₀*dot_n_r  )
    end
    
    intensity +=  sum( (1 .+ z)./4 )
    return 2real(intensity)
end

"""
    - use multi-threading of `Folds.mapreduce`
    - I still need to use another @view when slicing vrₙₘ
"""
function scattering_optimization_6(β, z, n_hat, r; k₀=1)
    N = length(β)

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

    intensity = ComplexF64(0)
    intensity = Folds.mapreduce(+, 1:number_configurations) do k
        (  aux_opt6(n_hat, @view(vrₙₘ[:, k]), vβₙₘ[k]) )
    end
    intensity +=  sum( (1 .+ z)./4 )
    return 2real(intensity)
end
function aux_opt6(n_hat, vrₙₘ, vβₘₙ; k₀=1)
    dot_n_r = cis( k₀*(n_hat[1]*vrₙₘ[1]
                    + n_hat[2]*vrₙₘ[2]
                    + n_hat[3]*vrₙₘ[3]))
   return dot_n_r*vβₘₙ
end


"""
    - do not call for external functions inside mapreduce
    - use @inbounds
"""
function scattering_optimization_7(β, z, n_hat, r; k₀=1)
    N = length(β)

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
    

    intensity = ComplexF64(0)
    intensity = Folds.mapreduce(+, 1:number_configurations) do k
        (
            begin 
                @inbounds dot_n_r = n_hat[1]*vrₙₘ[1, k] + n_hat[2]*vrₙₘ[2, k] + n_hat[3]*vrₙₘ[3, k];
                @inbounds vβₙₘ[k]*cis( k₀*dot_n_r) 
            end 
        )
    end
    intensity +=  sum( (1 .+ z)./4 )
    return 2real(intensity)
end


# initialize
N = 400
kL  = 10
s = 1e-6
Δ = 0.0
sensors = get_sensors_ring(;num_pts=90, kR=1.5kL, θ=5π/12)

Random.seed!(1301)
atoms = Atom(Cube(), N, kL)
laser = Laser(Gaussian3D(kL/8), s, Δ)
simulation = NonLinearOptics(MeanField(), atoms, laser)

β = rand(ComplexF64, N)
z = 2conj(β).*β .- 1
r = atoms.r
n = sensors[:,1]./norm(sensors[:,1])

print("naive:");@time  I_naive = scattering_version_naive(β,z, n, r);
print("opt1: "); @time I_opt_1 = scattering_optimization_1(β,z, n, r); @assert I_naive ≈ I_opt_1
print("opt2: "); @time I_opt_2 = scattering_optimization_2(β,z, n, r); @assert I_naive ≈ I_opt_2
print("opt3: "); @time I_opt_3 = scattering_optimization_3(β,z, n, r); @assert I_naive ≈ I_opt_3
print("opt4: "); @time I_opt_4 = scattering_optimization_4(β,z, n, r); @assert I_naive ≈ I_opt_4
print("opt5: "); @time I_opt_5 = scattering_optimization_5(β,z, n, r); @assert I_naive ≈ I_opt_5
print("opt6: "); @time I_opt_6 = scattering_optimization_6(β,z, n, r); @assert I_naive ≈ I_opt_6
print("opt7: "); @time I_opt_7 = scattering_optimization_7(β,z, n, r); @assert I_naive ≈ I_opt_7


# benchmarks
N = 6000

Random.seed!(1301)
atoms = Atom(Cube(), N, kL)
laser = Laser(Gaussian3D(kL/8), s, Δ)
simulation = NonLinearOptics(MeanField(), atoms, laser)

β = rand(ComplexF64, N)
z = 2conj(β).*β .- 1
r = atoms.r

print("naive:");@time  I_naive = scattering_version_naive(β,z, n, r);
print("opt1: "); @time I_opt_1 = scattering_optimization_1(β,z, n, r); @assert I_naive ≈ I_opt_1
print("opt2: "); @time I_opt_2 = scattering_optimization_2(β,z, n, r); @assert I_naive ≈ I_opt_2
print("opt3: "); @time I_opt_3 = scattering_optimization_3(β,z, n, r); @assert I_naive ≈ I_opt_3
print("opt4: "); @time I_opt_4 = scattering_optimization_4(β,z, n, r); @assert I_naive ≈ I_opt_4
print("opt5: "); @time I_opt_5 = scattering_optimization_5(β,z, n, r); @assert I_naive ≈ I_opt_5
print("opt6: "); @time I_opt_6 = scattering_optimization_6(β,z, n, r); @assert I_naive ≈ I_opt_6
print("opt7: "); @time I_opt_7 = scattering_optimization_7(β,z, n, r); @assert I_naive ≈ I_opt_7



"""
    Code for production is ideal up to 90 sensors,
        because the function takes advantage of CPU frequency burst
"""
function scattering_current_best(β, z, r, sensors; k₀=1)
    N = length(β)

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
@profview scattering_current_best(β, z, atoms.r, sensors)

for i=1:10
    @time scattering_current_best(β, z, atoms.r, sensors)
end

