using CoupledDipoles, Revise, LinearAlgebra
using Random, ThreadsX


#=
This file contain steps to compute improve the computation of
    Light scatterd on the far field approximation:

I(θ, φ) = ∑ₙ∑ₘ conj(βₙ)βₘ exp( im*k₀*n⋅(rₙ - rₘ) )

∑ₙ∑ₘ = n=1:N; m=1:N; n≠m

βⱼ = atom expected value of σ⁻ (j ∈ [n,m])
k₀ = unit of length of the system (k₀ = 1)
n = versor of the sensor ( n = r_sensor./norm(r_sensor) )
rⱼ = atom position (j ∈ [n,m])

=#

"""
    Following the definition of I(θ, φ) without any optimization
"""
function scattering_version_naive(β, n_hat, r; k₀=1)
    N = length(β)
    intensity = ComplexF64(0)

    for n=1:N
        for m=1:N
            if n≠m
                intensity += conj(β[n])*β[m]*exp( im*k₀*dot(n_hat, r[:, n]-r[:, m])   )
            end
        end
    end
    return intensity
end
function scattering_multithreading_version(β, n_hat, r; k₀=1)
    intensity = ThreadsX.sum( conj(β[n])*β[m]*exp( im*k₀*dot(n_hat, r[:,n]-r[:,m])) for n = 1:N, m = 1:N if n ≠ m )
    return intensity
end

"""
   - I compute `conj(β[n])*β[m]` once
   - avoid accessing position rₙ inside loop m (r_n = r[:, n])
   - put @views to access all collumns without producing SubArray
"""
@views function scattering_optimization_1(β, n_hat, r; k₀=1)
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
        r_n = r[:, n]
        for m=1:N
            if n≠m
                intensity += βₙₘ[n,m]*exp( im*k₀*dot(n_hat, r_n - r[:,m]  )   )
            end
        end
    end
    return intensity
end




"""
   - `zeros(N,N)` --> `Array{eltype(β)}(undef, N,N)`
   - expand the dot product
   - change `exp(im*x)` by `cis(x)`
"""
@views function scattering_optimization_2(β, n_hat, r; k₀=1)
    N = length(β)
    
    βₙₘ = Array{eltype(β)}(undef, N,N)
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
        r_n = r[:, n]
        for m=1:N
            if n≠m
                r_nm .= r_n - r[:,m]
                dot_n_r = n_hat[1]*r_nm[1] + n_hat[2]*r_nm[2] + n_hat[3]*r_nm[3]
                intensity += βₙₘ[n,m]*cis( k₀*dot_n_r  )
            end
        end
    end
        
    return intensity
end


"""
   - remove the `if n≠m` condition from `βₙₘ` and `rₙₘ`, and only retain in the intensity
   - stored in a matriz the vectors of `rₙ - rₘ`
"""
@views function scattering_optimization_3(β, n_hat, r; k₀=1)
    N = length(β)
    

    βₙₘ = Array{eltype(β)}(undef, N,N)
    for n=1:N
        for m=1:N            
            βₙₘ[n,m] = conj(β[n])*β[m]
        end
    end
    

    rₙₘ = Array{eltype(r)}(undef, 3, N^2)
    cont = 1
    for n=1:N
        r_n = r[:,n]
        for m=1:N
            rₙₘ[1,cont] = r_n[1] - r[1,m]
            rₙₘ[2,cont] = r_n[2] - r[2,m]
            rₙₘ[3,cont] = r_n[3] - r[3,m]
            cont += 1
        end
    end

    intensity = zero(eltype(β))
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
    
    return intensity
end


"""
    - using complex number relation: 2real(z) = z + conj(z)
        - Then I return: 2real(intensity)
        - Therefore, I don't need to compute `conj(z)` terms, which are the lower diagonal
            of βₙₘ and rₙₘ. 
    - Also, to avoid the `if n≠m`, I begin with upper diagonal  `m = (n+1):N`
    for n=1:N
        for m = (n+1):N  # <--- here is the change
            ...
        end
    end
"""
@views function scattering_optimization_4(β, n_hat, r; k₀=1)
    N = length(β)
    number_configurations = ((N^2)÷2 - N÷2)

    βₙₘ = Array{eltype(β)}(undef, number_configurations)
    cont = 1
    for n=1:N
        for m=(n+1):N            
            βₙₘ[cont] = conj(β[n])*β[m]
            cont += 1
        end
    end
    

    rₙₘ = Array{eltype(r)}(undef, 3, number_configurations)
    cont = 1
    for n=1:N
        r_n = r[:,n]
        for m=(n+1):N
            rₙₘ[1,cont] = r_n[1] - r[1,m]
            rₙₘ[2,cont] = r_n[2] - r[2,m]
            rₙₘ[3,cont] = r_n[3] - r[3,m]
            cont += 1
        end
    end

    intensity = zero(eltype(β))
    cont = 1    
    for n=1:N
        for m=(n+1):N
            dot_n_r = n_hat[1]*rₙₘ[1, cont] + n_hat[2]*rₙₘ[2, cont] + n_hat[3]*rₙₘ[3, cont]
            intensity += βₙₘ[cont]*cis( k₀*dot_n_r  )
            cont += 1
        end
    end
    
    return 2real(intensity)
end


"""
    - use multi-threading of `Folds.mapreduce`
    - add `@inbounds`
"""
@views function scattering_optimization_5(β, n_hat, r; k₀=1)
    N = length(β)
    number_configurations = ((N^2)÷2 - N÷2)

    βₙₘ = Array{eltype(β)}(undef, number_configurations)
    cont = 1
    for n=1:N
        for m=(n+1):N            
            βₙₘ[cont] = conj(β[n])*β[m]
            cont += 1
        end
    end
    

    rₙₘ = Array{eltype(r)}(undef, 3, number_configurations)
    cont = 1
    for n=1:N
        r_n = r[:,n]
        for m=(n+1):N
            rₙₘ[1,cont] = r_n[1] - r[1,m]
            rₙₘ[2,cont] = r_n[2] - r[2,m]
            rₙₘ[3,cont] = r_n[3] - r[3,m]
            cont += 1
        end
    end

    intensity = ThreadsX.mapreduce(+, 1:number_configurations) do cont
        (  
           begin 
                @inbounds dot_n_r = n_hat[1]*rₙₘ[1, cont] + n_hat[2]*rₙₘ[2, cont] + n_hat[3]*rₙₘ[3, cont]
                @inbounds βₙₘ[cont]*cis( k₀*dot_n_r  )
           end 
        )
    end
    
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

β = rand(ComplexF64, N)
r = atoms.r
n = sensors[:,1]./norm(sensors[:,1])

print("naive:") ;@time  I_naive = scattering_version_naive(β, n, r);
print("threX:"); @time I_threads =scattering_multithreading_version(β, n, r); @assert I_naive ≈ I_threads
print("opt1: "); @time I_opt_1 = scattering_optimization_1(β, n, r); @assert I_naive ≈ I_opt_1
print("opt2: "); @time I_opt_2 = scattering_optimization_2(β, n, r); @assert I_naive ≈ I_opt_2
print("opt3: "); @time I_opt_3 = scattering_optimization_3(β, n, r); @assert I_naive ≈ I_opt_3
print("opt4: "); @time I_opt_4 = scattering_optimization_4(β, n, r); @assert I_naive ≈ I_opt_4
print("opt5: "); @time I_opt_5 = scattering_optimization_5(β, n, r); @assert I_naive ≈ I_opt_5



# benchmarks
N = 6000

Random.seed!(1301)
atoms = Atom(Cube(), N, kL)

β = rand(ComplexF64, N)
r = atoms.r

print("naive:");@time  I_naive = scattering_version_naive(β, n, r);
print("threX:"); @time I_threads =scattering_multithreading_version(β, n, r); @assert I_naive ≈ I_threads
print("opt1: "); @time I_opt_1 = scattering_optimization_1(β, n, r); @assert I_naive ≈ I_opt_1
print("opt2: "); @time I_opt_2 = scattering_optimization_2(β, n, r); @assert I_naive ≈ I_opt_2
print("opt3: "); @time I_opt_3 = scattering_optimization_3(β, n, r); @assert I_naive ≈ I_opt_3
print("opt4: "); @time I_opt_4 = scattering_optimization_4(β, n, r); @assert I_naive ≈ I_opt_4
print("opt5: "); @time I_opt_5 = scattering_optimization_5(β, n, r); @assert I_naive ≈ I_opt_5
