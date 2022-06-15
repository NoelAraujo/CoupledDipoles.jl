function _OnePoint_Intensity(physic::Union{Scalar,Vectorial}, laser, atoms, sensor, β, scattering_func)
    E_L = laser_field(laser, sensor)
    E_scatt = scattering_func(atoms, β, sensor)
    return abs2(E_L + E_scatt)
end

@views function get_intensity_over_an_angle(problem::LinearOptics{Scalar}, atoms_states::AbstractVector, θ::Number)
    N = problem.atoms.N
    r = problem.atoms.r
    βₙ = view(atoms_states, 1:N)
    βₘ = conj.(βₙ)

    k₀sinθ = abs(k₀*sin(θ))
    k₀cosθ = k₀*cos(θ)
    number_configurations = ((N^2) ÷ 2 - N ÷ 2)

    xₙₘ = Array{Float64}(undef, number_configurations)
    yₙₘ, zₙₘ = similar(xₙₘ), similar(xₙₘ)
    count = 1
    for m in 1:N    
        rm = r[:,m]
        for n in (m+1):N
            xₙₘ[count] = r[1,n] - rm[1]
            yₙₘ[count] = r[2,n] - rm[2]
            zₙₘ[count] = r[3,n] - rm[3]
            count += 1
        end
    end

    βₙₘ = Array{ComplexF64}(undef, number_configurations)
    count = 1
    for m in 1:N
        b_m = βₘ[m]
        for n in (m+1):N
            b_n = βₙ[n]
            βₙₘ[count] = b_n*b_m
            count += 1
        end
    end
    
    total_intensity = ThreadsX.mapreduce(+, 1:number_configurations) do ii
        βₙₘ[ii]*exp(-im*zₙₘ[ii]*k₀cosθ)*besselj(0,k₀sinθ*sqrt(xₙₘ[ii]^2+yₙₘ[ii]^2))
    end
    total_intensity += sum(abs2, βₙ)/2
    return 2real(total_intensity)

end

@views function get_intensity_over_an_angle(problem::LinearOptics{Scalar}, atoms_states::AbstractVector, θ_range::AbstractVector)
    N = problem.atoms.N
    r = problem.atoms.r
    xₙ, yₙ, zₙ = r[1,:], r[2,:], r[3,:]

    βₙ = view(atoms_states, 1:N)
    βₘ = conj.(βₙ)
    
    number_configurations = ((N^2) ÷ 2 - N ÷ 2)

    xₙₘ = Array{Float64}(undef, number_configurations)
    yₙₘ, zₙₘ = similar(xₙₘ), similar(xₙₘ)
    count = 1
    for m in 1:N    
        for n in (m+1):N
            xₙₘ[count] = xₙ[n] - xₙ[m]
            yₙₘ[count] = yₙ[n] - yₙ[m]
            zₙₘ[count] = zₙ[n] - zₙ[m]
            count += 1
        end
    end

    βₙₘ = Array{ComplexF64}(undef, number_configurations)
    count = 1
    for m in 1:N
        b_m = βₘ[m]
        for n in (m+1):N
            b_n = βₙ[n]
            βₙₘ[count] = b_n*b_m
            count += 1
        end
    end

    intensities = zeros(length(θ_range))
    for (idx_θ, θ) in enumerate(θ_range)
        k₀sinθ = abs(k₀*sin(θ))
        k₀cosθ = k₀*cos(θ)

        _intensity = ThreadsX.mapreduce(+, 1:number_configurations) do ii
            βₙₘ[ii]*exp(-im*zₙₘ[ii]*k₀cosθ)*besselj(0,k₀sinθ*sqrt(xₙₘ[ii]^2+yₙₘ[ii]^2))
        end
        _intensity += sum(abs2, βₙ)/2

        intensities[idx_θ] = 2real(_intensity)
    end
    
    return intensities
end

