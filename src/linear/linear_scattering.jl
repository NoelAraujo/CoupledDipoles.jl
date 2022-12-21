# function _OnePoint_Intensity(physic::Union{Scalar,Vectorial}, laser, atoms, sensor, β, scattering_func)
#     E_L = laser_field(laser, sensor)
#     E_scatt = scattering_func(atoms, β, sensor)
#     return abs2(E_L + E_scatt)
# end

# @views function get_intensity_over_an_angle(problem::LinearOptics{Scalar}, atoms_states::AbstractVector, θ::Number)
#     N = problem.atoms.N
#     r = problem.atoms.r
#     βₙ = view(atoms_states, 1:N)
#     βₘ = conj.(βₙ)

#     k₀sinθ = abs(k₀*sin(θ))
#     k₀cosθ = k₀*cos(θ)
#     number_configurations = ((N^2) ÷ 2 - N ÷ 2)

#     xₙₘ = Array{Float64}(undef, number_configurations)
#     yₙₘ, zₙₘ = similar(xₙₘ), similar(xₙₘ)
#     count = 1
#     for m in 1:N
#         rm = r[:,m]
#         for n in (m+1):N
#             xₙₘ[count] = r[1,n] - rm[1]
#             yₙₘ[count] = r[2,n] - rm[2]
#             zₙₘ[count] = r[3,n] - rm[3]
#             count += 1
#         end
#     end

#     βₙₘ = Array{ComplexF64}(undef, number_configurations)
#     count = 1
#     for m in 1:N
#         b_m = βₘ[m]
#         for n in (m+1):N
#             b_n = βₙ[n]
#             βₙₘ[count] = b_n*b_m
#             count += 1
#         end
#     end

#     total_intensity = ThreadsX.mapreduce(+, 1:number_configurations) do ii
#         βₙₘ[ii]*exp(-im*zₙₘ[ii]*k₀cosθ)*Bessels.besselj0(k₀sinθ*sqrt(xₙₘ[ii]^2+yₙₘ[ii]^2))
#     end
#     total_intensity += sum(abs2, βₙ)/2
#     return 2real(total_intensity)

# end

"""
    used for the same state and many angles (suitable for cbs)
"""
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
            βₙₘ[ii]*exp(-im*zₙₘ[ii]*k₀cosθ)*Bessels.besselj0(k₀sinθ*sqrt(xₙₘ[ii]^2+yₙₘ[ii]^2))
        end
        _intensity += sum(abs2, βₙ)/2

        intensities[idx_θ] = 2real(_intensity)
    end

    return intensities
end




"""
    used for the single angle and single single state (most probably user case)
"""
function get_intensity_over_an_angle(problem::LinearOptics{Scalar}, atoms_states::Vector, θ::Number; tol=exp10(-7.4))
    if problem.atoms.N < 1000
        xⱼₘ2, yⱼₘ2, zⱼₘ = _rⱼₘ_distances(problem)
        return _intensity_angle_exact_parallel(problem, atoms_states, θ, xⱼₘ2, yⱼₘ2, zⱼₘ)
    else
        return _intensity_angle_approx_quadradure(problem, atoms_states, θ; tol=tol)
    end
end
"""
    used for the same angle and different states (for example, the output of `time_evolution`)
"""
function get_intensity_over_an_angle(problem::LinearOptics{Scalar}, atoms_states::Vector{Vector{ComplexF64}}, θ::Number; tol=exp10(-7.4), exact_solution=false)

    if exact_solution # O((N^2 - N)/2)
        xⱼₘ2, yⱼₘ2, zⱼₘ = _rⱼₘ_distances(problem)
        return map(atoms_states) do β
            _intensity_angle_exact_parallel(problem, β, θ, xⱼₘ2, yⱼₘ2, zⱼₘ)
        end
    end

    return ThreadsX.map(atoms_states) do β # O(N)
        _intensity_angle_approx_quadradure(problem, β, θ; tol=tol)
    end
end
@views function _rⱼₘ_distances(problem)
    N = problem.atoms.N
    r = problem.atoms.r

    xⱼ, yⱼ, zⱼ = r[1,:], r[2,:], r[3,:]

    number_configurations = ((N^2) ÷ 2 - N ÷ 2)

    xⱼₘ2 = Array{Float64}(undef, number_configurations)
    yⱼₘ2, zⱼₘ = similar(xⱼₘ2), similar(xⱼₘ2)
    count = 1
    for m in 1:N
        for n in (m+1):N
            xⱼₘ2[count] = (xⱼ[n] - xⱼ[m])^2
            yⱼₘ2[count] = (yⱼ[n] - yⱼ[m])^2
            zⱼₘ[count] = zⱼ[n] - zⱼ[m]
            count += 1
        end
    end

    return xⱼₘ2, yⱼₘ2, zⱼₘ
end



@views function _intensity_angle_exact_parallel(problem::LinearOptics{Scalar}, atoms_states::AbstractVector, θ::Number, xⱼₘ2, yⱼₘ2, zⱼₘ)
    N = problem.atoms.N
    βₙ = view(atoms_states, 1:N)
    βₘ = conj.(βₙ)

    number_configurations = ((N^2) ÷ 2 - N ÷ 2)

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

    k₀sinθ = abs(k₀*sin(θ))
    k₀cosθ = k₀*cos(θ)
    total_intensity = ThreadsX.mapreduce(+, 1:number_configurations) do ii
        βₙₘ[ii]*cis(-zⱼₘ[ii]*k₀cosθ)*Bessels.besselj0(k₀sinθ*sqrt(xⱼₘ2[ii] + yⱼₘ2[ii]))
    end
    total_intensity += sum(abs2, βₙ)/2

    k₀R = k₀*how_far_is_farField(problem)

    return (2π)*Γ^2/(4*k₀R^2)*2real(total_intensity)

end





@views function _intensity_angle_approx_quadradure(problem::LinearOptics{Scalar}, atoms_states::AbstractVector, θ::Number; tol=exp10(-7.4))
    R = how_far_is_farField(problem)

    sensor = Array{Float64}(undef,3)
    (intensity, _e) = hcubature([0.0], [2π], rtol = tol, atol=tol) do ϕ
        sensor[1] = R*sin(θ)*cos(ϕ[1])
        sensor[2] = R*sin(θ)*sin(ϕ[1])
        sensor[3] = R*cos(θ)
        scattered_intensity(problem, atoms_states, sensor; regime=:far_field)
    end
    return  intensity[1]
end