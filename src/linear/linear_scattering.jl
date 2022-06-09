function _OnePoint_Intensity(physic::Union{Scalar,Vectorial}, laser, atoms, sensor, β, scattering_func)
    E_L = laser_field(laser, sensor)
    E_scatt = scattering_func(atoms, β, sensor)
    return abs2(E_L + E_scatt)
end

function get_intensity_over_an_angle(problem::LinearOptics{Scalar}, atoms_states::AbstractVector, θ::Number)
    N = problem.atoms.N
    β = view(atoms_states, 1:N)

    ϕ_range = range(0, 2π; length=30)
    vr = view(problem.atoms.r, :, :)

    complex_intensity = zeros(ComplexF64, N)
    total_intensity = 0.0

    rx = vr[1, :]
    ry = vr[2, :]
    rz = vr[3, :]

    for ϕ ∈ ϕ_range
        complex_intensity .= cis.(-k₀ .* (rx .* sin(θ) .* cos(ϕ) + ry .* sin(θ) * sin(ϕ) + rz.*cos(θ))) .* β
        total_intensity += abs2(sum(complex_intensity))
    end
    return total_intensity / length(ϕ_range)
end

function get_intensity_over_an_angle_shared(problem::LinearOptics{Scalar}, β::AbstractVector, θ::Number, r_shared)
    N = problem.atoms.N

    ϕ_range = range(0, 2π; length=30)

    complex_intensity = zeros(ComplexF64, N)
    total_intensity = 0.0

    for k in 1:length(ϕ_range)
        ϕ = ϕ_range[k]
        for j in 1:N
            complex_intensity[j] = cis(-k₀ * (r_shared[1, j] * sin(θ) * cos(ϕ) + r_shared[2, j] * sin(θ) * sin(ϕ) + r_shared[3, j] * cos(θ))) * β[j]
        end
        total_intensity += abs2(sum(complex_intensity))
    end
    return total_intensity / length(ϕ_range)
end

function get_intensity_over_an_angle(problem::LinearOptics{Scalar}, atoms_states::Matrix, θ::Number)
    timeSteps = size(atoms_states, 2)
    intensities = zeros(timeSteps)

    # r_shared = SharedArray{Float64,2}(3, problem.atoms.N)
    r_shared = problem.atoms.r

    Threads.@threads for i in 1:timeSteps
        oneState = view(atoms_states, :, i)
        intensities[i] = get_intensity_over_an_angle_shared(problem, oneState, θ, r_shared)
    end

    return intensities
end
