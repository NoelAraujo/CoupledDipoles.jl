function _OnePoint_Intensity(physic::Union{Scalar,Vectorial}, laser, atoms, sensor, β, scattering_func)
    E_L = laser_field(laser, sensor)
    E_scatt = scattering_func(atoms, β, sensor)
    return abs2(E_L + E_scatt)
end
using QuadGK
@views function get_intensity_over_an_angle(problem::LinearOptics{Scalar}, atoms_states::AbstractVector, θ::Number)
    N = problem.atoms.N
    r = problem.atoms.r
    β = view(atoms_states, 1:N)
    βp = conj.(β)

    x, y, z = r[1,:], r[2,:], r[3,:]
    intensity = zero(ComplexF64)
    for j=1:N
        for jp=1:N
            xjj = x[j] - x[jp]
            yjj = y[j] - y[jp]
            zjj = z[j] - z[jp]

            intensity += myKernel2(β[j], βp[jp], xjj, yjj, zjj, θ)
        end
    end
    return real(intensity)
 
    return real(total_intensity)
end
function myKernel2(beta_j, beta_jp, xjj, yjj, zjj, θ)
    k₀sinθ = abs(k₀*sin(θ))
    k₀cosθ = k₀*cos(θ)

    v = beta_j*beta_jp*exp(-im*zjj*k₀cosθ)*besselj(0,k₀sinθ*sqrt(xjj^2+yjj^2))
    return v
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
