## necessary for laser intensity
function _get_intensity(problem::NonLinearOptics{MeanField}, field::AbstractMatrix)
    I_meanField = abs2.(field)
    return vec(real(I_meanField))
end


function _get_intensity_near_field(problem::NonLinearOptics{MeanField}, field, atomic_states, R, r; inelasticPart=true)
    N = problem.atoms.N
    σ⁻ = atomic_states[1:N]
    σᶻ = atomic_states[N+1:2N]

    I_meanField = abs2.(field)
    if inelasticPart
        I_meanField = I_meanField + real((Γ/(2k₀))^2*sum( (- abs2(σ⁻[j]) + (1 + σᶻ[j])/2)/sum(abs2, R-r[:,j]) for j=1:N))
    end
    return I_meanField
end


## Matrix case happens for single sensor field
function _get_intensity_far_field(problem::NonLinearOptics{MeanField}, field::AbstractMatrix, atomic_states::AbstractVector, R::Number; inelasticPart=true)
    σ⁻ = atomic_states[1]
    σᶻ = atomic_states[2]

    I_meanField = abs2.(field)
    if  inelasticPart
        I_meanField = I_meanField + (Γ/(2k₀*R))^2*( - abs2(σ⁻) + (1 + σᶻ)/2)
    end
    return vec(real(I_meanField))
end
function _get_intensity_far_field(problem::NonLinearOptics{MeanField}, field, atomic_states, R; inelasticPart=true)
    N = problem.atoms.N
    σ⁻ = atomic_states[1:N]
    σᶻ = atomic_states[N+1:2N]

    I_meanField = abs2.(field)
    if  inelasticPart
        population_correction = real((Γ/(2k₀*R))^2*sum( - abs2(σ⁻[j]) + (1 + σᶻ[j])/2 for j=1:N))
        I_meanField = I_meanField + population_correction
    end

    return I_meanField
end