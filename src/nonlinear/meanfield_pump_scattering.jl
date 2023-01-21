#=
            LASER FIELD
=#
function laser_field(problem::NonLinearOptics{MeanField}, sensor::AbstractVector)
    Ω₀ = raby_frequency(problem.laser)
    return Matrix(transpose([LASER_FACTOR*Ω₀*_scalar_laser_field(problem.laser, sensor)]))
end
function laser_field(problem::NonLinearOptics{MeanField}, sensors::AbstractMatrix)
    Ω₀ = raby_frequency(problem.laser)

    _laser_electric_fields = ThreadsX.map(eachcol(sensors)) do sensor
            LASER_FACTOR*Ω₀*_scalar_laser_field(problem.laser, sensor)
    end
    laser_electric_fields::Matrix{ComplexF64} = hcat(_laser_electric_fields...)
    return laser_electric_fields
end

## necessary for laser intensity
function _get_intensity(problem::NonLinearOptics{MeanField}, field::AbstractMatrix)
    I_meanField = abs2.(field)
    return vec(real(I_meanField))
end


function _get_intensity_near_field(problem::NonLinearOptics{MeanField}, field, atomic_states, R, r)
    N = problem.atoms.N
    σ⁻ = atomic_states[1:N]
    σᶻ = atomic_states[N+1:2N]
    proposed_intensity = abs2.(field) + real((Γ/(2k₀))^2*sum( (- abs2(σ⁻[j]) + (1 + σᶻ[j])/2)/sum(abs2, R-r[:,j]) for j=1:N))

     ## --> This topic needs further reasoning <--
    # for low saturation the values of  σᶻ maybe smaller than σ⁻, leading to negative values of intensity
    if proposed_intensity > 0
        I_meanField = proposed_intensity
    else
        I_meanField = abs2.(field)
    end
    return I_meanField
end


## Matrix case happens for single sensor field
function _get_intensity_far_field(problem::NonLinearOptics{MeanField}, field::AbstractMatrix, atomic_states::AbstractVector, R::Number)
    σ⁻ = atomic_states[1]
    σᶻ = atomic_states[2]
    I_meanField = abs2.(field) + (Γ/(2k₀*R))^2*( - abs2(σ⁻) + (1 + σᶻ)/2)
    return vec(real(I_meanField))
end
function _get_intensity_far_field(problem::NonLinearOptics{MeanField}, field, atomic_states, R)
    N = problem.atoms.N
    σ⁻ = atomic_states[1:N]
    σᶻ = atomic_states[N+1:2N]
    population_correction = real((Γ/(2k₀*R))^2*sum( - abs2(σ⁻[j]) + (1 + σᶻ[j])/2 for j=1:N))
    proposed_intensity = abs2.(field) + population_correction

    ## --> This topic needs further reasoning <--
    # for low saturation the values of  σᶻ maybe smaller than σ⁻, leading to negative values of intensity
    if proposed_intensity > 0
        I_meanField = proposed_intensity
    else
        I_meanField = abs2.(field)
    end

    return I_meanField
end







#=
            SCATTERED FIELD: :near_field
=#
function scattering_near_field(problem::NonLinearOptics{MeanField}, β, sensor)
    _scalar_scattering_near_field(problem.atoms, β, sensor)
end


function _meanfield_scattering_near_field(atoms::Atom{T}, β, sensor) where T <: TwoD
    # TO DO
    return nothing
end
function _meanfield_scattering_near_field(atoms::Atom{T}, β, sensor) where T <: ThreeD
    σ⁻ = β[1:N]
    E_scatt = _scalar_scattering_near_field(atoms, σ⁻, sensor)
    return E_scatt
end

#=
            SCATTERED FIELD: :far_field
=#
function scattering_far_field(problem::NonLinearOptics{MeanField}, β, sensor)
    _scalar_scattering_near_field(problem.atoms, β, sensor)
end

function _meanfield_scattering_far_field(atoms::Atom{T}, β, sensor) where T <: TwoD
    # TO DO
    return nothing
end
function _meanfield_scattering_far_field(atoms::Atom{T}, β, sensor) where T <: ThreeD
    σ⁻ = β[1:N]
    E_scatt = _scalar_scattering_far_field(atoms, σ⁻, sensor)
    return E_scatt
end
