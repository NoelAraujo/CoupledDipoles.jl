"""
    IPR = Inverse Participation Ratio
"""
function get_IPRs(problem)
    n_modes = get_number_modes(problem)
    IPRs = zeros(n_modes)
    Threads.@threads for n in 1:n_modes
        Ψ² = get_ψ²(problem, n)
        Ψ⁴ = Ψ² .^ 2
        IPRs[n] = sum(Ψ⁴) / sum(Ψ²)
    end
    return IPRs
end
function get_PRs(problem)
    IPRs = get_IPRs(problem)
    return 1.0 ./ IPRs
end

function get_number_modes(problem::T) where {T<:Scalar}
    return problem.atoms.N
end
"""
    get_spectrum(problem; forceComputation=false)
returns tuple (λ, ψ)
"""
function get_spectrum(problem; forceComputation=false)
    if is_spectrum_NOT_available(problem) || forceComputation
        H = get_interaction_matrix(problem)
        spectrum = eigen(H)

        problem.λ = spectrum.values
        problem.ψ = spectrum.vectors
    end
    return (λ=problem.λ, ψ=problem.ψ)
end

function is_spectrum_NOT_available(problem)
    if isnothing(problem.λ) || isnothing(problem.ψ)
        return true
    else
        return false
    end
end

"""
    get_all_ξ_and_R1(problem; probability_threshold=0.999)

Main function to obtain the `Localization Length`, ξ, and its `Coefficient of Determination`, R1 (R¹).
"""
function get_all_ξ_and_R1(problem; probability_threshold=0.999)
    if !isnothing(problem.ξ) && !isnothing(problem.R1)
        return problem.ξ, problem.R1
    end
    N = problem.atoms.N
    ξ = zeros(N)
    R1 = zeros(N)

    # create spectrum if neeeded
    get_spectrum(problem)

    Threads.@threads for n in 1:N
        DCM, ψ² = get_spatial_profile_single_mode(problem, n)
        ξ[n], R1[n] = get_single_ξ_and_R1(DCM, ψ²)
    end
    problem.ξ = ξ
    problem.R1 = R1
    return ξ, R1
end

function get_spatial_profile_single_mode(problem, mode_index::Integer)
    r = problem.atoms.r
    try
        ψ²ₙ = get_ψ²(problem, mode_index)
        r_cm = get_coordinates_of_center_of_mass(r, ψ²ₙ)

        # Specific for the problem (ones needs to define this function)
        # DCM = get_Distances_from_r_to_CM(convert_StaticArray_to_matrix(r), r_cm)
        DCM = get_Distances_from_r_to_CM(r, r_cm)

        sort_spatial_profile!(DCM, ψ²ₙ)
        return DCM, ψ²ₙ
    catch
        @error(
            "Probably you tried to access a spectrum not created (run `get_spectrum(problem)`). Or `mode_index` is not valid"
        )
    end
end
"""
    get_coordinates_of_center_of_mass(r, Ψ²_mode)
"""
function get_coordinates_of_center_of_mass(r, Ψ²_mode)
    ## Equivalent code, but slower
    # N = length(Ψ²_mode)
    # dimensions = size(r[1], 1) # index 1, because at least one element exist
    # r_cm = zeros(dimensions)
    # for n = 1:N
    #     r_cm .+= r[n].*Ψ²_mode[n]
    # end
    # return r_cm./ sum(Ψ²_mode)

    r_cm = sum(r .* Ψ²_mode)
    return r_cm ./ sum(Ψ²_mode)
end

"""
    get_Distances_from_r_to_CM(r, r_CM)

YOU need to specify for how to get the distance from points.
By default, I compute the Canonical Euclidian Distance
"""
get_Distances_from_r_to_CM(r, r_CM) = get_Distance_A_to_b(r, r_CM)

function sort_spatial_profile!(DCM, ψ²ₙ)
    idx = sortperm(DCM)
    DCM[:] = DCM[idx]
    ψ²ₙ[:] = ψ²ₙ[idx]
    return nothing
end

### --------------- Classification r---------------
"""
    classify_modes(problem::T) where {T<:Scalar}
returns a tuple `(loc, sub, super)` with indices
"""
function classify_modes(problem::T) where {T<:Scalar}
    ωₙ, Γₙ = get_energy_shift_and_linewith(problem)
    ξ, R1 = get_all_ξ_and_R1(problem)

    localized_modes = findall((Γₙ .< Γ) .* (R1 .≥ R1_threshold))
    sub_radiant_modes = findall((Γₙ .< Γ) .* (R1 .< R1_threshold))
    super_radiant_modes = findall(Γₙ .> Γ)
    return (loc=localized_modes, sub=sub_radiant_modes, super=super_radiant_modes)
end

