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

function get_number_modes(problem::LinearOptics{Scalar})
    return problem.atoms.N
end
"""
    get_spectrum(problem; forceComputation=false)
return ωₙ, Γₙ
"""
function get_spectrum(problem; forceComputation=false)
    @debug "start: get spectrum"
    ωₙ, Γₙ = zeros(problem.atoms.N), zeros(problem.atoms.N)

    if is_spectrum_NOT_available(problem) || forceComputation
        H = get_interaction_matrix(problem)
        spectrum = eigen(H)

        problem.spectrum[:λ] = spectrum.values
        problem.spectrum[:ψ] = spectrum.vectors
        
        ωₙ, Γₙ = imag.(spectrum.values), -real.(spectrum.values)
        if any(Γₙ .< 0 )
            @warn "some Γₙ were negatives and ignored without further investigation"
            Γₙ = abs.(Γₙ)
        end
    else
        ωₙ, Γₙ = imag.(problem.spectrum[:λ]), -real.(problem.spectrum[:λ])
    end

    @debug "end  : get spectrum"
    return ωₙ, Γₙ
end

function is_spectrum_NOT_available(problem)
    if iszero(problem.spectrum[:λ] ) || iszero(problem.spectrum[:ψ] )
        return true
    else
        return false
    end
end

function get_ψ(problem::LinearOptics, n::Integer)
    ψ = view(problem.spectrum[:ψ],:,n)
    return ψ
end
function get_ψ²(problem::LinearOptics, n::Integer)
    ψ = problem.spectrum[:ψ][:,n]
    return abs2.(ψ)
end



"""
    get_localization_length(LinearOptics; forceComputation=false)

Main function to obtain the `Localization Length`, ξ, and its `Coefficient of Determination`, R1 (R¹).
"""
function get_localization_length(problem::LinearOptics; forceComputation=false)
    if is_localization_NOT_available(problem) && forceComputation==false
        return problem.data[:ξ], problem.data[:R1]
    end
    N = problem.atoms.N
    ξₙ = zeros(N)
    R¹ₙ = zeros(N) 

    # create spectrum if neeeded
    get_spectrum(problem)

    for n in 1:N # faster without Threads.@threads
        DCM, ψ² = get_spatial_profile_single_mode(problem, n)
        ξₙ[n], R¹ₙ[n] = get_single_ξ_and_R1(DCM, ψ²)
    end
    problem.data[:ξ], problem.data[:R1] = ξₙ, R¹ₙ
    return ξₙ, R¹ₙ
end

function is_localization_NOT_available(problem)
    if haskey(problem.data, :ξ) && haskey(problem.data, :R1)
        return true
    else
        return false
    end
    
end

function get_spatial_profile_single_mode(problem, mode_index::Integer)
    r = problem.atoms.r
    try
        ψ²ₙ = get_ψ²(problem, mode_index)
        r_cm = get_coordinates_of_center_of_mass(r, ψ²ₙ)
        
        # Specific for the problem (ones needs to define this function)
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
    N = length(Ψ²_mode)
    dimensions = size(r, 1)
    r_cm = zeros(dimensions)
    j = 1
    for c in eachcol(r)
        r_cm .+= c.*Ψ²_mode[j]
        j += 1
    end
    return r_cm./ sum(Ψ²_mode)
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
    classify_modes(problem)
returns a tuple `(loc, sub, super)` with indices
"""
function classify_modes(problem; forceComputation=false) 
    ωₙ, Γₙ = get_spectrum(problem;             forceComputation=forceComputation)
    ξₙ, R¹ₙ = get_localization_length(problem; forceComputation=forceComputation)

    localized_modes = findall((Γₙ .< Γ) .* (R¹ₙ .≥ R1_threshold))
    sub_radiant_modes = findall((Γₙ .< Γ) .* (R¹ₙ .< R1_threshold))
    super_radiant_modes = findall(Γₙ .> Γ)
    return (loc=localized_modes, sub=sub_radiant_modes, super=super_radiant_modes)
end

