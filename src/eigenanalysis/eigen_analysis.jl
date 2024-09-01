"""
    get_IPRs(problem)

Inverse Participation Ratio is computed with `∑Ψ⁴ / ( ∑Ψ² )²` for each mode
"""
function get_IPRs(problem)
    n_modes = get_number_modes(problem)
    IPRs = zeros(n_modes)
    Threads.@threads for n in 1:n_modes
        Ψ² = get_ψ²(problem, n)
        Ψ⁴ = abs2.(Ψ²)
        IPRs[n] = sum(Ψ⁴) / sum(Ψ²)^2
    end
    return IPRs
end
"""
    get_PRs(problem)

Participation Ratio is computed with `( ∑Ψ² )² / ∑Ψ⁴ ` for each mode
"""
function get_PRs(problem)
    IPRs = get_IPRs(problem)
    return 1.0 ./ IPRs
end

function get_number_modes(problem::LinearOptics{Scalar})
    return problem.atoms.N
end
function get_number_modes(problem::LinearOptics{Vectorial})
    return 3problem.atoms.N
end
"""
    spectrum(problem; forceComputation=false)

Returns `ωₙ` and  `Γₙ` which are imag(λ) and -real(λ) - where λ is eigenvalues of interaction matrix

Values are cached, unless `forceComputation=true`
"""
function spectrum(problem; forceComputation=false)
    if is_spectrum_NOT_available(problem) || forceComputation
        λ = eigenvalues(problem; forceComputation=forceComputation)
        ωₙ, Γₙ = imag.(λ), -real.(λ)
        if any(Γₙ .< 0)
            @warn "some Γₙ were negatives and ignored without further investigation"
            Γₙ = abs.(Γₙ)
        end        
    else
        λ = problem.spectrum["λ"]
        ωₙ, Γₙ = imag.(λ), -real.(λ)
    end

    return ωₙ, Γₙ
end

function is_spectrum_NOT_available(problem)
    if  haskey(problem.spectrum, "isSpectrumAvailable")
        return false
    else
        return true
    end
end
function make_spectrum_available(problem)
    return problem.spectrum["isSpectrumAvailable"] = true
end

function get_ψ(problem::LinearOptics, n::Integer)
    ψ = view(problem.spectrum["ψ"], :, n)
    return ψ
end
function get_ψ²(problem::LinearOptics, n::Integer)
    ψ = problem.spectrum["ψ"][:, n]
    return abs2.(ψ)
end

"""
localization_length(problem::LinearOptics; forceComputation=false, regression_method=satman2015, probability_threshold=0.999, showprogress=false, speculative=true)

'speculative = true' compute robust fitting on likely localized modes.

Set 'speculative = false' to have a better localization length estimatives for extended modes
"""
function localization_length(problem::LinearOptics; forceComputation=false, regression_method=satman2015, probability_threshold=0.999, showprogress=false, speculative=true)
    if is_localization_NOT_available(problem) && forceComputation == false
        return problem.data[:ξ], problem.data[:R2]
    end
    N = problem.atoms.N
    ξₙ = zeros(N)
    R²ₙ = zeros(N)

    # create spectrum if neeeded
    eigenvectors(problem)

    pp = Progress(N)
    Threads.@threads for n in 1:N
        DCM, ψ² = spatial_profile_single_mode(problem, n)
        ξₙ[n], R²ₙ[n] = get_single_ξ_and_R2(DCM, ψ²; regression_method=regression_method, probability_threshold=probability_threshold, speculative=speculative)

        if showprogress
            next!(pp)
        end
    end
    problem.data[:ξ], problem.data[:R2] = ξₙ, R²ₙ
    return ξₙ, R²ₙ
end

function is_localization_NOT_available(problem)
    if haskey(problem.data, :ξ) && haskey(problem.data, :R2)
        return true
    else
        return false
    end
end
"""
    spatial_profile_single_mode(problem, mode_index::Integer)

Returns `DCM, ψ²ₙ`, that is, the Distance of atoms to the center of mass, and the absolute value of the mode.
"""
function spatial_profile_single_mode(problem, mode_index::Integer)
    r = problem.atoms.r
    eigenvectors(problem) # compute diagonzalization if not available
    try
        ψ²ₙ = get_ψ²(problem, mode_index)
        r_cm = coordinates_of_center_of_mass(r, ψ²ₙ)

        # Specific for the problem (ones needs to define this function)
        DCM = [norm(rj .- r_cm) for rj in eachcol(r)]# get_Distances_from_r_to_CM

        sort_spatial_profile!(DCM, ψ²ₙ)
        return DCM, ψ²ₙ
    catch
        @error("Probably you tried to access a spectrum not created (run `spectrum(problem)`). Or `mode_index` is not valid")
    end
end
"""
    coordinates_of_center_of_mass(r, Ψ²_mode)
"""
function coordinates_of_center_of_mass(r, Ψ²_mode)
    N = length(Ψ²_mode)
    dimensions = size(r, 1)
    r_cm = zeros(dimensions)
    j = 1
    for c in eachcol(r)
        r_cm .+= c .* Ψ²_mode[j]
        j += 1
    end
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

"""
    eigenvalues(problem::LinearOptics;forceComputation=false)

Computes eigenvalues of the interaction matrix

# Example:

```julia
using CoupledDipoles
N, ρk⁻³ = 1200, 1.0
atoms = Atom(CoupledDipoles.Sphere(), sphere_inputs(N, ρk⁻³)...)

s, Δ = exp10(-5), 1.0
laser = Laser(PlaneWave3D(), s, Δ)

prob = LinearOptics(Scalar(), atoms, laser)
λ = eigenvalues(prob)
```
"""
function eigenvalues(problem::LinearOptics; forceComputation=false)
    if is_eigvals_available(problem) && !forceComputation
        return problem.spectrum["λ"]
    end

    H = interaction_matrix(problem)
    λ = eigvals(H) 
    problem.spectrum["λ"] = λ # cache results
    return λ
end

function is_eigvals_available(problem)
    if  haskey(problem.spectrum, "λ")
        return true
    else
        return false
    end
end
 
"""
    eigenvectors(problem::LinearOptics; forceComputation=false)

Computes eigenvectors of the interaction matrix

# Example:

```julia
using CoupledDipoles
N, ρk⁻³ = 1200, 1.0
atoms = Atom(CoupledDipoles.Sphere(), sphere_inputs(N, ρk⁻³)...)

s, Δ = exp10(-5), 1.0
laser = Laser(PlaneWave3D(), s, Δ)

prob = LinearOptics(Scalar(), atoms, laser)
Ψ = eigenvectors(prob)
```
"""
function eigenvectors(problem::LinearOptics;  forceComputation=false)
    if is_spectrum_NOT_available(problem) || forceComputation
        H = interaction_matrix(problem)
        eigen_values_vectors = eigen(H)
        
        problem.spectrum["λ"] = eigen_values_vectors.values
        problem.spectrum["ψ"] = eigen_values_vectors.vectors
        make_spectrum_available(problem)

        return eigen_values_vectors.vectors    
    end        
    return problem.spectrum["ψ"]
end
### --------------- Classification r---------------
"""
    classify_modes(problem)

Returns a tuple `(loc, sub, super)` with indices.
"""
function classify_modes(problem; forceComputation=false, fitting_threshold=0.5, showprogress=false, speculative=true)
    ωₙ, Γₙ = spectrum(problem; forceComputation=forceComputation)
    ξₙ, R²ₙ = localization_length(problem; forceComputation=forceComputation, showprogress=showprogress, speculative=speculative)

    localized_modes = findall((Γₙ .< Γ) .* (R²ₙ .≥ fitting_threshold))
    sub_radiant_modes = findall((Γₙ .< Γ) .* (R²ₙ .< fitting_threshold))
    super_radiant_modes = findall(Γₙ .> Γ)
    return (loc=localized_modes, sub=sub_radiant_modes, super=super_radiant_modes)
end
