"""
    get_atoms(dimensions, N, rₘᵢₙ; kwargs...)

`dimensions` : 1D/2D/3D == `1`/`2`/`3`
`N` : number of atoms
`rₘᵢₙ` : radius of exclusion
"""
function get_atoms(dimensions, N, rₘᵢₙ; kwargs...)
    get_single_atom = kwargs[:createFunction]
    r = get_single_atom(; kwargs...)
    if N == 1
        return r
    end
    temp_atom = zeros(eltype(r), dimensions)
    nValid = 1
    for iLoop in 1:1_000_000
        temp_atom = get_single_atom(; kwargs...)

        if is_atom_valid(temp_atom, r, rₘᵢₙ)
            nValid = nValid + 1
            r = hcat(r, temp_atom)
        end
        if nValid == N
            break
        end
    end
    if nValid < N
        @error "Could not generate all data after 10^6 interactions"
    end

    return r
end

"""
    Checks if the `new atom` has a `minimum distance` for `all the other atoms`
"""
function is_atom_valid(new_atom, r, rₘᵢₙ)
    A = r
    b = new_atom
    allDistances = get_Distance_A_to_b(A, b)
    return all(allDistances .≥ rₘᵢₙ)
end

function get_Distance_A_to_b(A, b)
    N = size(A, 2)
    distanceAb = zeros(N)
    i = 1
    for c in eachcol(A)
        distanceAb[i] = Distances.evaluate(Euclidean(), c, b)
        i += 1
    end
    return distanceAb
end

"""
    get_rₘᵢₙ(ρ)

AD HOC value that changes many sometimes to fix different problems
"""
get_rₘᵢₙ(ρ) = (ρ^(-1 / 3)) / (0.8π^2)

"""
    get_pairwise_matrix(r)
r = matrix with `N` atoms positions

returns a `NxN` Float64 matrix, with zeros at diagonal
"""
function get_pairwise_matrix(r)
    r_matrix = transpose(r)
    R_jk = Distances.pairwise(Euclidean(), r_matrix, r_matrix; dims=1)
    R_jk[diagind(R_jk)] .= 0
    return R_jk
end

### -----SHAPE CONSTRUCTOR GIVEN MATRIX -------
function Atom(geometry::T, r::Matrix, kR::Union{Real,Integer}) where {T<:Dimension}
    N = size(r, 2) # remember to use each collum as a atom position
    return Atom(geometry, r, N, Float64(kR))
end

### --------------- CONVERSIONS ---------------
"""
    cube_inputs(N::Integer, ρ::Real)
"""
function cube_inputs(N::Integer, ρ::Real)
    kL = (N / ρ)^(1 / 3)
    return (N, kL)
end
"""
    sphere_inputs(N::Integer, ρ::Real)
"""
function sphere_inputs(N::Integer, ρ::Real)
    kR = (N / ((4π / 3) * ρ))^(1 / 3)
    return (N, kR)
end
"""
    cylinder_inputs(N::Integer, ρ::Real; R::Real, h::Real)
    
- if the use does not specify, R or h, than, it is assumed that both should be equal  
- if use only specify R or h, the other variable will change to match the density required
"""
function cylinder_inputs(N::Integer, ρ::Real; kwargs...)
    R = get(kwargs, :R, NaN)
    h = get(kwargs, :h, NaN)

    if (R < 0) || (h < 0)
        @error "Invalid inputs, R and h must be positives"
        R = NaN
        h = NaN
    elseif isnan(R) && isnan(h) # R and h not provided
        # assume that R=h ⇢ (h*πR^2 = πR^3) and solve for R
        R = cbrt(N / π * ρ)
        h = R
    elseif (h > 0) && isnan(R) # R not provided
        R = sqrt(N / (π * ρ * h))
    elseif (R > 0) && isnan(h) # h not provided
        h = N / (π * ρ * R^2)
    end
    return (N, R, h)
end
