"""
    get_atoms(dimensions, N, rₘᵢₙ; kwargs...)

`dimensions` : 1D/2D/3D == `1`/`2`/`3`
`N` : number of atoms
`rₘᵢₙ` : radius of exclusion
"""
function get_atoms(dimensions, N, rₘᵢₙ; kwargs...)
    new_atom = zeros(dimensions)
    
    get_single_atom = kwargs[:createFunction]
    r = get_single_atom(;kwargs...)
    if N == 1
        return r
    end
    temp_atom = zeros(eltype(r), dimensions)
    nValid = 1
    for iLoop in 1:1_000_000
        temp_atom = get_single_atom(;kwargs...)
    
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
    N = size(A,2)
    distanceAb = zeros(N)
    i = 1
    for c = eachcol(A)
        distanceAb[i] = Distances.evaluate(Euclidean(), c, b)
        i += 1
    end
    return distanceAb
end

"""
    get_rₘᵢₙ(ρ) = (ρ^(-2 / 3)) / π^2
"""
get_rₘᵢₙ(ρ) = (ρ^(-1 / 3)) / (0.8π^2)  # AD HOC

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
function Atom(geometry::T, r::Matrix, kR::Union{Real, Integer}) where T <: Dimension
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

### --------------- IO ---------------
"""
If atomic position is `atomic.r` of `size(N, dimensions)`.

Select one column, effectively as : `atoms.r[:, direction]`
"""
function select_atoms_axes(atoms, direction)
    return select_matrix_axes(atoms.r, direction)
end
