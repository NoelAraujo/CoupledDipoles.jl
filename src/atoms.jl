"""
    get_atoms(dimensions, N, rₘᵢₙ; kwargs...)

`dimensions` : 1D/2D/3D == `1`/`2`/`3`
`N` : number of atoms
`rₘᵢₙ` : radius of exclusion
"""
function get_atoms(dimensions, N, rₘᵢₙ; kwargs...)
    new_atom = zeros(dimensions)
    r = SArray[]

    ftn_sampling! = kwargs[:createFunction]

    ftn_sampling!(new_atom; kwargs...)
    push!(r, SVector{dimensions}(new_atom))
    if N == 1
        return r
    end
    nValid = 1
    for iLoop in 1:1_000_000
        ftn_sampling!(new_atom; kwargs...)

        if is_atom_valid(new_atom, r, rₘᵢₙ)
            nValid = nValid + 1
            push!(r, SVector{dimensions}(new_atom))
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
    N = length(A)
    distanceAb = zeros(N)
    Threads.@threads for i in 1:N
        distanceAb[i] = Distances.evaluate(Euclidean(), A[i], b)
    end
    return distanceAb
end
"""
    get_rₘᵢₙ(ρ) = (ρ^(-2 / 3)) / π^2
"""
get_rₘᵢₙ(ρ) = (ρ^(-2 / 3)) / π^2  # AD HOC
"""
    get_pairwise_matrix(r)
r = matrix with `N` atoms positions

returns a `NxN` Float64 matrix, with zeros at diagonal
"""
function get_pairwise_matrix(r)
    r_matrix = convert_StaticArray_to_matrix(r)
    R_jk = Distances.pairwise(Euclidean(), r_matrix, r_matrix; dims=1)
    R_jk[diagind(R_jk)] .= 0
    return R_jk
end

function convert_StaticArray_to_matrix(r::Vector{StaticArrays.SArray})
    N = length(r)
    r_matrix = zeros(N, 3)
    for n in 1:N
        r_matrix[n, :] = [r[n][1] r[n][2] r[n][3]]
    end
    return r_matrix
end
function convert_StaticArray_to_matrix(r::Matrix)
    return r
end

function convert_matrix_to_StaticArray(r::Matrix)
    N = size(r, 1)
    dimensions = size(r, 2)

    rs = SArray[]
    for i=1:N
        push!(rs, SVector{dimensions}(r[i,:]))
    end
    return rs
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

### --------------- SHAPES ---------------
function ftn_AtomsOnCube!(new_atom; kwargs...)
    kL = kwargs[:kL]
    new_atom[1] = -kL * rand() + (kL / 2)
    new_atom[2] = -kL * rand() + (kL / 2)
    new_atom[3] = -kL * rand() + (kL / 2)
    return nothing
end

function ftn_AtomsOnSphere!(new_atom; kwargs...)
    kR = kwargs[:kR]
    new_atom[1] = 2π * rand() # azimuth
    new_atom[2] = asin(2 * rand() - 1) # elevation
    new_atom[3] = kR * (rand()^(1 ./ 3.0)) # radii
    new_atom[:] = sph2cart(new_atom)
    return nothing
end
"""
    spherical_coordinate=[azimuth, elevation, r]

    azimuth = projection on XY-plane (in radians)
    "elevation" or "Polar" = projection on Z-axis (in radians)
    r = radius

	ref: https://www.mathworks.com/help/matlab/ref/sph2cart.html
"""
function sph2cart(spherical_coordinate)
    azimuth = spherical_coordinate[1]
    elevation = spherical_coordinate[2]
    radius = spherical_coordinate[3]
    x = radius * cos(elevation) * cos(azimuth)
    y = radius * cos(elevation) * sin(azimuth)
    z = radius * sin(elevation)
    return [x, y, z]
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

"""
If atomic position is `atomic.r` of `size(N, dimensions)`.

Select one row, effectively as : `atoms.r[n, :]`
"""
function get_one_atom(atoms::T where {T<:ThreeD}, n)
    atoms.r[n]
end
