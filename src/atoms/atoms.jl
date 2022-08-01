"""
    get_atoms(dimensions, N, rₘᵢₙ; kwargs...)

`dimensions` : 1D/2D/3D == `1`/`2`/`3`
`N` : number of atoms
`rₘᵢₙ` : radius of exclusion
"""
function get_atoms(dimensions, N, rₘᵢₙ; kwargs...)
    get_single_atom = kwargs[:createFunction]
    r_single = get_single_atom(; kwargs...)
    if N == 1
        return r_single
    end
    r = Array{eltype(rₘᵢₙ)}(undef, dimensions, N)
    r[:, 1] .= r_single

    temp_atom = zeros(eltype(r), dimensions)
    nValid = 1
    for iLoop in 1:1_000_000
        temp_atom = get_single_atom(; kwargs...)

        if is_atom_valid(temp_atom, r, rₘᵢₙ, nValid)
            nValid = nValid + 1
            r[:, nValid] .= temp_atom
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
function is_atom_valid(new_atom, r, rₘᵢₙ, N)
    Nr = length(new_atom)
    for i = 1:N
        c = @view r[:, i]
        if Nr == 3
            distanceAb = sqrt((c[1] - new_atom[1])^2 + (c[2] - new_atom[2])^2 + (c[3] - new_atom[3])^2)
            if distanceAb < rₘᵢₙ
                return false
            end
        elseif Nr == 2
            distanceAb = sqrt((c[1] - new_atom[1])^2 + (c[2] - new_atom[2])^2)
            if distanceAb < rₘᵢₙ
                return false
            end
        end
    end
    return true
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
function get_pairwise_matrix!(r, R_jk)
    r_matrix = transpose(r)
    Distances.pairwise!(R_jk, Euclidean(), r_matrix, r_matrix; dims=1)
    R_jk[diagind(R_jk)] .= 0
    return nothing
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
