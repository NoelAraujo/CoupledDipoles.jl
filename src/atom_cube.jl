"""
    Cube Geoemtry

The cube goes from [-kL/2, kL/2] (with homogeneous distribution)
"""
function Atom(
    geometry::Cube,
    N::Int64,
    kL::Union{Real,Integer};
    createFunction = ftn_AtomsOnCube::Function,
)
    @debug "start: Shape - Cube"

    dimensions = 3
    ρ = N / kL^3
    rₘᵢₙ = get_rₘᵢₙ(ρ)
    if rₘᵢₙ ≥ kL / 10
        rₘᵢₙ = kL / 100
    end

    r = get_atoms(dimensions, N, rₘᵢₙ; createFunction, kL)

    @debug "end  : Shape - Cube"
    return Atom(Cube(), r, N, Float64(kL))
end

function Atom(geometry::Cube, r::Matrix, kL::Union{Real,Integer})
    N = size(r, 2) # remember to use each collum as a atom position
    return Atom(Cube(), r, N, Float64(kL))
end

function ftn_AtomsOnCube(; kwargs...)
    kL = kwargs[:kL]

    x = -kL * rand() + (kL / 2)
    y = -kL * rand() + (kL / 2)
    z = -kL * rand() + (kL / 2)
    return [x, y, z]
end