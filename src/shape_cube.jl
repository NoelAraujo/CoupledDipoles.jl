#=
Examples:
    Shape(Cube(), 3, 3.6)
    Shape(Cube(), 2,3,[1,12])
    Shape(Cube(), 2,3,[1 2; 3 2])
    Shapes([Cube(), Cube()], [2,2], [2, 3], [1.0, [1 2; 3 2]])
=#
"""
    Cube Geoemtry

The cube goes from [-kL/2, kL/2] (with homogeneous distribution)
"""
function Shape(geometry::Cube, N::Int64, kL::Union{Real, Integer}; createFunction=ftn_AtomsOnCube::Function)
    @debug "creating Cube with N and kL"
    dimensions = 3
    ρ = N / kL^3
    rₘᵢₙ = get_rₘᵢₙ(ρ)
    if rₘᵢₙ ≥ kL / 10
        rₘᵢₙ = kL / 100
    end
    
    r = get_atoms(dimensions, N, rₘᵢₙ; createFunction, kL)
    return Shape(Cube(), r, N, Float64(kL))
end


function ftn_AtomsOnCube(;kwargs...)
    kL = kwargs[:kL]

    x = -kL * rand() + (kL / 2)
    y = -kL * rand() + (kL / 2)
    z = -kL * rand() + (kL / 2)
    return [x, y, z]
end