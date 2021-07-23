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
function Shape(geometry::Cube, N::Integer, kL::Union{Real, Integer})#; createFunction=ftn_AtomsOnCube!::Function)
    @debug "creating Cube with N and kL"
    dimensions = 3
    ρ = N / kL^3
    rₘᵢₙ = 0.1#get_rₘᵢₙ(ρ)
    if rₘᵢₙ ≥ kL / 10
        rₘᵢₙ = kL / 100
    end
    
    r = 0.5#get_atoms(dimensions, N, rₘᵢₙ; createFunction, kL)
    return Shape(Cube(), r, N, kL)
end
"""
    Cube(r::Matrix{T}, kL::T)

I expect that: size(r) = (N atoms, 3) - 3 == N dimensions
"""
function Shape(geometry::Type{Cube}, r::Matrix{T}, kL::T) where T <: Number
    @debug "creating Cube given atoms positions and cube side"
    N = size(r, 1)
    dimensions = 3

    r_Static = SArray[]
    for n in 1:N
        push!(r_Static, SVector{dimensions}(r[n, :]))
    end
    return Shape(Cube(), r_Static, N, kL)
end
"""
    Cube(r::Vector{StaticArrays.SArray})

I expect that: length(r) = N atoms; length(r[1]) =  N dimensions
"""
function Cube(r::Vector{StaticArrays.SArray}, kL)
    N = length(r)
    return Shape(Cube(), r, N, kL)
end

function highlight(s, colour)
    io = IOBuffer()
    printstyled(IOContext(io, :color => true), s, color = colour)
    return io |> take! |> String
end

highlight("noel", :yellow) |> printstyled