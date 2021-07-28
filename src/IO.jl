Base.size(x::Shape{T}) where T = x.sizes

import Base.eltype
eltype(laser::Laser{Gaussian3D}) = "Gaussian3D"
eltype(problem::LinearProblem{Scalar}) = "Scalar"
eltype(problem::LinearProblem{Vectorial}) = "Vectorial"

function Base.show(io::IO, atoms::Shape{Cube})
    geometry_text = "Atoms on a $( highlight("Cube", :yellow) )"
    N_text = " with N=$(highlight(atoms.N, :yellow)) and "
    size_text = "side kL=$(highlight(make_short(atoms.sizes), :yellow))"
    return printstyled(io, geometry_text*N_text*size_text )
end


function Base.show(io::IO, laser::Laser{Gaussian3D})
    l_t = "$( highlight("Gaussian 3D", :yellow) )"
    w_t = " laser with waist w₀=$(highlight(laser.pump.w₀, :yellow)), "
    s_t = "s=$(highlight(make_short(laser.s), :yellow)) and "
    d_t = "Δ=$(highlight(make_short(laser.Δ), :yellow))"
    return printstyled(io, l_t*w_t*s_t*d_t )
end


function Base.show(io::IO, problem::LinearProblem{T}) where T <: Linear
    s_t = "$( highlight(eltype(problem), :yellow) ) problem "
    N_t = "with N=$(highlight(problem.atoms.N, :yellow)) atoms, "
    l_t = "and $(highlight(eltype(problem.laser), :yellow)) laser"
    return printstyled(io, s_t*N_t*l_t )
end

function highlight(s, colour)
    io = IOBuffer()
    printstyled(IOContext(io, :color => true), s, color = colour)
    return io |> take! |> String
end
function make_short(x)
    if (abs(x) ≈ 0)
        return x
    elseif (abs(x) ≥ 1e3) || (abs(x) ≤ 1e-3)
        return @sprintf("%.2E", x)
    else
        return round(x; digits=2)
    end
end
