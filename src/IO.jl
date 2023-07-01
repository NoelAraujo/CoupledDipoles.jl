Base.size(x::Atom{T}) where {T} = x.sizes
Base.size(x::Atom{Cylinder}) = x.sizes[:h]

import Base.eltype
eltype(laser::Laser{PlaneWave3D}) = "PlaneWave3D"
eltype(laser::Laser{Gaussian3D}) = "Gaussian3D"

eltype(problem::LinearOptics{Scalar}) = "Scalar"
eltype(problem::LinearOptics{Vectorial}) = "Vectorial"

eltype(problem::NonLinearOptics{MeanField}) = "MeanField"
eltype(problem::NonLinearOptics{PairCorrelation}) = "PairCorrelation"

function Base.show(io::IO, atoms::Atom{Cube})
    geometry_text = "Atoms on a $( highlight("Cube", :yellow) )"
    N_text = " with N=$(highlight(atoms.N, :yellow)) and "
    size_text = "side kL=$(highlight(make_short(atoms.sizes), :yellow))"
    return printstyled(io, geometry_text * N_text * size_text)
end
function Base.show(io::IO, atoms::Atom{Sphere})
    geometry_text = "Atoms on a $( highlight("Sphere", :yellow) )"
    N_text = " with N=$(highlight(atoms.N, :yellow)) and "
    size_text = "side kR=$(highlight(make_short(atoms.sizes), :yellow))"
    return printstyled(io, geometry_text * N_text * size_text)
end
function Base.show(io::IO, atoms::Atom{Cylinder})
    geometry_text = "Atoms on a $( highlight("Cylinder", :yellow) )"
    N_text = " with N=$(highlight(atoms.N, :yellow)) and "
    size_text = "[R=$(highlight(make_short(atoms.sizes[:R]), :yellow)), h=$(highlight(make_short(atoms.sizes[:h]), :yellow))]"
    return printstyled(io, geometry_text * N_text * size_text)
end

function Base.show(io::IO, laser::Laser{PlaneWave3D})
    l_t = "$( highlight("PlaneWave 3D", :yellow) ) "
    dr_t = " laser at direction=$(highlight(laser.direction, :yellow)), "
    s_t = " s=$(highlight(make_short(laser.s), :yellow)) and "
    d_t = "Δ=$(highlight(make_short(laser.Δ), :yellow))"
    return printstyled(io, l_t * dr_t * s_t * d_t)
end

function Base.show(io::IO, laser::Laser{Gaussian3D})
    l_t = "$( highlight("Gaussian 3D", :yellow) )"
    w_t = " laser with waist w₀=$(highlight(laser.pump.w₀, :yellow)), "
    s_t = "s=$(highlight(make_short(laser.s), :yellow)) and "
    d_t = "Δ=$(highlight(make_short(laser.Δ), :yellow))"
    return printstyled(io, l_t * w_t * s_t * d_t)
end

function Base.show(io::IO, problem::LinearOptics{T}) where {T<:Linear}
    s_t = "$( highlight("LinearOptics", :blue) ): $( highlight(eltype(problem), :yellow) ) problem "
    N_t = "with N=$(highlight(problem.atoms.N, :yellow)) atoms, "
    l_t = "and $(highlight(eltype(problem.laser), :yellow)) laser"
    return printstyled(io, s_t * N_t * l_t)
end
function Base.show(io::IO, problem::NonLinearOptics{T}) where {T<:NonLinear}
    s_t = "$( highlight("NonLinearOptics", :blue) ): $( highlight(eltype(problem), :yellow) ) problem "
    N_t = "with N=$(highlight(problem.atoms.N, :yellow)) atoms, "
    l_t = "and $(highlight(eltype(problem.laser), :yellow)) laser"
    return printstyled(io, s_t * N_t * l_t)
end

function highlight(s, colour)
    io = IOBuffer()
    printstyled(IOContext(io, :color => true), s; color=colour)
    return String(take!(io))
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
