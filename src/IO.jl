Base.size(x::Shape{T}) where T = x.sizes

function Base.show(io::IO, atoms::Shape{Cube})
    geometry_text = "Atoms on a $( highlight("Cube", :yellow) )"
    N_text = " with N=$(highlight(atoms.N, :yellow)) and "
    size_text = "side kL=$(highlight(make_short(atoms.sizes), :yellow))"
    return printstyled(io, geometry_text*N_text*size_text )
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
