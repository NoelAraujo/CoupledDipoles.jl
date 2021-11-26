function LinearOptics(physic::Union{Scalar,Vectorial}, atoms, laser)
    @debug "start: LinearOptics - $( typeof(physic) )"

    kernelFunction = get_kernelFunction(physic, atoms)
    spectrum = 
        Dict(:λ => Array{ComplexF64}(undef, atoms.N), 
             :ψ => Array{ComplexF64}(undef, atoms.N, atoms.N),
             :isSpectrumAvailable => false)
    data = Dict()
    
    @debug "end  : LinearOptics - $( typeof(physic) )"
    return LinearOptics(physic, atoms, laser, kernelFunction, spectrum, data)
end

# get_kernelFunction(physics::Scalar,    dimension::Atom{<:TwoD}) = sin
# get_kernelFunction(physics::Vectorial, dimension::Atom{<:TwoD}) = cos

get_kernelFunction(physics::Scalar, dimension::Atom{<:ThreeD}) = green_scalar!
get_kernelFunction(physics::Vectorial, dimension::Atom{<:ThreeD}) = green_vectorial!
