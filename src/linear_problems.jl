function LinearProblem(physic::Union{Scalar, Vectorial}, atoms, laser)
    @debug "start: LinearProblem - $( typeof(physic) )"
    
    kernelFunction = get_kernelFunction(physic, atoms)
    spectrum = Dict(:λ=>Missing, :ψ=>Missing)
    data = Dict()
    
    @debug "end  : LinearProblem - $( typeof(physic) )"
    return LinearProblem(physic, atoms, laser, kernelFunction, spectrum, data)
end

# get_kernelFunction(physics::Scalar,    dimension::Shape{<:TwoD}) = sin
# get_kernelFunction(physics::Vectorial, dimension::Shape{<:TwoD}) = cos

get_kernelFunction(physics::Scalar,    dimension::Shape{<:ThreeD}) = green_scalar!
get_kernelFunction(physics::Vectorial, dimension::Shape{<:ThreeD}) = green_vectorial!
