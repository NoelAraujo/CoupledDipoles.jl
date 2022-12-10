function LinearOptics(physic::Scalar, atoms, laser)
    @debug "start: LinearOptics - $( typeof(physic) )"
    laser_copy = deepcopy(laser) # i don't want change the original laser
    if !all(laser_copy.polarization .==  0)
        laser_copy.polarization = [0,0,0]
    end
    kernelFunction = get_kernelFunction(physic, atoms)
    spectrum = Dict()
    data = Dict()

    @debug "end  : LinearOptics - $( typeof(physic) )"
    return LinearOptics(physic, atoms, laser_copy, kernelFunction, spectrum, data)
end
function LinearOptics(physic::Vectorial, atoms, laser)
    @debug "start: LinearOptics - $( typeof(physic) )"
    laser_copy = deepcopy(laser)
    # user does not specified the polarization direction
    if all(laser_copy.polarization .==  0)
        new_polarization = _produce_orthogonal_polarization(laser)
        @warn "Laser did not have a `polarization`. You got a random polarization: $(new_polarization)" maxlog = 2
        laser_copy.polarization = new_polarization
    end
    kernelFunction = get_kernelFunction(physic, atoms)
    spectrum = Dict()
    data = Dict()

    @debug "end  : LinearOptics - $( typeof(physic) )"
    return LinearOptics(physic, atoms, laser_copy, kernelFunction, spectrum, data)
end


# get_kernelFunction(physics::Scalar,    dimension::Atom{<:TwoD}) = sin
# get_kernelFunction(physics::Vectorial, dimension::Atom{<:TwoD}) = cos

get_kernelFunction(physics::Scalar, dimension::Atom{<:ThreeD}) = green_scalar!
get_kernelFunction(physics::Vectorial, dimension::Atom{<:ThreeD}) = green_vectorial!
