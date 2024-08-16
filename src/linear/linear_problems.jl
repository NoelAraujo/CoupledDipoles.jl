"""
    LinearOptics(physic::Scalar, atoms, laser)
"""
function LinearOptics(physic::Scalar, atoms, laser)
    @debug "start: LinearOptics - $( typeof(physic) )"
    laser_copy = deepcopy(laser) # i don't want change the original laser
    if !all(laser_copy.polarization .==  0)
        laser_copy.polarization = [0,0,0]
    end
    kernelFunction! = get_kernelFunction!(physic, atoms)
    spectrum = Dict()
    data = Dict()

    @debug "end  : LinearOptics - $( typeof(physic) )"
    return LinearOptics(physic, atoms, laser_copy, kernelFunction!, spectrum, data)
end
"""
    LinearOptics(physic::Vectorial, atoms, laser)
"""
function LinearOptics(physic::Vectorial, atoms, laser)
    @debug "start: LinearOptics - $( typeof(physic) )"
    laser_copy = deepcopy(laser)
    # user does not specified the polarization direction
    if all(laser_copy.polarization .==  0)
        new_polarization = _produce_orthogonal_polarization(laser)
        @warn "Laser did not have a `polarization`. You got a random polarization: $(new_polarization)" maxlog = 2
        laser_copy.polarization = new_polarization
    end
    kernelFunction! = get_kernelFunction!(physic, atoms)
    spectrum = Dict()
    data = Dict()

    @debug "end  : LinearOptics - $( typeof(physic) )"
    return LinearOptics(physic, atoms, laser_copy, kernelFunction!, spectrum, data)
end


# get_kernelFunction!(physics::Scalar,    dimension::Atom{<:TwoD}) = sin
# get_kernelFunction!(physics::Vectorial, dimension::Atom{<:TwoD}) = cos

get_kernelFunction!(physics::Scalar, dimension::Atom{<:ThreeD}) = green_scalar!
get_kernelFunction!(physics::Vectorial, dimension::Atom{<:ThreeD}) = green_vectorial!

"""
    β_eff  = [all X - all Y - all Z]
"""
function _vecAux_Matrix_into_longArray(β)
    ## Ωₙ_eff  = [all X - all Y - all Z]
    β_array = vcat(view(β, 1, :), view(β, 2, :), view(β, 3, :))

    return β_array
end

"""
    β_eff = vcat(view(Ωₙ, 1, :), view(Ωₙ, 2, :), view(Ωₙ, 3, :))
"""
function _vecAux_longArray_into_Matrix(N, β)
    # transpose and NOT transpose conjugated
    # because i am just changing the array format to create
    # an effetive result
    β_x = β[1:N]
    β_y = β[(N+1):2N]
    β_z = β[(2N+1):3N]
    β_matrix = copy(transpose(hcat(β_x, β_y, β_z)))
    return β_matrix
end
