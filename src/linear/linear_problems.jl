function LinearOptics(physic::Scalar, atoms, laser)
    @debug "start: LinearOptics - $( typeof(physic) )"

    kernelFunction = get_kernelFunction(physic, atoms)
    spectrum = Dict()
    data = Dict()

    @debug "end  : LinearOptics - $( typeof(physic) )"
    return LinearOptics(physic, atoms, laser, kernelFunction, spectrum, data)
end
function LinearOptics(physic::Vectorial, atoms, laser)
    @debug "start: LinearOptics - $( typeof(physic) )"

    if all(laser.polarization .==  0)
        laser.polarization = [1,0,0]
    end
    kernelFunction = get_kernelFunction(physic, atoms)
    spectrum = Dict()
    data = Dict()

    @debug "end  : LinearOptics - $( typeof(physic) )"
    return LinearOptics(physic, atoms, laser, kernelFunction, spectrum, data)
end


"""
    Laser(laser::U, s::T , Δ::T; direction=[0,0,1], polarization=[0, 0, 0])

    polarization=[0,0,0] means Scalar Physics
"""
function Laser(laser::U, s::T , Δ::T; direction=[0,0,1], polarization=[0, 0, 0]) where {U<: Pump, T}
    if dot(direction, polarization) ≈ 0
        return Laser(laser, s, Δ, direction, polarization)
    else
        @error "`direction` and `polarization` must be orthogonal vectors."
    end
end

function turn_off!(laser::U) where U<: Pump
    laser.s = zero(eltype(laser.s))
end
# get_kernelFunction(physics::Scalar,    dimension::Atom{<:TwoD}) = sin
# get_kernelFunction(physics::Vectorial, dimension::Atom{<:TwoD}) = cos

get_kernelFunction(physics::Scalar, dimension::Atom{<:ThreeD}) = green_scalar!
get_kernelFunction(physics::Vectorial, dimension::Atom{<:ThreeD}) = green_vectorial!
