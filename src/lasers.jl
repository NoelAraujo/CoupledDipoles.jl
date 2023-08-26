"""
    Laser(laser::U, s::T , Δ::T; direction=[0,0,1], polarization=[0, 0, 0])
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
function turn_laser_off!(problem)
    problem.laser.s = zero(eltype(problem.laser.s))
    nothing
end


function _produce_orthogonal_polarization(laser)
    direction = laser.direction
    T = eltype(direction)

    #=
     The default direction is on z-axis,
     thus, a good guess for polarization is on x-axis.

     If is not the case, i create some "random number":

     1. For reproducebility, the seed will be proportional to the
     laser direction.
        Because i don't need something complex for the seed:
        - it must be positive (this explains "^4" and not "^3")
        -  i need to be different than zero when rounding, otherwise
            all values will be same (this explains the "10")

    2. I need to orthonormalize it with Gram-Schmidt.
    =#
    polarization = [one(T), zero(T), zero(T)]
    if direction ⋅ polarization ≈ 0
         return polarization
    else
        aSeed = round(Int64, sum(10direction)^4)
        random_polarization = rand(MersenneTwister(aSeed), T, length(direction))

        v = random_polarization - (direction'*random_polarization)*direction
        new_polarization = v./norm(v)
        return new_polarization
    end
 end



 """
 laser_intensity(laser, points)

 |laser_field(laser, points)|^2
 """
 function laser_intensity(problem, sensors)
     E_laser = laser_field(problem, sensors)
     intesities = mapreduce(vcat, eachcol(E_laser)) do field
         _get_intensity(problem, field)
     end
     return intesities
 end

"""
    raby_frequency(laser) = Γ √(s / 2) # saturation is on ressonance
"""
function raby_frequency(laser)
    return Γ * √(laser.s / 2)
end