struct SensorType
    createSensor_func::Function
    domain::Tuple
end

function transmission(problem, β; regime=:near_field, rtol=exp10(-10))
    # For future, use this link to create integration points in any direction: 
    # https://fasiha.github.io/post/spherical-cap/
    # https://github.com/fasiha/sphere-cap-random/blob/gh-pages/src/capRandom.js
    if problem.laser.direction ≉ [0,0,1]
        @error "`transmission`  only works for lasers pointing into z-direction."
    end
    # these values were ad hoc, but do not come from a systematic anaylis
    θₘₐₓ = problem.atoms.N < 50 ? deg2rad(50) : deg2rad(25)
    integral_domain = ((0.0, 0.0), (θₘₐₓ, 2π))
    ## The values the intensity are small
    ## making the integration convergence slow
    ## One solution is to multiply all values for some factor
    scaling_factor = 25*size(problem.atoms)

    (I_scattered, _e) = hcubature(integral_domain..., rtol=rtol) do x
        sensor = getSensor_on_Sphere(x, problem)
        scaling_factor*laser_and_scattered_intensity(problem, β, sensor; regime=regime)
    end

    (I_laser, _e) = hcubature(integral_domain..., rtol=rtol) do x
        sensor = getSensor_on_Sphere(x, problem)
        scaling_factor*laser_intensity(problem, sensor)[1]
    end

    T = I_scattered / I_laser
    return T[1] # T is a matrix with 1 element
end

@inline function getSensor_on_Sphere(x, problem)
    θ, ϕ = x[1], x[2]
    new_R = 25*size(problem.atoms) # near field values are faster to integrate
    spherical_coordinate = [θ, ϕ, new_R]
    sensor = sph2cart(spherical_coordinate)
    return sensor
end
@inline function _create_plane_sensor(x, problem)
    x, y = x[1], x[2]
    z = 5 * size(problem.atoms)
    new_z = z + 0.5 * z
    sensor = [x, y, new_z]
    return sensor
end
