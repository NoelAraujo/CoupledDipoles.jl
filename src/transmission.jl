struct SensorType
    createSensor_func::Function
    domain::Tuple
end

function transmission(problem, β; regime=:near_field)
    if problem.laser.direction ≉ [0,0,1]
        @error "`transmission`  only works for lasers pointing into z-direction."
    end

    θₘₐₓ = problem.atoms.N < 50 ? deg2rad(45) : deg2rad(25)
    integral_domain = ((0.0, 0.0), (θₘₐₓ, 2π))

    (I_scattered, _e) = hcubature(integral_domain...) do x # , rtol=1e-12
        sensor = getSensor_on_Sphere(x, problem)
        laser_and_scattered_intensity(problem, β, sensor; regime=regime)
    end

    (I_laser, _e) = hcubature(integral_domain...) do x # , rtol=1e-12
        sensor = getSensor_on_Sphere(x, problem)
        laser_intensity(problem, sensor)
    end

    T = I_scattered / I_laser
    return T[1] # T is a matrix with 1 element
end

@inline function getSensor_on_Sphere(x, problem)
    θ, ϕ = x[1], x[2]
    new_R = 5 * size(problem.atoms) # near field values are faster to integrate
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
