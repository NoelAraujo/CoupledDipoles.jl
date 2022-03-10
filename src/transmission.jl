struct SensorType
    createSensor_func::Function
    domain::Tuple
end

function transmission(problem, β, scatt_func)

    θₘₐₓ = problem.atoms.N < 50 ? deg2rad(45) : deg2rad(15)
    # integral_domain = get(kwargs, :domain, ((0.0, 0.0), (θₘₐₓ, 2π)))
    integral_domain = ((0.0, 0.0), (θₘₐₓ, 2π))
    
    (I_scattered, _e) = hcubature(integral_domain...) do x # , rtol=1e-12
        sensor = getSensor_on_Sphere(x, problem)
        scattering_intensity(problem, β, sensor, scatt_func)
    end

    (I_laser, _e) = hcubature(integral_domain...) do x # , rtol=1e-12
        sensor = getSensor_on_Sphere(x, problem)
        laser_intensity(problem.laser, sensor)
    end

    T = I_scattered / I_laser
    return T
end


@inline function getSensor_on_Sphere(x, problem)
    θ, ϕ = x[1], x[2]
    new_R = how_far_is_FarField(problem.atoms)
    spherical_coordinate = [θ, ϕ, new_R]
    sensor = sph2cart(spherical_coordinate)
    return sensor
end
@inline function _create_plane_sensor(x, problem)
    x, y = x[1], x[2]
    z = size(problem.atoms)
    new_z = z + 0.5 * z
    sensor = [x, y, new_z]
    return sensor
end

#=
    The value of "farField_factor*size(atoms)" is necessary to avoid
    transmission above 100% for small number of particles - numerical erros 
    goes to zero, and we get invalid results, if I allow the Far Field
    be to far, i induce such mistakes.

    The "5.1*(size(atoms)^2)" expression is based on math arguments.
    See the explanation inside benchmarks folder ('benchmarks/benchmark_farField.jl')
=#
@inline function how_far_is_FarField(atoms)
    if atoms.N < 50
        return farField_factor * size(atoms)
    else
        return 5.01 * (size(atoms)^2)
    end
end