struct SensorType
    createSensor_func::Function
    domain::Tuple
end

function get_transmission(problem, β; kwargs...)

    create_sensors_func = get(kwargs, :create_sensors_func, _create_sphere_sensor)
    θₘₐₓ = problem.atoms.N < 50 ? deg2rad(45) : deg2rad(15)
    integral_domain = get(kwargs, :domain, ((0.0, 0.0), (θₘₐₓ, 2π)))

    sensor_setup = SensorType(create_sensors_func, integral_domain)

    scattering = get(kwargs, :scattering, :farField)
    I_scattered = _get_intensity_laser_plus_scattering(problem,β; scattering = scattering, sensor_setup = sensor_setup)
    I_laser =     _get_intensity_laser(problem; sensor_setup)
    T = I_scattered / I_laser
    return T
end


function _get_intensity_laser_plus_scattering(
    problem,
    β;
    scattering::Symbol,
    sensor_setup::SensorType,
)
    if scattering == :nearField
        scattering_func = _scattering_nearField
    elseif scattering == :farField
        scattering_func = _scattering_farField
    else
        @error "Invalid scattering Value: should be `:nearField` or `:farField` "
    end

    view_of_atoms = view(problem.atoms.r, :, :)

    _toCall_I_total(x) = _I_total(
        x,
        problem,
        view_of_atoms,
        β,
        scattering_func,
        sensor_setup.createSensor_func,
    )
    (int_scatt, _e) = hcubature(_toCall_I_total, sensor_setup.domain...)
    return int_scatt
end
function _I_total(
    x,
    problem::LinearOptics,
    v_r,
    β::Vector{ComplexF64},
    scattering_func::Function,
    getSensor::Function,
)

    sensor = getSensor(x, problem)

    return _get_intensity_over_sensor(
        problem.atoms.shape,
        problem.laser,
        v_r,
        sensor,
        β;
        scattering_func = scattering_func,
    )
end


function _get_intensity_laser(problem; sensor_setup::SensorType)
    integral_args = sensor_setup.domain
    myargs = (:problem => problem, :createSensor_func => sensor_setup.createSensor_func)

    (int_laser, _e) = hcubature(x -> _I_laser(x; myargs...), integral_args...)
    return int_laser
end
function _I_laser(x; kwargs...)
    problem = kwargs[:problem]
    getSensor = kwargs[:createSensor_func]
    sensor = getSensor(x, problem)

    return abs2(apply_laser_over_oneSensor(problem.laser, sensor))
end


@inline function _create_sphere_sensor(x, problem)
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