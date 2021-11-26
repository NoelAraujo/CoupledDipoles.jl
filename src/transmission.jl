"""
    get_transmission(problem, β)
"""
function get_transmission(problem, β)
    new_R = how_far_is_FarField(problem.atoms)
    v_r = view(problem.atoms.r, :, :)

    (int_scatt, err) =
        hcubature(x -> _func_total(x, new_R, problem, v_r, β), (0, 0.0), (deg2rad(15), 2π))
    (int_laser, err) =
        hcubature(x -> _func_laser(x, new_R, problem, v_r, β), (0, 0.0), (deg2rad(15), 2π))
    return int_scatt / int_laser
end

function _func_total(x, new_R, problem, v_r, β)
    θ, ϕ = x[1], x[2]
    spherical_coordinate = [θ, ϕ, new_R]
    sensor = sph2cart(spherical_coordinate)
    return _get_intensity_over_sensor(problem.atoms.shape, problem.laser, v_r, sensor, β)
end
function _func_laser(x, new_R, problem, v_r, β)
    θ, ϕ = x[1], x[2]
    spherical_coordinate = [θ, ϕ, new_R]
    sensor = sph2cart(spherical_coordinate)
    return abs2(apply_laser_over_oneSensor(problem.laser, sensor))
end


function how_far_is_FarField(atoms)
    return farField_factor*size(atoms)
end