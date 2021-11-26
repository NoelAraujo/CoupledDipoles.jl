"""
    get_transmission(problem, β)
"""
function get_transmission(problem, β)
    new_R = how_far_is_FarField(problem.atoms)
    v_r = view(problem.atoms.r, :, :)
    #= 
        For small number of atoms, I need to integrate more surface area.
        Larger number of particles have stable results, and I can speed up 
        the computation using a smaller area.
        "Small number" means 50 atoms, but is a number ad hoc.
    =#
    θₘₐₓ = problem.atoms.N < 50 ? deg2rad(45) : deg2rad(15)
    ϕₘₐₓ = 2π
    (int_scatt, err) =
        hcubature(x -> _func_total(x, new_R, problem, v_r, β), (0, 0.0), (θₘₐₓ, ϕₘₐₓ))
    (int_laser, err) =
        hcubature(x -> _func_laser(x, new_R, problem, v_r, β), (0, 0.0), (θₘₐₓ, ϕₘₐₓ))
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

#=
    The value of "farField_factor*size(atoms)" is necessary to avoid
    transmission above 100% for small number of particles - numerical erros 
    goes to zero, and we get invalid results, if I allow the Far Field
    be to far, i induce such mistakes.

    The "5.1*(size(atoms)^2)" expression is based on math arguments.
    See the explanation inside benchmarks folder ('benchmarks/benchmark_farField.jl')
=#
function how_far_is_FarField(atoms)
    if atoms.N < 50
        return farField_factor*size(atoms)
    else
        return 5.01*(size(atoms)^2)
    end
end