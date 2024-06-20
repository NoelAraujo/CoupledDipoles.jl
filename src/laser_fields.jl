# --------------------------------- GENERAL FUNCTIONS ---------------------------------
"""
	laser_field(problem, sensors)

Compute -(im/2)*Ω.(sensors), where Ω is the laser from the problem
"""
function laser_field(problem, sensors::AbstractMatrix)
    Ω₀ = raby_frequency(problem.laser)
    if should_return_value(problem) # that is, can return a 'Value' (Scalar) or a 'Vector' (Vectorial)
        return _VALUE_laser_field(problem, sensors, Ω₀)
    else
        return _VECTOR_laser_field(problem, sensors, Ω₀)
    end
end
## if more vectorial models will be created, change this lines
should_return_value(problem) = true
should_return_value(problem::LinearOptics{Vectorial}) = false





# ************** EDGE CASES **************
"""
	laser_field(laser::Laser{T}, sensors::AbstractMatrix) where T <: Pump

Assumes no polarizations. If they exist, should use `laser_field(problem, sensors)`
"""
function laser_field(laser::Laser{T}, sensors::AbstractMatrix) where {T<:Pump}
    Ω₀ = raby_frequency(laser)

    _dummy_size = maximum(sensors)
    _dummy_atoms = Atom(Cube(), sensors, _dummy_size)
    _dummy_problem = LinearOptics(Scalar(), _dummy_atoms, laser)
    return _VALUE_laser_field(_dummy_problem, sensors, Ω₀)
end
function laser_field(laser::Laser{T}, atoms::Atom{U}) where {T<:Pump,U<:Dimension}
    # user may want field at atoms positions
    return laser_field(laser, atoms.r)
end
function laser_field(laser::Laser{T}, sensor::AbstractVector) where {T<:Pump}
    # user may want field at single points. but we need to use Matrices and not Vectors
    @warn "should use Matrices and not Vectors" maxlog = 1
    return laser_field(laser, Array(Matrix(sensor')'))
end
function laser_field(problem, atoms::Atom{U}) where U<:Dimension
    return laser_field(problem, atoms.r)
end
# ****************************************


function _VALUE_laser_field(problem, sensors, Ω₀)
    _laser_electric_fields = ThreadsX.map(eachcol(sensors)) do sensor
        Ω₀ * laser_geometry(problem.laser, sensor)
    end

    laser_electric_fields::Matrix{ComplexF64} = hcat(_laser_electric_fields...)

    return laser_electric_fields
end

function _VECTOR_laser_field(problem::LinearOptics{Vectorial}, sensors, Ω₀)
    polarization = problem.laser.polarization
    laser_electric_fields = mapreduce(hcat, eachcol(sensors)) do sensor
        Ω₀ * polarization .* laser_geometry(problem.laser, sensor)
    end

    return laser_electric_fields
end





# --------------------------------- MODEL SPECIFIC ---------------------------------

function laser_geometry(laser::Laser{PlaneWave3D}, sensor)
    direction = view(laser.direction, :)
    Eᵢ = cis(+k₀ * (sensor[1] * direction[1] + sensor[2] * direction[2] + sensor[3] * direction[3]))
    return Eᵢ
end
function laser_geometry(laser::Laser{Gaussian3D}, sensor)
    w₀ = laser.pump.w₀
    k = laser.direction

    k_dot_r = dot(sensor, k)
    distance_from_axes = (abs2(norm(sensor)) - abs2(dot(sensor, k)) / abs2(norm(k)))
    denominator_factor = 1 + 2im * k_dot_r / (k₀ * w₀^2)

    Eᵢ = cis(+k₀ * k_dot_r)
    Eᵢ = Eᵢ * exp(-distance_from_axes / (denominator_factor * w₀^2))
    Eᵢ = Eᵢ / denominator_factor
    return Eᵢ
end


"""
	`_scalar_laser_field` is legacy code. Use `laser_geometry` instead.
"""
function _scalar_laser_field(laser, sensor)
    @warn "Legacy code. Use `laser_geometry` intead" maxlog = 1
    laser_geometry(laser, sensor)
end
