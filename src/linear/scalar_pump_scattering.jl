#=
            LASER FIELD
=#
function laser_field(problem::LinearOptics{Scalar}, sensor::AbstractVector)
    Ω₀ = raby_frequency(problem.laser)
    return Matrix(transpose([LASER_FACTOR*Ω₀*_scalar_laser_field(problem.laser, sensor)]))
end
function laser_field(problem::LinearOptics{Scalar}, sensors::AbstractMatrix)
    Ω₀ = raby_frequency(problem.laser)

    _laser_electric_fields = ThreadsX.map(eachcol(sensors)) do sensor
            LASER_FACTOR*Ω₀*_scalar_laser_field(problem.laser, sensor)
    end
    laser_electric_fields::Matrix{ComplexF64} = hcat(_laser_electric_fields...)
    return laser_electric_fields
end

## Matrix case happens for single sensor field
function _get_intensity(problem::LinearOptics{Scalar}, field::AbstractMatrix)
    return vec(abs2.(field))
end
function _get_intensity(problem::LinearOptics{Scalar}, fields::AbstractVector)
    return abs2.(fields)
end



"""
    _scalar_laser_field(laser::Laser{PlaneWave3D}, sensor)

Compute Plane Wave in Any Direction: exp( im* laser.direction ⋅ sensor )
"""
function _scalar_laser_field(laser::Laser{PlaneWave3D}, sensor)
    direction = view(laser.direction, :)
    Eᵢ = cis(+k₀ * (sensor[1] * direction[1] + sensor[2] * direction[2] + sensor[3] * direction[3]))
    return Eᵢ
end

"""
    _scalar_laser_field(laser::Laser{Gaussian3D}, sensor)

Compute Paraxial Gaussian Beam in Any Direction
"""
function _scalar_laser_field(laser::Laser{Gaussian3D}, sensor)
    w₀ = laser.pump.w₀
    k = laser.direction

    k_dot_r = dot(sensor, k)
    distance_from_axes = (abs2(norm(sensor)) - abs2( dot(sensor, k) )/abs2(norm(k)))
    denominator_factor = 1 + 2im * k_dot_r / (k₀ * w₀^2)

    Eᵢ = cis(+k₀ * k_dot_r)
    Eᵢ = Eᵢ * exp(- distance_from_axes / (denominator_factor * w₀^2))
    Eᵢ = Eᵢ / denominator_factor
    return Eᵢ
end


#=
            SCATTERED FIELD: :near_field
=#
"""
    Eletric Field: +i * (Γ/2) * ∑ⱼ exp(i*k₀* |R-R⃗ⱼ|)/(k₀*|R-R⃗ⱼ|)
"""
function scattering_near_field(problem::LinearOptics{Scalar}, β, sensor)
    _scalar_scattering_near_field(problem.atoms, β, sensor)
end

function _scalar_scattering_near_field(atoms::Atom{T}, β, sensor)  where T <: TwoD
    ## TO DO
    return nothing
end
function _scalar_scattering_near_field(atoms::Atom{T}, β, sensor) where T <: ThreeD
    atoms = atoms.r

    E_scatt = mapreduce(+, pairs(eachcol(atoms))) do x
        j, atom = x
        d_SensorAtom = k₀*sqrt((sensor[1] - atom[1])^2 + (sensor[2] - atom[2])^2 + (sensor[3] - atom[3])^2)
        (cis(d_SensorAtom) / d_SensorAtom)*β[j]
    end

    E_scatt = +im*(Γ / 2) * E_scatt
    return E_scatt
end

#=
            SCATTERED FIELD: :far_field
=#
"""
    Eletric Field: +i * (Γ/2) * (exp(ikR) / kR) * ∑ⱼ exp(-i*k₀* n̂⋅R⃗ⱼ)
"""
function scattering_far_field(problem::LinearOptics{Scalar}, β, sensor)
    _scalar_scattering_far_field(problem.atoms, β, sensor)
end

function _scalar_scattering_far_field(atoms::Atom{T}, β, sensor)  where T <: TwoD
    ## TO DO
    return nothing
end
function _scalar_scattering_far_field(atoms::Atom{T}, β, sensor) where T <: ThreeD
    r = atoms.r
    sensor_norm = k₀*norm(sensor)
    n̂ = sensor / sensor_norm

    E_scatt = mapreduce(+, pairs(eachcol(r))) do x
        j, atom = x
        cis(-(n̂[1]*atom[1] + n̂[2]*atom[2] + n̂[3]*atom[3])) * β[j]
    end
    R = how_far_is_farField(atoms)
    E_scatt = +im*(Γ / 2 ) * (cis(R) / R) * E_scatt
    return E_scatt
end