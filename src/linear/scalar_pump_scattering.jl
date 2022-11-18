function laser_field(problem::LinearOptics{Scalar}, atoms::Atom{T}) where T <: Dimension
    return laser_field(problem, atoms.r)
end

"""
    laser_field(problem, sensor) = (-im/2)*Ω₀*[laser geometry]
"""
function laser_field(problem::LinearOptics{Scalar}, sensor::AbstractVector)
    Ω₀ = raby_frequency(problem.laser)
    return (-im/2)*Ω₀*_scalar_laser_field(problem.laser, sensor)
end
function laser_field(problem::LinearOptics{Scalar}, sensors::AbstractMatrix)
    Ω₀ = raby_frequency(problem.laser)
    return map(eachcol(sensors)) do sensor
                (-im/2)*Ω₀*_scalar_laser_field(problem.laser, sensor)
           end
end




"""
    raby_frequency(laser) = √(s * (1 + 4(Δ / Γ)^2) / 2)
"""
function raby_frequency(laser)
    return √(laser.s * (1 + 4(laser.Δ / Γ)^2) / 2)
end




## PUMP FIELD
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




## SCATTERING FIELD: FAR FIELD
"""
    Eletric Field: -(Γ/2) * (exp(ikr) / ikr) * ∑ⱼ exp(-i*k₀* n̂⋅R⃗ⱼ)
"""
function scattering_far_field(problem::LinearOptics{Scalar}, β, sensor)
    _scattering_far_field(problem.atoms, β, sensor)
end

function _scattering_far_field(atoms::Atom{T}, β, sensor)  where T <: TwoD
    ## TO DO
    return nothing
end
function _scattering_far_field(atoms::Atom{T}, β, sensor) where T <: ThreeD
    atoms = atoms.r
    n̂ = sensor / norm(sensor)

    E_scatt = mapreduce(+, pairs(eachcol(atoms))) do x
        j, atom = x
        cis(-k₀ * n̂[1] * atom[1] + n̂[2] * atom[2] + n̂[3] * atom[3]) * β[j]
    end

    ikr = im * k₀ * norm(sensor)
    E_scatt = -(Γ / 2) * E_scatt * exp(ikr) / ikr
    return E_scatt
end




## SCATTERING FIELD: NEAR FIELD
"""
    Eletric Field: +(Γ/2) * ∑ⱼ exp(-i*k₀* n̂⋅R⃗ⱼ)/(k₀* sensor⋅R⃗ⱼ)
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
        d_SensorAtom = sqrt((sensor[1] - atom[1])^2 + (sensor[2] - atom[2])^2 + (sensor[3] - atom[3])^2)
        cis(k₀ * d_SensorAtom) * (β[j] / d_SensorAtom)
    end
    E_scatt = +(Γ / 2) * im * E_scatt
    return E_scatt
end





## INTENSITY
function scattered_intensity(problem::LinearOptics{Scalar}, atomic_states, sensor_positions; regime=:far_field)
    fields = scattered_electric_field(problem, atomic_states, sensor_positions; regime = regime)
    intesities = abs2.(fields)
    return intesities
end