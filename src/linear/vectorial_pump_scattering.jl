function laser_field(problem::LinearOptics{Vectorial}, atoms::Atom{T}) where T <: Dimension
    return laser_field(problem, atoms.r)
end

"""
    laser_field(::LinearOptics{Vectorial}, sensor) = (-im/2)*Ω₀*[laser geometry]

    returns a matrix (rows, cols) = (atoms, xyz electric-field components)
"""
function laser_field(problem::LinearOptics{Vectorial}, sensor::AbstractVector)
    Ω₀ = raby_frequency(problem.laser)
    polarization = problem.laser.polarization
    return (-im/2)*Ω₀*_vectorial_laser_field(problem.laser, polarization, sensor)
 end

 function laser_field(problem::LinearOptics{Vectorial}, sensors::AbstractMatrix)
    Ω₀ = raby_frequency(problem.laser)
    polarization = problem.laser.polarization
    return mapreduce(hcat, eachcol(sensors)) do sensor
                (-im/2)*Ω₀*_vectorial_laser_field(problem.laser, polarization, sensor)
            end
 end


 ## PUMP FIELD
function _vectorial_laser_field(laser::Laser{PlaneWave3D}, polarization, sensor)
    Eᵢ = _scalar_laser_field(laser, sensor)

    return Eᵢ.*polarization
end
function _vectorial_laser_field(laser::Laser{Gaussian3D}, polarization, sensor)
    Eᵢ = _scalar_laser_field(laser, sensor)

    return Eᵢ.*polarization
end



## SCATTERING FIELD
function scattering_far_field(problem::LinearOptics{Vectorial}, β, sensor)
    _vectorial_scattering_far_field(problem.atoms, β, sensor)
end

function _vectorial_scattering_far_field(atoms::Atom{T}, β, sensor) where T <: TwoD
    return nothing
end
function _vectorial_scattering_far_field(atoms::Atom{T}, β, sensor) where T <: ThreeD
    r = atoms.r
    N = size(r,2)
    sensor_distance = norm(sensor)
    r_versor = sensor./sensor_distance
    E = zeros(Complex{eltype(r)}, 3)
    for μ = 1:3
        for j=1:N, η=1:3
            term1 = cis(sensor_distance)/(im*sensor_distance)
            term2 = float(μ == η) - r_versor[μ]*(r_versor[η]')
            term3 = cis(- r_versor[1]*r[1] - r_versor[2]*r[2] - r_versor[3]*r[3] )
            term4 = β[η, j]
            E[μ] += term1*term2*term3*term4
        end
    end
    return -im*E
end



function scattering_near_field(problem::LinearOptics{Vectorial}, β, sensor)
    _vectorial_scattering_near_field(problem.atoms, β, sensor)
end

function _vectorial_scattering_near_field(atoms::Atom{T}, β, sensor) where T <: TwoD
    return nothing
end
function _vectorial_scattering_near_field(atoms::Atom{T}, β, sensor) where T <: ThreeD
    r = atoms.r
    E = mapreduce(+, pairs(eachcol(r))) do atom
        j, rⱼ = atom
        _vectorial_3D_green_kernel(sensor - rⱼ)*view(β, :, j)
    end
    return -im*(k₀^3/(6π))*E
end
function _vectorial_3D_green_kernel(r_jm::Vector)
    r = k₀*norm(r_jm)
    r2 = r^2

    term1 = (3/2)cis(r)/(im*r)
    term2 = 1 + im/r - 1/r2
    term3 = -1 - 3im/r + 3/r2

    term4 = reshape(kron(r_jm, r_jm), 3, 3)./r2
    G = term1*(term2*I(3) + term3*term4)

    return G
end









## INTENSITY
function scattered_intensity(problem::LinearOptics{Vectorial}, atomic_states, sensor_positions; regime=:far_field)
    fields = scattered_electric_field(problem, atomic_states, sensor_positions; regime = regime)
    intesities = map(eachcol(fields)) do field
                    norm(field)
                 end
    return intesities
end