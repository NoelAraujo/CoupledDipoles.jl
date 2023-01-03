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
    return LASER_FACTOR*Ω₀*_vectorial_laser_field(problem.laser, polarization, sensor)
 end

 function laser_field(problem::LinearOptics{Vectorial}, sensors::AbstractMatrix)
    Ω₀ = raby_frequency(problem.laser)
    polarization = problem.laser.polarization

    _laser_electric_fields = ThreadsX.map(eachcol(sensors)) do sensor
        LASER_FACTOR*Ω₀*_vectorial_laser_field(problem.laser, polarization, sensor)
    end
    laser_electric_fields::Matrix{ComplexF64} = hcat(_laser_electric_fields...)
    return laser_electric_fields

 end



function _get_intensity(problem::LinearOptics{Vectorial}, field::AbstractVector)
    return sum(abs2, field)
end
function _get_intensity(problem::LinearOptics{Vectorial}, fields::AbstractMatrix)
    return vec(sum(abs2, fields; dims=1))
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

    norm_sensor = norm(sensor)
    n_hat = sensor./norm_sensor
    E_scatt = zeros(Complex{eltype(r)}, 3)
    terms = zeros(ComplexF64, 3)
    for μ = 1:3
        for (j, rⱼ) in enumerate(eachcol(r))
            for η=1:3
                terms[1] = (float(μ==η) - n_hat[μ]*n_hat[η]')
                terms[2] = cis( -dot(n_hat, rⱼ ))
                terms[3] = β[η, j]
                E_scatt[μ] += prod(terms)
            end
        end
    end
    return -im*(3Γ/4)*(cis(norm_sensor)/norm_sensor)*E_scatt
end



function scattering_near_field(problem::LinearOptics{Vectorial}, β, sensor)
    _vectorial_scattering_near_field(problem.atoms, β, sensor)
end

function _vectorial_scattering_near_field(atoms::Atom{T}, β, sensor) where T <: TwoD
    return nothing
end
function _vectorial_scattering_near_field(atoms::Atom{T}, β, sensor) where T <: ThreeD
    r = atoms.r
    N = atoms.N

    E_x = zero(eltype(β))
    E_y, E_z = zero(eltype(β)), zero(eltype(β))
    r_jm = zeros(eltype(atoms.r), 3)
    G = zeros(ComplexF64, 3, 3)
    for (j, rⱼ) = enumerate(eachcol(r))
        r_jm[1] = sensor[1] - rⱼ[1]
        r_jm[2] = sensor[2] - rⱼ[2]
        r_jm[3] = sensor[3] - rⱼ[3]

        βⱼ = view(β, :, j)

        _vectorial_3D_green_kernel!(r_jm, G)
        for η=1:3
            E_x += G[1,η]*βⱼ[η]
        end
        for η=1:3
            E_y += G[2,η]*βⱼ[η]
        end
        for η=1:3
            E_z += G[3,η]*βⱼ[η]
        end
    end

    return -im*(Γ/2)*[E_x, E_y, E_z]
end
function _vectorial_3D_green_kernel(r_jm::Vector)
    G = Array{Complex{eltype(r_jm)}}(undef, 3,3)
    _vectorial_3D_green_kernel!(r_jm, G)
    return G
end

@inline function _vectorial_3D_green_kernel!(r_jm::Vector, G::Matrix)
    r = k₀*norm(r_jm)
    r2 = r^2

    P = (3/2)*(cis(r)/r)*(1 + im/r - 1/r2)
    Q = (3/2)*(cis(r)/r)*(-1 - 3im/r + 3/r2)/r2

    x, y, z = r_jm[1], r_jm[2], r_jm[3]
    G[1] = P+Q*x^2
    G[4] = Q*x*y
    G[7] = Q*x*z

    G[2] = Q*y*x
    G[5] = P + Q*y^2
    G[8] = Q*y*z

    G[3] = Q*z*x
    G[6] = Q*z*y
    G[9] = P + Q*z^2

    return nothing
end







## INTENSITY
function scattered_intensity(problem::LinearOptics{Vectorial}, atomic_states, sensor_positions; regime=:far_field)
    fields = scattered_electric_field(problem, atomic_states, sensor_positions; regime = regime)
    intesities = map(eachcol(fields)) do field
                    # sum(abs2, field) # == mapreduce(abs2, +, field) == norm(field)^2
                    _get_intensity(problem, field)
                 end
    return intesities
end