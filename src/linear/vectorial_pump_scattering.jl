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

    _laser_electric_fields = ThreadsX.map(eachcol(sensors)) do sensor
        (-im/2)*Ω₀*_vectorial_laser_field(problem.laser, polarization, sensor)
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
    N = size(r,2)
    n_hat = sensor./norm(sensor)
    E = zeros(Complex{eltype(r)}, 3)
    for μ = 1:3
        for j=1:N, η=1:3
            term1 = (float(μ==η) - n_hat[μ]*n_hat[η]')
            term2 = exp(-im*dot(n_hat, view(r,:,j)  ))
            term3 = β[η, j]
            E[μ] += term1*term2*term3
        end
    end
    return -(1/4π)*(exp(im*norm(sensor))/norm(sensor))*E
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
    G = zeros(ComplexF64, 3, 3)
    for j = 1:N
        rⱼ = r[:, j]
        βⱼ = β[:, j]

        _vectorial_3D_green_kernel!(sensor - rⱼ, G)
        for η=1:3
            E_x += G[1,η]*βⱼ[η]
        end
        for η=1:3
            E_y += G[2,η]*βⱼ[η]
        end
        for η=1:3
            E_z += G[3,η]*βⱼ[η]
        end



        # Xoj, Yoj, Zoj = sensor[1] - rⱼ[1], sensor[2] - rⱼ[2], sensor[3] - rⱼ[3]
        # Roj = k₀*sqrt( Xoj^2 + Yoj^2 + Zoj^2 )
        # c1 = 3*cis(Roj)/(2im*Roj^3)
        # c2 = 1im/(Roj) - 1/(Roj)^2
        # E_x = E_x + c1*(
        #      ( (Roj^2 - Xoj^2) + (Roj^2 - 3*Xoj^2)*c2)*βⱼ[1]
        #     +( ( - Xoj*Yoj) + ( - 3*Xoj*Yoj)*c2)*βⱼ[2]
        #     +( ( - Xoj*Zoj) + ( - 3*Xoj*Zoj)*c2)*βⱼ[3]
        #     )
        # E_y = E_y + c1*(
        #      ( ( - Yoj*Xoj) + ( - 3*Yoj*Xoj)*c2)*βⱼ[1]
        #     +( (Roj^2 - Yoj^2) + (Roj^2 - 3*Xoj^2)*c2)*βⱼ[2]
        #     +( ( - Xoj*Zoj) + ( - 3*Yoj*Zoj)*c2)*βⱼ[3]
        #     )
        # E_z = E_z + c1*(
        #     ( ( - Zoj*Xoj) + ( - 3*Zoj*Xoj)*c2)*βⱼ[1]
        #     +( ( - Zoj*Yoj) + ( - 3*Zoj*Yoj)*c2)*βⱼ[2]
        #     +( (Roj^2 - Zoj^2) + (Roj^2 - 3*Zoj^2)*c2)*βⱼ[3]
        # )
    end

    return [E_x, E_y, E_z]
end
function _vectorial_3D_green_kernel(r_jm::Vector)
    G = Array{Complex{eltype(r_jm)}}(undef, 3,3)
    _vectorial_3D_green_kernel!(r_jm, G)
    return G
end

function _vectorial_3D_green_kernel!(r_jm::Vector, G::Matrix)
    r = k₀*norm(r_jm)
    r2 = r^2

    #= ----- Equivalent code, but slower -----

    term1 = (3/2)cis(r)/(r)
    term2 = 1 + im/r - 1/r2
    term3 = -1 - 3im/r + 3/r2
    term4 = reshape(kron(r_jm, r_jm), 3, 3)./r2

    G = term1*(term2*I(3) + term3*term4)
    =#

    P = (3/2)*(cis(r)/(im*r))*(1 + im/r - 1/r2)
    Q = (3/2)*(cis(r)/(im*r))*(-1 - 3im/r + 3/r2)/r2

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