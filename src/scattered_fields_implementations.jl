#=
            SCATTERED FIELD: :near_field
=#
# --------------------------------- SCALAR ---------------------------------
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



# --------------------------------- VECTORIAL ---------------------------------
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
        # _E = G*βⱼ
        # E_x += _E[1]
        # E_y += _E[2]
        # E_z += _E[3]
        for η=1:3
            E_x += G[1, η]*βⱼ[η]
        end
        for η=1:3
            E_y += G[2, η]*βⱼ[η]
        end
        for η=1:3
            E_z += G[3, η]*βⱼ[η]
        end

        # r = norm(r_jm)
        # r̂ = r_jm./r

        # P = (3/2)*(cis(r)/r)
        # Q = im/r - 1/r^2

        # x, y, z = r̂[1], r̂[2], r̂[3]
        # E_x +=  P*(r^2 - x^2 + Q*(r^2 - 3*x^2))*βⱼ[1]
        # E_x += P*( -x*y - Q*3*x*y )*βⱼ[2]
        # E_x += P*( -x*z - Q*3*x*z )*βⱼ[3]

        # E_y += P*( -x*y - Q*3*x*y )*βⱼ[1]
        # E_y += P*(r^2 - y^2 + Q*(r^2 - 3*y^2))*βⱼ[2]
        # E_y += P*( -x*z - Q*3*x*z )*βⱼ[3]

        # E_z += P*( -z*x - Q*3*z*x )*βⱼ[1]
        # E_z += P*( -z*y - Q*3*z*y )*βⱼ[2]
        # E_z += P*(r^2 - z^2 + Q*(r^2 - 3*z^2))*βⱼ[3]
    end

    return +im*(Γ/2)*[E_x, E_y, E_z]
end
function _vectorial_3D_green_kernel(r_jm::Vector)
    G = Array{Complex{eltype(r_jm)}}(undef, 3,3)
    _vectorial_3D_green_kernel!(r_jm, G)
    return G
end

@inline function _vectorial_3D_green_kernel!(r_jm::Vector, G::Matrix)
    ## alternative version
    # r = norm(r_jm)
    # r̂ = r_jm ./ r

    # P = (3 / 2) * (cis(r) / r)
    # Q = im / r - 1 / r^2

    # x, y, z = r̂[1], r̂[2], r̂[3]

    # G[1] = P * (1 - x^2 + (1 - 3 * x^2)*Q)
    # G[4] = P * (-x * y - Q * 3 * x * y)
    # G[7] = P * (-x * z - Q * 3 * x * z)

    # G[2] = P * (-y * x - Q * 3 * y * x)
    # G[5] = P * (1 - y^2 + Q * (1 - 3 * y^2))
    # G[8] = P * (-y * z - Q * 3 * y * z)

    # G[3] = P * (-z * x - Q * 3 * z * x)
    # G[6] = P * (-z * y - Q * 3 * z * y)
    # G[9] = P * (1 - z^2 + Q * (1 - 3 * z^2))

    r = k₀ * norm(r_jm)
    r2 = r^2

    P = 1 + (im / r) - (1 / r2)
    Q = (-1 - (3im / r) + (3 / r2)) / r2

    x, y, z = r_jm[1], r_jm[2], r_jm[3]
    G[1] = P + x^2 * Q
    G[4] = x * y * Q
    G[7] = x * z * Q

    G[2] = y * x * Q
    G[5] = P + y^2 * Q
    G[8] = y * z * Q

    G[3] = z * x * Q
    G[6] = z * y * Q
    G[9] = P + z^2 * Q

    G .*= (3 / 2) * (cis(r) / r)
    return nothing
end







#=
            SCATTERED FIELD: :far_field
=#
# --------------------------------- SCALAR ---------------------------------
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


# --------------------------------- VECTORIAL ---------------------------------
function _vectorial_scattering_far_field(atoms::Atom{T}, β, sensor) where T <: TwoD
    return nothing
end
function _vectorial_scattering_far_field(atoms::Atom{T}, β, sensor) where T <: ThreeD
    r = atoms.r
    β_matrix = _vecAux_longArray_into_Matrix(atoms.N, β)

    norm_sensor = norm(sensor)
    n_hat = sensor./norm_sensor
    E_scatt = zeros(Complex{eltype(r)}, 3)
    terms = zeros(ComplexF64, 3)
    for μ = 1:3
        for (j, rⱼ) in enumerate(eachcol(r))
            for η=1:3
                terms[1] = (float(μ==η) - n_hat[μ]*n_hat[η]')
                terms[2] = cis( -dot(n_hat, rⱼ ))
                terms[3] = β_matrix[η, j]
                E_scatt[μ] += prod(terms)
            end
        end
    end
    return +im*(3Γ/4)*(cis(norm_sensor)/norm_sensor)*E_scatt
end