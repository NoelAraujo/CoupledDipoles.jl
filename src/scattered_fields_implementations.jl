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
            E_x += G[1,η]*βⱼ[η]
        end
        for η=1:3
            E_y += G[2,η]*βⱼ[η]
        end
        for η=1:3
            E_z += G[3,η]*βⱼ[η]
        end
    end

    return +im*(Γ/2)*[E_x, E_y, E_z]
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