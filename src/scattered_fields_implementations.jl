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
    
    E_scatt = zeros(eltype(β), 3)
    r_jm = zeros(eltype(atoms.r), 3)
    G = zeros(ComplexF64, 3, 3)
    for (j, rⱼ) = enumerate(eachcol(r))       
        βⱼ = view(β, :, j)

        Xoj = sensor[1] .- rⱼ[1]
        Yoj = sensor[2] .- rⱼ[2]
        Zoj = sensor[3] .- rⱼ[3]
        Roj = sqrt(Xoj^2 + Yoj^2 + Zoj^2)

        # Vectorial field
        c1 = cis(Roj)./(1im*Roj^3)
        c2 = 1im/Roj - 1/(Roj^2)
        
        betax = βⱼ[1]
        betay = βⱼ[2]
        betaz = βⱼ[3]

        E_scatt[1] -= c1 .* (
                ( Roj.^2 - Xoj.^2 + c2.* (Roj.^2 - 3* Xoj.^2) ) .* betax +
                ( -Xoj.* Yoj.- 3* c2.* Xoj.* Yoj ) .* betay +
                ( -Xoj.* Zoj.- 3* c2.* Xoj.* Zoj ) .* betaz
            )

        E_scatt[2] -= c1 .* (
                ( Roj.^2 .- Yoj.^2 .+ c2.* (Roj.^2 .- 3* Yoj.^2) ) .* betay +
                ( -Xoj.* Yoj.- 3* c2.* Xoj.* Yoj ) .* betax +
                ( -Yoj.* Zoj.- 3* c2.* Yoj.* Zoj ) .* betaz
            )
            

        E_scatt[3] -= c1 .* (
                ( Roj.^2 .- Zoj.^2 .+ c2.* (Roj.^2 .- 3* Zoj.^2) ) .* betaz +
                ( -Xoj.* Zoj .- 3* c2.* Xoj.* Zoj ) .* betax +
                ( -Yoj.* Zoj .- 3* c2.* Yoj.* Zoj ) .* betay
            )            
    end
    return -im*3Γ/2*E_scatt
end
function _vectorial_3D_green_kernel(r_jm::Vector)
    G = Array{Complex{eltype(r_jm)}}(undef, 3,3)
    _vectorial_3D_green_kernel!(r_jm, G)
    return G
end


@inline function _vectorial_3D_green_kernel!(r_jm::Vector, G::Matrix)
#   reverse engineered (and adapted) from Ana Cipris
    r = k₀ * norm(r_jm)
    r2 = r^2
        
    ### v1
    # c1 = (3/2)*cis(r)/r
    # c2 = 1im/r - 1/(r^2)
    
    # n_jm = r_jm./r
    # n_x, n_y, n_z = n_jm[1], n_jm[2], n_jm[3]

    # G[1,1] = c1*( (1 - n_x*n_x)  + (1 - 3*n_x*n_x)*c2  ) # Gxx
    # G[1,2] = c1*( (0 - n_x*n_y)  + (0 - 3*n_x*n_y)*c2  ) # Gxy
    # G[1,3] = c1*( (0 - n_x*n_z)  + (0 - 3*n_x*n_z)*c2  ) # Gxz

    # G[2,1] = c1*( (0 - n_y*n_x)  + (0 - 3*n_y*n_x)*c2  ) # Gyx
    # G[2,2] = c1*( (1 - n_y*n_y)  + (1 - 3*n_y*n_y)*c2  ) # Gyy
    # G[2,3] = c1*( (0 - n_y*n_z)  + (0 - 3*n_y*n_z)*c2  ) # Gyz

    # G[3,1] = c1*( (0 - n_z*n_x)  + (0 - 3*n_z*n_x)*c2  ) # Gzx
    # G[3,2] = c1*( (0 - n_z*n_y)  + (0 - 3*n_z*n_y)*c2  ) # Gzy
    # G[3,3] = c1*( (1 - n_z*n_z)  + (1 - 3*n_z*n_z)*c2  ) # Gzz

    ### v2
    # x, y, z = r_jm[1], r_jm[2], r_jm[3]
    # c1 = (3/2)cis(r)
    # c2 = 1im/r - 1.0/r^2

    # G[1,1] = ( (1 + c2)/r - (1 + 3*c2) * x^2/r^3) # Gxx
    # G[1,2] = ( - (1 + 3*c2) * x*y/r^3) # Gxy
    # G[1,3] = ( - (1 + 3*c2) * x*z/r^3) # Gxz

    # G[2,1] = (-(1 + 3*c2) *x*y/r^3) # Gyx
    # G[2,2] = ((1 + c2)/r - (1 + 3*c2) * y^2/r^3) # Gyy
    # G[2,3] = (-(1 + 3*c2) * y * z/r^3) # Gyz

    # G[3,1] = ( -(1 + 3*c2) * x * z/r^3) # Gzx
    # G[3,2] = ( -(1 + 3*c2) * y * z/r^3) # Gzy
    # G[3,3] = ((1 + c2)/r - (1 + 3*c2) * z^2/r^3) # Gzz

    # G .*= c1

    ### v3
    x, y, z = r_jm[1], r_jm[2], r_jm[3]
    c1 = (3/2)*cis(r)/(r^3)
    c2 = 1im/r - 1/r^2

    G[1,1] = c1*( (1 + c2)*r^2 - (1 + 3*c2) * x^2) # Gxx
    G[1,2] = c1*(-(1 + 3*c2) * x * y) # Gxy
    G[1,3] = c1*(-(1 + 3*c2) * x * z) # Gxz

    G[2,1] = c1*(-(1 + 3*c2) * x * y) # Gyx
    G[2,2] = c1*((1 + c2)*r^2 - (1 + 3*c2) * y^2) # Gyy
    G[2,3] = c1*(-(1 + 3*c2) * y * z) # Gyz

    G[3,1] = c1*(-(1 + 3*c2) * x * z) # Gzx
    G[3,2] = c1*(-(1 + 3*c2) * y * z) # Gzy
    G[3,3] = c1*((1 + c2)*r^2 - (1 + 3*c2) * z^2) # Gzz
    
    
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
    r_far_field = how_far_is_farField(atoms)
    
    norm_sensor = norm(sensor)
    n̂ = sensor / norm_sensor
    
    E_scatt = zeros(Complex{eltype(r)}, 3)

    r_jm = zeros(eltype(atoms.r), 3)
    r_jm[1] = sensor[1]
    r_jm[2] = sensor[2]
    r_jm[3] = sensor[3]

    G = zeros(ComplexF64, 3, 3)
    _vectorial_3D_green_kernel_far_field!(r_jm, G) # use 'r_jm' == 'sensor' and NOT its normalized version
    
    for (j, rⱼ) = enumerate(eachcol(r))
        βⱼ = view(β, :, j)
        E_scatt += (G*βⱼ).*cis(-dot(n̂, rⱼ))
    end

    return +im*(3Γ/4)*(cis(r_far_field)/r_far_field).*E_scatt
end


@inline function _vectorial_3D_green_kernel_far_field!(r_jm::Vector, G::Matrix)
    # this code comes from removing 'c2' from :near_field
    
    # x, y, z = r_jm[1], r_jm[2], r_jm[3]
    # r = norm(r_jm)

    # G[1,1] = r  - x^2 # Gxx
    # G[1,2] = -x * y # Gxy
    # G[1,3] = -x * z # Gxz

    # G[2,1] = -y * x # Gyx
    # G[2,2] = r  - y^2 # Gyy
    # G[2,3] = -y * z # Gyz

    # G[3,1] = -z * x # Gzx
    # G[3,2] = -z * y # Gzy
    # G[3,3] = r  - z^2 # Gzz

    r = k₀ * norm(r_jm)
    c1 = (3/2)*cis(r)/r
    c2 = 1im/r - 1/r^2
    
    n_jm = r_jm./r
    n_x, n_y, n_z = n_jm[1], n_jm[2], n_jm[3]

    G[1,1] = c1*( (1 - n_x*n_x)  ) # Gxx
    G[1,2] = c1*( (0 - n_x*n_y)  ) # Gxy
    G[1,3] = c1*( (0 - n_x*n_z)  ) # Gxz

    G[2,1] = c1*( (0 - n_y*n_x)  ) # Gyx
    G[2,2] = c1*( (1 - n_y*n_y)  ) # Gyy
    G[2,3] = c1*( (0 - n_y*n_z)  ) # Gyz

    G[3,1] = c1*( (0 - n_z*n_x)  ) # Gzx
    G[3,2] = c1*( (0 - n_z*n_y)  ) # Gzy
    G[3,3] = c1*( (1 - n_z*n_z)  ) # Gzz
        
    return nothing
end