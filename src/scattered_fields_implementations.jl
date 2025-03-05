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

    for (j, rⱼ) = enumerate(eachcol(r))
        βⱼ = view(β, :, j)

        Xoj = sensor[1] .- rⱼ[1]
        Yoj = sensor[2] .- rⱼ[2]
        Zoj = sensor[3] .- rⱼ[3]
        Roj = sqrt(Xoj^2 + Yoj^2 + Zoj^2)

        # Vectorial field
        c1 = (3/2)*cis(Roj)./(1im*Roj^3)
        c2 = 1im/Roj - 1/(Roj^2)

        betax = βⱼ[1]
        betay = βⱼ[2]
        betaz = βⱼ[3]

        E_scatt[1] += c1 .* (
                ( Roj.^2 - Xoj.^2 + c2.* (Roj.^2 - 3* Xoj.^2) ) .* betax +
                ( -Xoj.* Yoj.- 3* c2.* Xoj.* Yoj ) .* betay +
                ( -Xoj.* Zoj.- 3* c2.* Xoj.* Zoj ) .* betaz
            )

        E_scatt[2] += c1 .* (
                ( Roj.^2 .- Yoj.^2 .+ c2.* (Roj.^2 .- 3* Yoj.^2) ) .* betay +
                ( -Xoj.* Yoj.- 3* c2.* Xoj.* Yoj ) .* betax +
                ( -Yoj.* Zoj.- 3* c2.* Yoj.* Zoj ) .* betaz
            )


        E_scatt[3] += c1 .* (
                ( Roj.^2 .- Zoj.^2 .+ c2.* (Roj.^2 .- 3* Zoj.^2) ) .* betaz +
                ( -Xoj.* Zoj .- 3* c2.* Xoj.* Zoj ) .* betax +
                ( -Yoj.* Zoj .- 3* c2.* Yoj.* Zoj ) .* betay
            )
    end
    return -Γ*E_scatt # I don't know how to justify that this works    
end
function _vectorial_3D_green_single_atom(r_jm::Vector)
    G = Array{Complex{eltype(r_jm)}}(undef, 3,3)
    _vectorial_3D_green_single_atom(r_jm, G)
    return G
end


@inline function _vectorial_3D_green_single_atom(r_jm::Vector, G::Matrix)
#   reverse engineered (and adapted) from Ana Cipris
    r = k₀ * norm(r_jm)
    r2 = r^2

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
## NEEDS BENCHMARKS 
@views function _vectorial_scattering_far_field(atoms::Atom{T}, β, sensor) where T <: ThreeD
    N = atoms.N
    r = atoms.r[:,:]
    r_far_field = how_far_is_farField(atoms)

    norm_sensor = sqrt(sum(abs2, sensor))
    n̂ = sensor ./ norm_sensor
    δ = Diagonal([1,1,1])
    E_scatt = zeros(Complex{eltype(r)}, 3)

    for ζ = 1:3
        E_scatt[ζ] = sum( (δ[ζ,η]  - n̂[ζ]*n̂[η])*cis(-dot(n̂, r[:,j]))*β[η, j] for j=1:N, η=1:3 )
    end

    return +im*(3Γ/4)*(cis(norm_sensor)/(norm_sensor)).*E_scatt    
end