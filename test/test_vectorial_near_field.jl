function Escattered3(Xi, Yi, Zi, Beta3, Xobs, Yobs, Zobs, k0)
    N = length(Xi)
    Escat3x = zeros(ComplexF64, size(Xobs))
    Escat3y = zeros(ComplexF64, size(Yobs))
    Escat3z = zeros(ComplexF64, size(Zobs))

    for ii = 1:N
        Xoj = Xobs .- Xi[ii]
        Yoj = Yobs .- Yi[ii]
        Zoj = Zobs .- Zi[ii]
        Roj = sqrt.(Xoj.^2 .+ Yoj.^2 .+ Zoj.^2)

        # Vectorial field
        c1 = 3*exp.(1im*k0*Roj)./(2*1im*k0*Roj.^3)
        c2 = 1im./(k0*Roj) .- 1.0./(k0*Roj).^2

        Escat3x .-= c1 .* (
                ((1 .+ c2) .* Roj.^2 .- (1 .+ 3*c2) .* Xoj.^2) .* Beta3[ii] +
                (-(1 .+ 3*c2) .* Xoj .* Yoj) .* Beta3[ii + N] +
                (-(1 .+ 3*c2) .* Xoj .* Zoj) .* Beta3[ii + 2*N]
            )

        Escat3y .-= c1 .* (
                ((1 .+ c2) .* Roj.^2 .- (1 .+ 3*c2) .* Yoj.^2) .* Beta3[ii + N] +
                (-(1 .+ 3*c2) .* Xoj .* Yoj) .* Beta3[ii] +
                (-(1 .+ 3*c2) .* Yoj .* Zoj) .* Beta3[ii + 2*N]
            )

        Escat3z .-= c1 .* (
                ((1 .+ c2) .* Roj.^2 .- (1 .+ 3*c2) .* Zoj.^2) .* Beta3[ii + 2*N] +
                (-(1 .+ 3*c2) .* Xoj .* Zoj) .* Beta3[ii] +
                (-(1 .+ 3*c2) .* Yoj .* Zoj) .* Beta3[ii + N]
            )
    end

    return Escat3x, Escat3y, Escat3z
end

function Escattered3_v2(Xi, Yi, Zi, Beta3, Xobs, Yobs, Zobs, k0)
    N = length(Xi)
    Escat3x = zeros(ComplexF64, size(Xobs))
    Escat3y = zeros(ComplexF64, size(Yobs))
    Escat3z = zeros(ComplexF64, size(Zobs))

    for ii = 1:N
        Xoj = Xobs .- Xi[ii]
        Yoj = Yobs .- Yi[ii]
        Zoj = Zobs .- Zi[ii]
        Roj = sqrt.(Xoj.^2 .+ Yoj.^2 .+ Zoj.^2)

        # Vectorial field
        c1 = 3*exp.(1im*k0*Roj)./(2*1im*k0*Roj.^3)
        c2 = 1im./(k0*Roj) .- 1.0./(k0*Roj).^2
        Gxx = ( (1 .+ c2) .* Roj.^2 .- (1 .+ 3*c2) .* Xoj.^2)
        Gxy = (-(1 .+ 3*c2) .* Xoj .* Yoj)
        Gxz = (-(1 .+ 3*c2) .* Xoj .* Zoj)
        Escat3x .-= c1 .* (
                Gxx .* Beta3[ii] +
                Gxy .* Beta3[ii + N] +
                Gxz .* Beta3[ii + 2*N]
            )

        Gyy = ((1 .+ c2) .* Roj.^2 .- (1 .+ 3*c2) .* Yoj.^2)
        Gyx = (-(1 .+ 3*c2) .* Xoj .* Yoj)
        Gyz = (-(1 .+ 3*c2) .* Yoj .* Zoj)
        Escat3y .-= c1 .* (
                Gyy.* Beta3[ii + N] +
                Gyx.* Beta3[ii] +
                Gyz.* Beta3[ii + 2*N]
            )

        Gzz = ((1 .+ c2) .* Roj.^2 .- (1 .+ 3*c2) .* Zoj.^2)
        Gyx = (-(1 .+ 3*c2) .* Xoj .* Zoj)
        Gzy = (-(1 .+ 3*c2) .* Yoj .* Zoj)
        Escat3z .-= c1 .* (
                Gzz .* Beta3[ii + 2*N] +
                Gyx .* Beta3[ii] +
                Gzy.* Beta3[ii + N]
            )
    end

    return Escat3x, Escat3y, Escat3z
end


function Escattered3_v3(Xi, Yi, Zi, Beta3, Xobs, Yobs, Zobs, k0; idx_sensor=1)
    N = length(Xi)
    Escat3x = zero(ComplexF64)
    Escat3y = zero(ComplexF64)
    Escat3z = zero(ComplexF64)

    for ii = 1:N
        Xoj = Xobs[idx_sensor] .- Xi[ii]
        Yoj = Yobs[idx_sensor] .- Yi[ii]
        Zoj = Zobs[idx_sensor] .- Zi[ii]
        Roj = sqrt.(Xoj.^2 .+ Yoj.^2 .+ Zoj.^2)

        # Vectorial field
        c1 = 3*exp.(1im*k0*Roj)./(2*1im*k0*Roj.^3)
        c2 = 1im./(k0*Roj) .- 1.0./(k0*Roj).^2
        
        betax = Beta3[ii]
        betay = Beta3[ii + N]
        betaz = Beta3[ii + 2*N]

        Gxx = ( (1 .+ c2) .* Roj.^2 .- (1 .+ 3*c2) .* Xoj.^2)
        Gxy = (-(1 .+ 3*c2) .* Xoj .* Yoj)
        Gxz = (-(1 .+ 3*c2) .* Xoj .* Zoj)
        Escat3x -= c1 .* (
                Gxx * betax +
                Gxy * betay +
                Gxz * betaz
            )

        Gyy = ((1 .+ c2) .* Roj.^2 .- (1 .+ 3*c2) .* Yoj.^2)
        Gyx = (-(1 .+ 3*c2) .* Xoj .* Yoj)
        Gyz = (-(1 .+ 3*c2) .* Yoj .* Zoj)
        Escat3y -= c1 .* (
                Gyx * betax +
                Gyy * betay +
                Gyz * betaz
            )

        Gzz = ((1 .+ c2) .* Roj.^2 .- (1 .+ 3*c2) .* Zoj.^2)
        Gyx = (-(1 .+ 3*c2) .* Xoj .* Zoj)
        Gzy = (-(1 .+ 3*c2) .* Yoj .* Zoj)
        Escat3z -= c1 .* (
                Gyx * betax +
                Gzy *  betay +
                Gzz * betaz
            )
    end

    return Escat3x, Escat3y, Escat3z
end



@views function Escattered3_v4(atoms, Beta, sensor; k0=1, idx_sensor=1)
    Xi, Yi, Zi = atoms[1, :], atoms[2, :], atoms[3, :]
    N = length(Xi)
    
    Xobs = sensor[1, idx_sensor]
    Yobs = sensor[2, idx_sensor]
    Zobs = sensor[3, idx_sensor]

    Escat3x = zero(ComplexF64)
    Escat3y = zero(ComplexF64)
    Escat3z = zero(ComplexF64)
    for j = 1:N
		betax, betay, betaz = Beta[1,j], Beta[2,j], Beta[3,j]
         	
        Xoj = Xobs .- Xi[j]
        Yoj = Yobs .- Yi[j]
        Zoj = Zobs .- Zi[j]
        Roj = sqrt.(Xoj.^2 .+ Yoj.^2 .+ Zoj.^2)

        # Vectorial field
        c1 = (3/2)cis(k0*Roj)/(1im*k0*Roj^3)
        c2 = 1im/(k0*Roj) - 1.0/(k0*Roj)^2

        Gxx = ( (1 + c2) * Roj.^2 - (1 + 3*c2) * Xoj.^2)
        Gxy = (-(1 + 3*c2) * Xoj * Yoj)
        Gxz = (-(1 + 3*c2) * Xoj * Zoj)

		Gyx = (-(1 + 3*c2) * Xoj * Yoj)
		Gyy = ((1 + c2) * Roj.^2 - (1 + 3*c2) * Yoj^2)
        Gyz = (-(1 + 3*c2) * Yoj * Zoj)

		Gzx = (-(1 + 3*c2) * Xoj * Zoj)
        Gzy = (-(1 + 3*c2) * Yoj * Zoj)
		Gzz = ((1 + c2) * Roj.^2 - (1 + 3*c2) * Zoj^2)
		
		Escat3x -= c1 * ( Gxx * betax + Gxy * betay + Gxz * betaz)
        Escat3y -= c1 * ( Gyx * betax + Gyy * betay + Gyz * betaz)
        Escat3z -= c1 * ( Gzx * betax + Gzy * betay + Gzz * betaz)
    end

    return Escat3x, Escat3y, Escat3z
end


@views function Escattered3_v4(atoms, Beta, sensor; k0=1, idx_sensor=1)
    Xi, Yi, Zi = atoms[1, :], atoms[2, :], atoms[3, :]
    N = length(Xi)
    
    Xobs = sensor[1, idx_sensor]
    Yobs = sensor[2, idx_sensor]
    Zobs = sensor[3, idx_sensor]

    Escat3x = zero(ComplexF64)
    Escat3y = zero(ComplexF64)
    Escat3z = zero(ComplexF64)
	G = zeros(ComplexF64, 3, 3)
    for j = 1:N
		betax, betay, betaz = Beta[1,j], Beta[2,j], Beta[3,j]
         	
        Xoj = Xobs - Xi[j]
        Yoj = Yobs - Yi[j]
        Zoj = Zobs - Zi[j]
        Roj = sqrt(Xoj^2 + Yoj^2 + Zoj^2)

        # Vectorial field
        c1 = (3/2)cis(k0*Roj)/(1im*k0*Roj^3)
        c2 = 1im/(k0*Roj) - 1.0/(k0*Roj)^2

        G[1,1] = ( (1 + c2) * Roj.^2 - (1 + 3*c2) * Xoj.^2) # Gxx
        G[1,2] = (-(1 + 3*c2) * Xoj * Yoj) # Gxy
        G[1,3] = (-(1 + 3*c2) * Xoj * Zoj) # Gxz

		G[2,1] = (-(1 + 3*c2) * Xoj * Yoj) # Gyx
		G[2,2] = ((1 + c2) * Roj.^2 - (1 + 3*c2) * Yoj^2) # Gyy
        G[2,3] = (-(1 + 3*c2) * Yoj * Zoj) # Gyz

		G[3,1] = (-(1 + 3*c2) * Xoj * Zoj) # Gzx
        G[3,2] = (-(1 + 3*c2) * Yoj * Zoj) # Gzy
		G[3,3] = ((1 + c2) * Roj.^2 - (1 + 3*c2) * Zoj^2) # Gzz

		_Escat3 = c1.*(G*Beta[:,j])

		Escat3x -= _Escat3[1]
        Escat3y -= _Escat3[2]
        Escat3z -= _Escat3[3]
    end

    return Escat3x, Escat3y, Escat3z
end

Xi = [1, 2, 3, 11];
Yi = [4, 5, 6, 3];
Zi = [7, 8, 9, -.5];
Beta3 = [10, 11, 12, 13, 14, 15, 16, 17, 18, 2, 10, 3];
Xobs = [19, 20, 21, 30];
Yobs = [22, 23, 24, -32];
Zobs = [25, 26, 27, 24];
k0 = 1;

@time Escat3x, Escat3y, Escat3z = Escattered3(Xi, Yi, Zi, Beta3, Xobs, Yobs, Zobs, k0);
@time Escat3x_v2, Escat3y_v2, Escat3z_v2 = Escattered3_v2(Xi, Yi, Zi, Beta3, Xobs, Yobs, Zobs, k0);
@time Escat3x_v3, Escat3y_v3, Escat3z_v3 = Escattered3_v3(Xi, Yi, Zi, Beta3, Xobs, Yobs, Zobs, k0);

Escat3x[1] ≈ Escat3x_v3
Escat3y[1] ≈ Escat3y_v3


Escat3x ≈ Escat3x_v2
Escat3y ≈ Escat3y_v2
Escat3z ≈ Escat3z_v2


