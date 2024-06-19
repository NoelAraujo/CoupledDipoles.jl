function Escattered3(Xi, Yi, Zi, Beta3, Xobs, Yobs, Zobs, k0)
    N = length(Xi)
    Escat3x = zeros(ComplexF64, size(Xobs))
    Escat3y = zeros(ComplexF64, size(Xobs))
    Escat3z = zeros(ComplexF64, size(Xobs))

    for jN in 1:N
        Xoj = Xobs.- Xi[jN]
        Yoj = Yobs.- Yi[jN]
        Zoj = Zobs.- Zi[jN]
        Roj = sqrt.(Xoj.^2 .+ Yoj.^2 .+ Zoj.^2)

        # Vectorial field
        c1 = exp.(im * k0 * Roj)./ (im * k0.* Roj.^3)
        c2 = im./ (k0 * Roj) .- 1.0./ (k0 * Roj).^2

        # Calculate Escat3x, Escat3y, and Escat3z
        Escat3x .= Escat3x .- c1.* (Roj.^2 .- Xoj.^2 .+ c2.* (Roj.^2 .- 3* Xoj.^2)).* Beta3[jN] .-
                         c1.* (-Xoj.* Yoj .- 3* c2.* Xoj.* Yoj).* Beta3[jN + N] .-
                         c1.* (-Xoj.* Zoj .- 3* c2.* Xoj.* Zoj).* Beta3[jN + 2*N]

        Escat3y .= Escat3y .- c1.* (Roj.^2 .- Yoj.^2 .+ c2.* (Roj.^2 .- 3* Yoj.^2)).* Beta3[jN + N] .-
                         c1.* (-Xoj.* Yoj.- 3* c2.* Xoj.* Yoj).* Beta3[jN] .-
                         c1.* (-Yoj.* Zoj.- 3* c2.* Yoj.* Zoj).* Beta3[jN + 2*N]

        Escat3z .= Escat3z .- c1.* (Roj.^2 .- Zoj.^2 .+ c2.* (Roj.^2 .- 3* Zoj.^2)).* Beta3[jN + 2*N] .-
                         c1.* (-Xoj.* Zoj .- 3* c2.* Xoj.* Zoj).* Beta3[jN] .-
                         c1.* (-Yoj.* Zoj .- 3* c2.* Yoj.* Zoj).* Beta3[jN + N]
    end

    return Escat3x, Escat3y, Escat3z
end
