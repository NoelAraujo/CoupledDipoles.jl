using LinearAlgebra
function kernel3c!(Xt, Yt, Zt, Ker, c1, c2, k)
    # Calculate the number of particles
    N = length(Xt)

    # Calculate the relative positions between particles
    Xjm = -(Xt .* ones(1, N) .- ones(N, 1) .* Xt')
    Yjm = -(Yt .* ones(1, N) .- ones(N, 1) .* Yt')
    Zjm = -(Zt .* ones(1, N) .- ones(N, 1) .* Zt')

    # Calculate the distances between particles
    Rjm = sqrt.(Xjm.^2 .+ Yjm.^2 .+ Zjm.^2) .+ diagm(ones(N))

    # Calculate the coefficients
    for i in eachindex(c1)
        c1[i] = exp(im * k * Rjm[i]) / (im * k * Rjm[i]^3)
    end
    c2 .= im ./ (k * Rjm) .- 1 ./ (k * Rjm).^2

    # Calculate the kernel matrix elements
    Ker_xx = c1 .* (Rjm.^2 .- Xjm.^2 .+ c2 .* (Rjm.^2 .- 3 .* Xjm.^2))
    Ker_xy = c1 .* (-Xjm .* Yjm .- 3 .* c2 .* Xjm .* Yjm)
    Ker_xz = c1 .* (-Xjm .* Zjm .- 3 .* c2 .* Xjm .* Zjm)
    Ker_yy = c1 .* (Rjm.^2 .- Yjm.^2 .+ c2 .* (Rjm.^2 .- 3 .* Yjm.^2))
    Ker_yz = c1 .* (-Yjm .* Zjm .- 3 .* c2 .* Yjm .* Zjm)
    Ker_zz = c1 .* (Rjm.^2 .- Zjm.^2 .+ c2 .* (Rjm.^2 .- 3 .* Zjm.^2))

    # Assemble the kernel matrix
    Ker[1:N, 1:N] .= Ker_xx
    Ker[N+1:2N, N+1:2N] .= Ker_yy
    Ker[2N+1:3N, 2N+1:3N] .= Ker_zz
    Ker[1:N, N+1:2N] .= Ker_xy
    Ker[N+1:2N, 1:N] .= Ker_xy
    Ker[1:N, 2N+1:3N] .= Ker_xz
    Ker[2N+1:3N, 1:N] .= Ker_xz
    Ker[N+1:2N, 2N+1:3N] .= Ker_yz
    Ker[2N+1:3N, N+1:2N] .= Ker_yz

    # Remove artificially non-zero diagonal terms
    Ker[diagind(Ker)] .= 0im
        
    return Ker
end
