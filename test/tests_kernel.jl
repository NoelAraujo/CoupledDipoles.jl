using LinearAlgebra
function kernel3c(Xt, Yt, Zt)
    N = length(Xt)
    k = 2 * π / 461e-9
    Xjm = Xt * ones(1, N) - ones(N, 1) * Xt'
    Yjm = Yt * ones(1, N) - ones(N, 1) * Yt'
    Zjm = Zt * ones(1, N) - ones(N, 1) * Zt'
    Rjm = sqrt.(Xjm.^2 + Yjm.^2 + Zjm.^2) .+ diagm(ones(N))
    
    Ker = zeros(ComplexF64, 3N, 3N)
    c1 = exp.(1im * k * Rjm) ./ (1im * k * Rjm.^3)
    c2 = 1im ./ (k * Rjm) - 1.0 ./ (k * Rjm).^2
    
    
    Ker[1:N, 1:N] .= c1 .* (Rjm.^2 - Xjm.^2 + c2 .* (Rjm.^2 - 3 * Xjm.^2))
    Ker[N+1:2N, N+1:2N] .= c1 .* (Rjm.^2 - Yjm.^2 + c2 .* (Rjm.^2 - 3 * Yjm.^2))
    Ker[2N+1:3N, 2N+1:3N] .= c1 .* (Rjm.^2 - Zjm.^2 + c2 .* (Rjm.^2 - 3 * Zjm.^2))
    Ker[1:N, N+1:2N] .= c1 .* (-Xjm .* Yjm - 3 * c2 .* Xjm .* Yjm)
    Ker[N+1:2N, 1:N] .= Ker[1:N, N+1:2N]
    Ker[1:N, 2N+1:3N] .= c1 .* (-Xjm .* Zjm - 3 * c2 .* Xjm .* Zjm)
    Ker[2N+1:3N, 1:N] .= Ker[1:N, 2N+1:3N]
    Ker[N+1:2N, 2N+1:3N] .= c1 .* (-Yjm .* Zjm - 3 * c2 .* Yjm .* Zjm)
    Ker[2N+1:3N, N+1:2N] .= Ker[N+1:2N, 2N+1:3N]
    
    Ker = Ker .- diagm(diag(Ker))  # Removing artificially non-zero diagonal terms
    kn = Ker
end



Xt = [1, 2, 3]
Yt = [4, 5, 6]
Zt = [7, 8, 9]

kn_julia = kernel3c(Xt, Yt, Zt)


# matlab
Xt = [1; 2; 3];
Yt = [4; 5; 6];
Zt = [7; 8; 9];




# Xt, Yt, Zt = atoms.r[1, :], atoms.r[2, :], atoms.r[3, :]
# N = length(Xt)
# k = k₀

# Xjm = Xt * ones(1, N) - ones(N, 1) * Xt'
# Yjm = Yt * ones(1, N) - ones(N, 1) * Yt'
# Zjm = Zt * ones(1, N) - ones(N, 1) * Zt'
# Rjm = sqrt.(Xjm.^2 + Yjm.^2 + Zjm.^2) .+ diagm(ones(N))

# c1 = exp.(1im * k * Rjm) ./ (1im * k * Rjm.^3)
# c2 = 1im ./ (k * Rjm) - 1.0 ./ (k * Rjm).^2


# G[1:N, 1:N] .= c1 .* (Rjm.^2 - Xjm.^2 + c2 .* (Rjm.^2 - 3 * Xjm.^2))
# G[N+1:2N, N+1:2N] .= c1 .* (Rjm.^2 - Yjm.^2 + c2 .* (Rjm.^2 - 3 * Yjm.^2))
# G[2N+1:3N, 2N+1:3N] .= c1 .* (Rjm.^2 - Zjm.^2 + c2 .* (Rjm.^2 - 3 * Zjm.^2))
# G[1:N, N+1:2N] .= c1 .* (-Xjm .* Yjm - 3 * c2 .* Xjm .* Yjm)
# G[N+1:2N, 1:N] .= G[1:N, N+1:2N]
# G[1:N, 2N+1:3N] .= c1 .* (-Xjm .* Zjm - 3 * c2 .* Xjm .* Zjm)
# G[2N+1:3N, 1:N] .= G[1:N, 2N+1:3N]
# G[N+1:2N, 2N+1:3N] .= c1 .* (-Yjm .* Zjm - 3 * c2 .* Yjm .* Zjm)
# G[2N+1:3N, N+1:2N] .= G[N+1:2N, 2N+1:3N]