function mean_field_near_field(R, r, β, Ω; Γ=1, k₀=1)
    N = length(β)
    E_L = -im*(1/2).*Ω
    E_S = +im*(Γ/2)*sum( β[j]*exp(im*k₀*norm(R - r[:, j]))/(k₀*norm(R - r[:, j])) for j=1:N)
    E_scatt_1 = E_L + E_S
    I_scatt_1 = abs2(E_scatt_1)

    # term1 = abs2(E_L)
    # term2 = E_L*conj(E_S) + conj(E_L)*E_S
    # term2 = 2*real(conj(E_L)*E_S)
    # term3 = abs2(E_S)

    term1 = abs2(Ω) / 4
    # term2 = (-(Γ/4)*Ω)*conj(sum( β[j]*exp(im*k₀*norm(R - r[:, j]))/(k₀*norm(R - r[:, j])) for j=1:N)) - (Γ/4)*conj(Ω)*(sum( β[j]*exp(im*k₀*norm(R - r[:, j]))/(k₀*norm(R - r[:, j])) for j=1:N))
    term2 = -2real((Γ/4)*conj(Ω)*(sum( β[j]*exp(im*k₀*norm(R - r[:, j]))/(k₀*norm(R - r[:, j])) for j=1:N)))
    term3 = (Γ/2)^2*sum(β[j]*conj(β[m])*exp(im*k₀*norm(R-r[:,j]))*exp(-im*k₀*norm(R-r[:,m]))/( norm(R-r[:,j])*norm(R-r[:,m]) ) for j=1:N, m=1:N)

    I_scatt_2 = real(term1 + term2 + term3)

    return I_scatt_1 ≈ I_scatt_2
end

using LinearAlgebra
using Random; Random.seed!(444)
N = 10
R = rand(3)
r = rand(3, N)
β = rand(ComplexF64, N)
Ω = rand(ComplexF64)

mean_field_near_field(R, r, β, Ω; Γ=1, k₀=1)