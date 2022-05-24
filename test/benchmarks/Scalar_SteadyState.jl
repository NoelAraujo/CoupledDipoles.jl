using LinearAlgebra, Random
Random.seed!(55732)
const k₀ = Γ = 1
Δ = 2

# atoms
r = hcat([2; 1; 0], [-8; -1; 4])

# laser over atoms
Ω = rand(ComplexF64, 2)

# dipole interaction
G = zeros(ComplexF64, 2, 2)
for j in 1:2, m in 1:2
    r_jm = norm(r[:, j] - r[:, m])
    G[j, m] = -(Γ / 2) * cis(k₀ * r_jm) / (1im * k₀ * r_jm)
end
G[diagind(G)] .= 1im * Δ - Γ / 2

# steady state
βₛₛ = G \ Ω
