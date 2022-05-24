using LinearAlgebra, Random
Random.seed!(55732)
const k₀ = Γ = 1

# sensor
r⃗ = 200 * [1, 1, 1] # 200 units from the origin is a far field value for this example
r = norm(r⃗)
n̂ = r⃗ ./ r
Ω = rand(ComplexF64) # campo elétrico em r

# atoms
N = 1000
R⃗ = rand(3, N) # atomic positions
σ⁻ = rand(ComplexF64, N) # atomic excitations
σ⁺ = conj.(σ⁻)
σᶻ = 2σ⁺ .* σ⁻ .- 1

# intensities
E = -(im / 2) * Ω - (Γ / 2) * (exp(im * k₀ * r) / (im * k₀ * r)) * sum(σ⁻[j] * exp(-im * k₀ * (n̂ ⋅ R⃗[:, j])) for j in 1:N)
I_1 = real(E * conj(E))

term1 = abs2(Ω) / 4
term2 = real(im * Ω * (exp(-im * k₀ * r) / (im * k₀ * r)) * sum(σ⁺[j] * exp(+im * k₀ * (n̂ ⋅ R⃗[:, j])) for j in 1:N))
term3 = sum(σ⁻[j] * σ⁺[m] * exp(-im * k₀ * (n̂ ⋅ (R⃗[:, j] - R⃗[:, m]))) for j = 1:N, m = 1:N if j ≠ m)
term4 = sum((1 + σᶻ[j]) / 2 for j in 1:N)
I_2 = real(term1 - (Γ / 2) * term2 + (Γ^2 / 4) * (1 / (k₀ * r)^2) * (term3 + term4))

@show real(I_1), real(I_2), I_1 ≈ I_2;
