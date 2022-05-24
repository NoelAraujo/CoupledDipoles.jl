function test1!(du, u, p, t)
    G, Ωₙ, N, c = p

    βₙ = @view u[1:N]
    zₙ = @view u[(N + 1):end]

    G_βₙ = G * βₙ
    Wₙ = Ωₙ / 2 .+ im * G_βₙ
    temp1 = c * βₙ .+ im * Wₙ .* zₙ
    temp2 = -Γ * (1 .+ zₙ) - 4 * imag.(βₙ .* conj.(Wₙ))
    du[:] .= vcat(temp1, temp2)

    return nothing
end
function test2!(du, u, p, t)
    G, Ωₙ, N, c = p

    βₙ = @view u[1:N]
    zₙ = @view u[(N + 1):end]

    for j in 1:N
        du[j] = c * βₙ[j] + 0.5 * im * Ωₙ[j] * zₙ[j] - sum(G[j, m] * βₙ[m] for m = 1:N if j ≠ m) * zₙ[j]
    end
    for j in 1:N
        du[N + j] =
            (-im * Ωₙ[j] * conj(βₙ[j]) + im * conj(Ωₙ[j]) * βₙ[j]) - Γ * (1 + zₙ[j]) -
            (-2 / Γ) * (sum(G[j, m] * βₙ[m] * conj(βₙ[j]) for m = 1:N if m ≠ j) + conj.(sum(G[j, m] * βₙ[m] * conj(βₙ[j]) for m = 1:N if m ≠ j)))
    end

    return nothing
end

using LinearAlgebra
const Γ = 1
@time using DifferentialEquations

N = 5
rjk = rand(N, N)
G = -(Γ / 2) * (cos.(rjk) ./ rjk + im * (sin.(rjk) ./ rjk))
G[diagind(G)] .= 0

Ωₙ = rand(ComplexF64, N)
Δ = 2;
c = (im * Δ - Γ / 2);

p = G, Ωₙ, N, c

# TESTING
for _ in 1:5
    u0 = rand(ComplexF64, 2N)
    du1 = similar(u0)
    du2 = similar(u0)
    test1!(du1, u0, p, 0)
    test2!(du2, u0, p, 0)
    # [du1 du2 du1.≈du2] |> display # PRODUCE SAME OUTPUT
    display(all(du1 .≈ du2))
end

# INTEGRATING IN TIME
u0 = rand(ComplexF64, 2N)
prob1 = ODEProblem(test1!, u0, (0.0, 1.0), p)
prob2 = ODEProblem(test2!, u0, (0.0, 1.0), p)

sol1 = solve(prob1; saveat=range(0; stop=1, length=20))
sol2 = solve(prob2; saveat=range(0; stop=1, length=20))

[sol1.u[end] sol2.u[end] sol1.u[end] .≈ sol2.u[end]]
