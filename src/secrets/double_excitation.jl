using DifferentialEquations, Random, LinearAlgebra, MKL, TensorOperations, ProgressMeter, Plots, LsqFit

@views function doubleExcitation_memoryTuned(du, u, p, t)
    N, G, Ωₗ, Γₗ, Δₗ, dβₖₗ = p
    βₗ = u[1:N]
    βₖₗ = reshape(u[(N + 1):end], N, N)

    for l in 1:N
        du[l] = (im * Δₗ[l] - Γₗ[l] / 2) * βₗ[l] - 0.5 * im * Ωₗ[l]  # v2
    end
    mul!(du[1:N], -G, βₗ) # du[1:N] .= - G*βₗ

    for k in 1:N
        Ωₖ = Ωₗ[k]
        Gₖ = G[k, :]
        βₖ_linha = βₖₗ[k, :]
        βₖ = βₗ[k]
        for l in 1:N
            dβₖₗ[k, l] = (im * (Δₗ[k] + Δₗ[l]) - 0.5 * (Γₗ[k] + Γₗ[l])) * βₖₗ[k, l] - 0.5im * (Ωₗ[l] * βₖ + Ωₗ[k] * βₗ[l]) - dot(G[l, :], βₖ_linha) - dot(Gₖ, βₖₗ[:, l])
        end
    end
    du[(N + 1):end] .= dβₖₗ[:]

    return nothing
end
function doubleExcitation_speedTuned(du, u, p, t)
    N, G, Ωₗ, Γₗ, Δₗ, rightIndex, _temp = p
    βₗ = u[1:N]
    βₖₗ = reshape(u[(N + 1):end], N, N)

    for l in 1:N
        du[l] = (im * Δₗ[l] - Γₗ[l] / 2) * βₗ[l] - 0.5 * im * Ωₗ[l]
    end
    mul!(du[1:N], G, βₗ)

    @tensor begin
        _temp[k, l] = G[l, m] * βₖₗ[k, m] + G[k, m] * βₖₗ[m, l]
    end

    for k in 1:N
        for l in 1:N
            du[N + rightIndex[k, l]] = (im * (Δₗ[k] + Δₗ[l]) - 0.5 * (Γₗ[k] + Γₗ[l])) * βₖₗ[k, l] - 0.5im * (Ωₗ[l] * βₗ[k] + Ωₗ[k] * βₗ[l]) + _temp[k, l]
        end
    end

    return nothing
end

Random.seed!(2022)

N_range = round.(Int, 10.0 .^ range(log10(10), log10(150); length=20))
p = Progress(length(N_range); showspeed=true)
times_memory = Float64[]
times_speed = Float64[]

for N in N_range
    G = rand(ComplexF64, N, N)
    G[diagind(G)] .= 0 # avoid the case "m ≠ l"  for dot product
    G .= -G  # avoid an extra operation of negative ("-") inside the TensorOperations

    Ωₗ = rand(ComplexF64, N)
    Γₗ = rand(ComplexF64, N)
    Δₗ = rand(N)
    rightIndex = LinearIndices(G)
    _temp = similar(G)

    p_memory = N, G, Ωₗ, Γₗ, Δₗ, _temp
    p_speed = N, G, Ωₗ, Γₗ, Δₗ, rightIndex, _temp

    u0 = rand(ComplexF64, N + N^2)

    prob_memory = ODEProblem(doubleExcitation_memoryTuned, u0, (0.0, 5.0), p_memory)
    prob_speed = ODEProblem(doubleExcitation_speedTuned, u0, (0.0, 5.0), p_speed)

    elapsed_memory = @elapsed solve(prob_memory, Tsit5())
    elapsed_speed = @elapsed solve(prob_speed, Tsit5())

    push!(times_memory, elapsed_memory)
    push!(times_speed, elapsed_speed)
    ProgressMeter.next!(p)
end

plot(N_range, times_memory; label="Memory Optimized", scale=:log10, markershape=:circle)
plot!(N_range, times_speed; label="Speed Optimized", legend=:topleft, markershape=:circle)
xlabel!("Number of Atoms")
ylabel!(" Time [seconds] ")
@. model_fit(x, p) = p[1] * x^p[2]

fit_result = curve_fit(model_fit, N_range, times_memory, rand(2))
p1 = fit_result.param[1]
p2 = fit_result.param[2]
fit_memory = p1 * N_range .^ p2
plot!(N_range, fit_memory; label="y = $(round(p1, digits=8)) * x^$(round(p2, digits=3))", c=:black, linestyle=:dash)

fit_result = curve_fit(model_fit, N_range, times_speed, rand(2))
p1 = fit_result.param[1]
p2 = fit_result.param[2]
fit_speed = p1 * N_range .^ p2
plot!(N_range, fit_speed; label="y = $(round(p1, digits=8)) * x^$(round(p2, digits=3))", c=:black, linestyle=:dot, lw=3)
