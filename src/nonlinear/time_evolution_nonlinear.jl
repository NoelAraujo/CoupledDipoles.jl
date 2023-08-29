# --------------------------------- GENERAL FUNCTIONS ---------------------------------
function time_evolution(problem::NonLinearOptics{T}, u₀, tspan::Tuple; ode_solver=true, change_solver=false, kargs...) where {T<:NonLinear}
    if ode_solver == false
        @warn "NonLinearOptics does not have formal solution. Using numerical solution instead." maxlog = 1
    end

    G = copy(interaction_matrix(problem))
    Ωₙ = laser_field(problem.laser, problem.atoms)

    solution = time_evolution_ode_solver(problem, u₀, tspan, Ωₙ, G; change_solver=change_solver, kargs...)

    return solution
end
function time_evolution_ode_solver(problem::NonLinearOptics{T}, u₀, tspan::Tuple, Ωₙ, G::AbstractMatrix; change_solver=false, kargs...) where {T<:NonLinear}
    parameters = get_evolution_params(problem, G, Ωₙ)

    ### calls for solver
    problemFunction = get_evolution_function(problem)
    prob = ODEProblem(problemFunction, u₀, tspan, parameters)

    if ρ_of(problem.atoms) < 0.1 || change_solver == true
        # RDPK3SpFSAL35
        solution = OrdinaryDiffEq.solve(prob, VCABM3(); abstol=1e-10, kargs...)
    else
        solution = OrdinaryDiffEq.solve(prob, RDPK3Sp35(); abstol=1e-10, kargs...)
    end

    return solution
end



# --------------------------------- MODEL SPECIFIC ---------------------------------

## MeanField
function get_evolution_params(problem::NonLinearOptics{MeanField}, G, Ωₙ)
    #=
        I don't sum over diagonal elements during time evolution
        thus, to avoid an IF statement, I put a zero on diagonal
    =#
    G[diagind(G)] .= zero(eltype(G))

    # `laser_field = (-im/2)*Ω`, but I need only `Ω`
    Ωₙ_clean = vec(Ωₙ) ./ LASER_FACTOR
    Wₙ = similar(Ωₙ_clean)
    G_βₙ = similar(Ωₙ_clean)
    temp1 = similar(Ωₙ_clean)
    temp2 = similar(Ωₙ_clean)

    parameters = view(G, :, :), view(Ωₙ_clean, :), Wₙ, problem.laser.Δ, problem.atoms.N, G_βₙ, temp1, temp2
    return parameters
end

get_evolution_function(problem::NonLinearOptics{MeanField}) = MeanField!
function MeanField!(du, u, p, t)
    # parameters
    G, Ωₙ, Wₙ, Δ, N, G_βₙ, temp1, temp2 = p

    βₙ = @view u[1:N]
    zₙ = @view u[(N+1):end]

    mul!(G_βₙ, G, βₙ) # == G_βₙ = G*βₙ
    # Wₙ .= Ωₙ / 2 .+ im * (G_βₙ - diagG .* βₙ) # if I don't clean the diagonal in callee function
    @. Wₙ = Ωₙ / 2 + im * G_βₙ
    @. temp1 = (im * Δ - Γ / 2) * βₙ + im * Wₙ * zₙ
    @. temp2 = -Γ * (1 + zₙ) - 4 * imag(βₙ * conj(Wₙ))
    du[1:N] .= temp1
    du[N+1:end] .= temp2

    return nothing
end


## PariCorrelation
function get_evolution_params(problem::NonLinearOptics{PairCorrelation}, G, Ωₙ)
    # for the PairCorrelation, the interaction matrix does not contain 'im*(Γ/2)'
    G_c = (2 / Γ) .* G
    G_c .= G_c ./ im


    N = problem.atoms.N
    Gconj = conj.(G_c)
    Γⱼₘ = real.(G_c)

    # `laser_field = (-im/2)*Ω`, but I need only `Ω`
    Ω⁻ = vec(Ωₙ) ./ LASER_FACTOR
    Ω⁺ = conj.(Ω⁻)

    Δ = problem.laser.Δ

    parameters = (N, view(G_c, :, :), view(Gconj, :, :), view(Γⱼₘ, :, :), view(Ω⁻, :), view(Ω⁺, :), Δ)

    return parameters
end

get_evolution_function(problem::NonLinearOptics{PairCorrelation}) = PairCorrelation!
function PairCorrelation!(du, u, params, t)
    N = params[1]
    G = params[2]
    Gconj = params[3]
    Γⱼₘ = params[4]

    Ω⁻ = params[5]
    Ω⁺ = params[6]
    Δ = params[7]

    σ⁻ = u[1:N]
    σᶻ = u[N+1:2*N]
    σ⁺ = conj.(σ⁻)

    σᶻσ⁻ = reshape(u[2*N+1:2*N+N^2], (N, N))
    σ⁺σ⁻ = reshape(u[2*N+1+N^2:2*N+2*N^2], (N, N))
    σ⁻σ⁻ = reshape(u[2*N+1+2*N^2:2*N+3*N^2], (N, N))
    σᶻσᶻ = reshape(u[2*N+1+3*N^2:2*N+4*N^2], (N, N))

    σ⁻σᶻ = transpose(σᶻσ⁻)
    σ⁺σᶻ = transpose(conj.(σᶻσ⁻))
    σᶻσ⁺ = transpose(σ⁺σᶻ)

    σ⁺σ⁻σ⁻ = truncation(σ⁺, σ⁻, σ⁻, σ⁺σ⁻, σ⁺σ⁻, σ⁻σ⁻)
    σᶻσᶻσ⁻ = truncation(σᶻ, σᶻ, σ⁻, σᶻσᶻ, σᶻσ⁻, σᶻσ⁻)
    σ⁺σ⁻σᶻ = truncation(σ⁺, σ⁻, σᶻ, σ⁺σ⁻, σ⁺σᶻ, σ⁻σᶻ)
    σᶻσ⁻σ⁻ = truncation(σᶻ, σ⁻, σ⁻, σᶻσ⁻, σᶻσ⁻, σ⁻σ⁻)
    σ⁺σᶻσ⁻ = truncation(σ⁺, σᶻ, σ⁻, σ⁺σᶻ, σ⁺σ⁻, σᶻσ⁻)
    σᶻσ⁺σ⁻ = truncation(σᶻ, σ⁺, σ⁻, σᶻσ⁺, σᶻσ⁻, σ⁺σ⁻)

    dₜ_σ⁻_p1 = @tullio p1[j] := begin
        (im * Δ - Γ / 2) * σ⁻[j] + im * (Ω⁻[j] / 2) * σᶻ[j]
    end
    dₜ_σ⁻_p2 = @tullio p2[j] := begin
        if m ≠ j
            (Γ / 2) * G[j, m] * σᶻσ⁻[j, m]
        else
            zero(eltype(u))
        end
    end
    dₜ_σ⁻ = dₜ_σ⁻_p1 + dₜ_σ⁻_p2 # <-----

    dₜ_σᶻ_p1 = @tullio p3[j] := begin
        im * (Ω⁺[j] * σ⁻[j] - conj(Ω⁺[j] * σ⁻[j])) - Γ * (1 + σᶻ[j])
    end
    dₜ_σᶻ_p2 = @tullio p4[j] := begin
        if m ≠ j
            -Γ * (G[j, m] * σ⁺σ⁻[j, m] + conj(G[j, m] * σ⁺σ⁻[j, m]))
        else
            zero(eltype(u))
        end
    end
    dₜ_σᶻ = dₜ_σᶻ_p1 + dₜ_σᶻ_p2 # <-----


    dₜ_σᶻσ⁻_p1 = @tullio p6[j, m] := begin
        if j ≠ m
            (im * Δ - 3Γ / 2) * σᶻσ⁻[j, m] - Γ * σ⁻[m] + im * (Ω⁺[j] * σ⁻σ⁻[j, m] - Ω⁻[j] * σ⁺σ⁻[j, m] + 0.5 * Ω⁻[m] * σᶻσᶻ[j, m]) - Γ * Γⱼₘ[j, m] * σ⁻σᶻ[j, m] - (Γ / 2) * Gconj[j, m] * σ⁻[j]
        else
            zero(eltype(u))
        end
    end
    dₜ_σᶻσ⁻_p2 = @tullio p7[j, m] := begin
        if (k ≠ j) && (k ≠ m) && (j ≠ m)
            -Γ * (G[j, k] * σ⁺σ⁻σ⁻[j, m, k] + Gconj[j, k] * σ⁺σ⁻σ⁻[k, m, j]) + 0.5Γ * (G[m, k] * σᶻσᶻσ⁻[m, j, k])
        else
            zero(eltype(u))
        end
    end
    dₜ_σᶻσ⁻ = dₜ_σᶻσ⁻_p1 + dₜ_σᶻσ⁻_p2 # <-----


    dₜ_σ⁺σ⁻_p1 = @tullio p8[j, m] := begin
        if j ≠ m
            -Γ * σ⁺σ⁻[j, m] - 0.5im * (Ω⁺[j] * σᶻσ⁻[j, m] - Ω⁻[m] * σ⁺σᶻ[j, m]) + (Γ / 4) * (G[j, m] * σᶻ[m] + Gconj[j, m] * σᶻ[j]) + (Γ / 2) * Γⱼₘ[j, m] * σᶻσᶻ[j, m]
        else
            zero(eltype(u))
        end
    end
    dₜ_σ⁺σ⁻_p2 = @tullio p9[j, m] := begin
        if (k ≠ j) && (k ≠ m) && (j ≠ m)
            +(Γ / 2) * (Gconj[j, k] * σ⁺σ⁻σᶻ[k, m, j] + G[m, k] * σ⁺σ⁻σᶻ[j, k, m])
        else
            zero(ComplexF64)
        end
    end
    dₜ_σ⁺σ⁻ = dₜ_σ⁺σ⁻_p1 + dₜ_σ⁺σ⁻_p2 # <-----

    dₜ_σ⁻σ⁻_p1 = @tullio p10[j, m] := begin
        if j ≠ m
            (2im * Δ - Γ) * σ⁻σ⁻[j, m] + 0.5im * (Ω⁻[j] * σᶻσ⁻[j, m] + Ω⁻[m] * σᶻσ⁻[m, j])
        else
            zero(ComplexF64)
        end
    end
    dₜ_σ⁻σ⁻_p2 = @tullio p11[j, m] := begin
        if (k ≠ j) && (k ≠ m) && (j ≠ m)
            +(Γ / 2) * (G[j, k] * σᶻσ⁻σ⁻[j, m, k] + G[m, k] * σᶻσ⁻σ⁻[m, j, k])
        else
            zero(ComplexF64)
        end
    end
    dₜ_σ⁻σ⁻ = dₜ_σ⁻σ⁻_p1 + dₜ_σ⁻σ⁻_p2 # <-----

    dₜ_σᶻσᶻ_p1 = @tullio p12[j, m] := begin
        if j ≠ m
            -Γ * (σᶻ[j] + σᶻ[m] + 2σᶻσᶻ[j, m]) + im * (Ω⁺[j] * σᶻσ⁻[m, j] + Ω⁺[m] * σᶻσ⁻[j, m] - conj(Ω⁺[j] * σᶻσ⁻[m, j] + Ω⁺[m] * σᶻσ⁻[j, m])) + 2Γ * Γⱼₘ[j, m] * (σ⁺σ⁻[j, m] + conj(σ⁺σ⁻[j, m]))
        else
            zero(ComplexF64)
        end
    end
    dₜ_σᶻσᶻ_p2 = @tullio p13[j, m] := begin
        if (k ≠ j) && (k ≠ m) && (j ≠ m)
            -2Γ * real(G[j, k] * σ⁺σᶻσ⁻[j, m, k] + G[m, k] * σᶻσ⁺σ⁻[j, m, k])
        else
            zero(ComplexF64)
        end
    end
    dₜ_σᶻσᶻ = dₜ_σᶻσᶻ_p1 + dₜ_σᶻσᶻ_p2 # <-----

    du[:] .= vcat(dₜ_σ⁻, dₜ_σᶻ,
        reshape(dₜ_σᶻσ⁻, (N^2)),
        reshape(dₜ_σ⁺σ⁻, (N^2)),
        reshape(dₜ_σ⁻σ⁻, (N^2)),
        reshape(dₜ_σᶻσᶻ, (N^2))
    )
    dₜ_σᶻ_p1 = dₜ_σᶻ_p2 = dₜ_σ⁻_p1 = dₜ_σ⁻_p2 = 1
    dₜ_σᶻσ⁻_p1 = dₜ_σᶻσ⁻_p2 = 1
    dₜ_σ⁺σ⁻_p1 = dₜ_σ⁺σ⁻_p2 = dₜ_σ⁻σ⁻_p1 = dₜ_σ⁻σ⁻_p2 = dₜ_σᶻσᶻ_p1 = dₜ_σᶻσᶻ_p2 = 1
    σ⁺σ⁻σ⁻ = σᶻσᶻσ⁻ = σ⁺σ⁻σᶻ = σᶻσ⁻σ⁻ = σ⁺σᶻσ⁻ = 1
    nothing
end
function truncation(A, B, C, AB, AC, BC)
    @tullio ABC[j, m, k] := begin
        AB[j, m] * C[k] + AC[j, k] * B[m] + BC[m, k] * A[j] - 2A[j] * B[m] * C[k]
    end
end

