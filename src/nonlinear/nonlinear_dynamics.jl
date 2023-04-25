"""
    default_evolution_initial_condition(NonLinearOptics{MeanField})
β₀ = zeros(ComplexF64, atoms.N)
z₀ = 2β₀.*conj.(β₀) .- 1
"""
function default_initial_condition(problem::NonLinearOptics{MeanField})
    β₀ = zeros(ComplexF64, problem.atoms.N)
    z₀ = -ones(ComplexF64, problem.atoms.N)
    u₀ = vcat(β₀, z₀)
    return u₀
end

"""
    steady_state(problem::NonLinearOptics{MeanField})
"""
function steady_state(problem::NonLinearOptics{MeanField}; tmax=250.0, reltol=1e-7, abstol=1e-7, m=75, bruteforce=false)
    G = interaction_matrix(problem)
    if bruteforce
        tspan = (0.0, tmax)
        u₀ = default_initial_condition(problem)
        return time_evolution(problem, u₀, tspan, G; reltol=reltol, abstol=abstol, save_on=false).u[end] # evolve a little bit
    else
        #=
            I don't sum over diagonal elements during time evolution
        thus, to avoid an IF statement, I put a zero on diagonal
        =#
        saveDiag = diagind(G)
        G[diagind(G)] .= zero(eltype(G))

        # `laser_field = (-im/2)*Ω`, but I need only `Ω`
        Ωₙ::Vector{ComplexF64} = laser_field(problem.laser, problem.atoms) / LASER_FACTOR
        Wₙ = similar(Ωₙ)
        G_βₙ = similar(Ωₙ)
        temp1 = similar(Ωₙ)
        temp2 = similar(Ωₙ)
        parameters = view(G, :, :), view(Ωₙ, :), Wₙ, problem.laser.Δ, problem.atoms.N, G_βₙ, temp1, temp2

        ## `nlsolve` convergence is senstive to initial conditions
        ## therefore, i decided to make a small time evoltuion, and use the result as initial condition
        ## NOTE: this trick is usefull only for small N, for now, is N < 1200
        ## but this number was NOT obtained by systematic tests, and could be improved
        u₀ = default_initial_condition(problem)
        if problem.atoms.N < 1200
            tspan = (0.0, tmax)
            u₀::Vector{ComplexF64} = time_evolution(problem, u₀, tspan, G; reltol=reltol, abstol=abstol, save_on=false).u[end] # evolve a little bit
        end

        try
            solution = nlsolve((du,u)->MeanField!(du, u, parameters, 0.0), u₀, method = :anderson, m=m, autodiff = :forward)

            # !!!! restore diagonal !!!!
            G[diagind(G)] .= saveDiag
            return solution.zero
        catch

            ## For lower N (N < 100 ?), nlsolve does not converge (i don't know why)
            ## Instead of returning an error, I return the result from time evolution.
            @warn "Steady State may not be accurate. Consider increasing number of particles."
            if problem.atoms.N < 1200
                return u₀
            else
                tspan = (0.0, tmax)
                u₀ = default_initial_condition(problem)
                u₀_attempt = time_evolution(problem, u₀, tspan, G; reltol=reltol, abstol=abstol, save_on=false).u[end] # evolve a little bit
                # !!!! restore diagonal !!!!
                G[diagind(G)] .= saveDiag
                return u₀_attempt
            end
        end
    end
end

function time_evolution(problem::NonLinearOptics{MeanField}, u₀, tspan::Tuple; kargs...)

    G = interaction_matrix(problem)
    solution = time_evolution(problem, u₀, tspan, G; kargs...)

    return solution
end
function time_evolution(problem::NonLinearOptics{MeanField}, u₀, tspan::Tuple, G::AbstractMatrix; kargs...)
    #=
        I don't sum over diagonal elements during time evolution
     thus, to avoid an IF statement, I put a zero on diagonal
    =#
    saveDiag = diagind(G)
    G[diagind(G)] .= zero(eltype(G))

    # `laser_field = (-im/2)*Ω`, but I need only `Ω`
    Ωₙ = laser_field(problem.laser, problem.atoms) / LASER_FACTOR
    Wₙ = similar(Ωₙ)
    G_βₙ = similar(Ωₙ)
    temp1 = similar(Ωₙ)
    temp2 = similar(Ωₙ)

    parameters = view(G, :, :), view(Ωₙ, :), Wₙ, problem.laser.Δ, problem.atoms.N, G_βₙ, temp1, temp2

    ### calls for solver
    problemFunction = get_evolution_function(problem)
    prob = ODEProblem(problemFunction, u₀, tspan, parameters)
    solution = OrdinaryDiffEq.solve(prob, RDPK3SpFSAL35(); kargs...)

    # # !!!! restore diagonal !!!!
    # G[diagind(G)] .= saveDiag

    return solution
end
get_evolution_function(problem::NonLinearOptics{MeanField}) = MeanField!

function MeanField!(du, u, p, t)
    # parameters
    G, Ωₙ, Wₙ, Δ, N, G_βₙ, temp1, temp2 = p

    βₙ = @view u[1:N]
    zₙ = @view u[(N + 1):end]

    #= Code below is equivalento to =#
    # Wₙ[:] .= Ωₙ/2 + im*(G*βₙ - diagG.*βₙ) # don't forget the element wise multiplication of "diagG.*βₙ"
    # @. du[1:N] = (im*Δ - Γ/2)*βₙ + im*Wₙ*zₙ
    # @. du[N+1:end] = -Γ*(1 + zₙ) - 4*imag(βₙ*conj(Wₙ))

    # mul!(G_βₙ, G, βₙ) # == G_βₙ = G*βₙ
    # @simd for i ∈ eachindex(Wₙ)
    #     @inbounds Wₙ[i] = Ωₙ[i]/2 + im*(G_βₙ[i] - diagG[i]*βₙ[i])
    # end
    # @simd for i ∈ eachindex(βₙ)
    #     @inbounds du[i] = (im*Δ - Γ/2)*βₙ[i] + im*Wₙ[i]*zₙ[i]
    # end
    # @simd for i ∈ eachindex(zₙ)
    #     @inbounds du[i+N] = -Γ*(1 + zₙ[i]) - 4*imag(βₙ[i]*conj(Wₙ[i]))
    # end

    mul!(G_βₙ, G, βₙ) # == G_βₙ = G*βₙ
    # Wₙ .= Ωₙ / 2 .+ im * (G_βₙ - diagG .* βₙ) # if I don't clean the diagonal in callee function
    @. Wₙ = Ωₙ / 2 + im * G_βₙ
    @. temp1 = (im * Δ - Γ / 2) * βₙ + im * Wₙ * zₙ
    @. temp2 = -Γ * (1 + zₙ) - 4 * imag(βₙ * conj(Wₙ))
    du[1:N] .= temp1
    du[N+1:end] .= temp2

    return nothing
end
"""
    not optimized code
"""
function MeanField!_v2(du, u, p, t)
    # parameters
    G, Ωₙ, Wₙ, Δ, N, G_βₙ, temp1, temp2 = p

    βₙ = @view u[1:N]
    zₙ = @view u[(N + 1):end]

    for j in 1:N
        temp1[j] = (im * Δ - Γ / 2) * βₙ[j] + 0.5 * im * Ωₙ[j] * zₙ[j] + (-2 / Γ) * (Γ / 2) * sum(G[j, m] * βₙ[m] for m = 1:N if j ≠ m) * zₙ[j]
    end
    for j in 1:N
        temp2[j] =
            (im * conj(Ωₙ[j]) * βₙ[j] - im * Ωₙ[j] * conj(βₙ[j])) - Γ * (1 + zₙ[j]) -
            (-2 / Γ) * (sum(G[j, m] * βₙ[m] * conj(βₙ[j]) for m = 1:N if m ≠ j) + conj.(sum(G[j, m] * βₙ[m] * conj(βₙ[j]) for m = 1:N if m ≠ j)))
    end

    du[:] .= vcat(temp1, temp2)

    return nothing
end
