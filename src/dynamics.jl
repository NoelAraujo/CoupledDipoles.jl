function steady_state(problem::LinearOptics{Scalar})
    @debug "start: get steady state"

    G = interaction_matrix(problem)
    Ωₙ = laser_field(problem.laser, problem.atoms)
    if problem.atoms.N > 1
        βₛ = -(G \ Ωₙ)
    else
        # The negative sign was NOT forgotten
        # after some math, I verified that it does not exist
        βₛ = Ωₙ / G[1]
    end

    @debug "end  : get steady state"
    return βₛ
end
function steady_state(problem::NonLinearOptics{MeanField})
    @debug "start: stedy state - NonLinearOptics"

    u₀ = default_initial_condition(problem)
    tspan = (0.0, 500)
    steady_state = time_evolution(problem, u₀, tspan; save_on=false)

    @debug "end  : stedy state - NonLinearOptics"
    return steady_state.u[end]
end

function time_evolution(problem::LinearOptics{T}, u₀, tspan::Tuple; kargs...) where {T<:Linear}
    ### use default G and Ωₙ
    G = copy(interaction_matrix(problem))
    Ωₙ = laser_field(problem.laser, problem.atoms)
    solution = time_evolution(problem, u₀, tspan, Ωₙ, G; kargs...)

    return solution
end
function time_evolution(problem::LinearOptics{T}, u₀, tspan::Tuple, Ωₙ::Vector; kargs...) where {T<:Linear}
    ### use default G
    G = copy(interaction_matrix(problem))
    solution = time_evolution(problem, u₀, tspan, Ωₙ, G; kargs...)

    return solution
end
function time_evolution(problem::LinearOptics{T}, u₀, tspan::Tuple, Ωₙ::Vector, G::Matrix; kargs...) where {T<:Linear}
    @debug "start: time evolution - LinearOptics"

    ### parameters == constant vector and matrices
    parameters = view(G, :, :), view(Ωₙ, :)

    ### calls for solver
    problemFunction = get_evolution_function(problem)
    prob = ODEProblem(problemFunction, u₀, tspan, parameters)
    solution = DifferentialEquations.solve(prob, VCABM(); dt=1e-10, abstol=1e-10, reltol=1e-10, kargs...)

    @debug "end  : time evolution - LinearOptics"
    return solution
end

"""
    default_evolution_initial_condition(::Scalar) = zeros(ComplexF64, N)
"""
function default_evolution_initial_condition(problem::LinearOptics{Scalar})
    return zeros(ComplexF64, problem.atoms.N) # I must use "zeros" and NOT an undef array - with trash data inside
end

get_evolution_function(problem::LinearOptics{Scalar}) = Scalar!

function Scalar!(du, u, p, t)
    G, Ωₙ = p

    #=
    Equivalent to:
        "du[:] = G*u + Ωₙ"
    But using inplace operation
    =#
    mul!(du, G, u)
    du .-= Ωₙ

    return nothing
end

# ### --------------- MEAN FIELD ---------------
function time_evolution(problem::NonLinearOptics{MeanField}, u₀, tspan::Tuple; kargs...)
    @debug "start: time evolution - NonLinearOptics"
    G = interaction_matrix(problem)

    #= 
        I don't sum over diagonal elements during time evolution
     thus, to avoid an IF statement, I put a zero on diagonal 
    =#
    saveDiag = diagind(G)
    G[diagind(G)] .= zero(eltype(G))

    # laser_field returns `(-im/2)*Ω`, but I need only `Ω`
    Ωₙ = laser_field(problem.laser, problem.atoms) / (-im / 2)
    Wₙ = zeros(eltype(Ωₙ), size(Ωₙ))
    G_βₙ = zeros(eltype(Ωₙ), size(Ωₙ))
    temp1 = similar(Ωₙ)
    temp2 = similar(Ωₙ)

    parameters = view(G, :, :), view(Ωₙ, :), Wₙ, problem.laser.Δ, problem.atoms.N, G_βₙ, temp1, temp2

    ### calls for solver
    problemFunction = get_evolution_function(problem)
    prob = ODEProblem(problemFunction, u₀, tspan, parameters)
    solution = DifferentialEquations.solve(prob, OwrenZen3(); abstol=1e-8, kargs...)

    # !!!! restore diagonal !!!!
    G[diagind(G)] .= saveDiag

    @debug "end  : time evolution - NonLinearOptics"
    return solution
end
get_evolution_function(problem::NonLinearOptics{MeanField}) = MeanField!

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
    Wₙ .= Ωₙ / 2 .+ im * G_βₙ
    temp1 .= (im * Δ - Γ / 2) * βₙ .+ im * Wₙ .* zₙ
    temp2 .= -Γ * (1 .+ zₙ) - 4 * imag.(βₙ .* conj.(Wₙ))
    du[:] .= vcat(temp1, temp2)

    return nothing
end

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
