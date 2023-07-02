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
function default_initial_condition(problem::NonLinearOptics{PairCorrelation})
    β₀ = zeros(ComplexF64, problem.atoms.N)
    z₀ = -ones(ComplexF64, problem.atoms.N)
    u₀ = vcat(β₀, z₀, zeros(ComplexF64, 4*problem.atoms.N^2))
    # u₀ = zeros(ComplexF64, 2*problem.atoms.N + 4*problem.atoms.N^2)
    return u₀
end


"""
    steady_state(problem::NonLinearOptics{MeanField})
"""
function steady_state(problem::NonLinearOptics{MeanField}; tmax=250.0, reltol=1e-10, abstol=1e-10, m=90, bruteforce=false)
    G = interaction_matrix(problem)
    if bruteforce
        tspan = (0.0, tmax)
        u₀ = default_initial_condition(problem)
        return time_evolution(problem, u₀, tspan, G; reltol=reltol, abstol=abstol, save_on=false).u[end]
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


function time_evolution(problem::NonLinearOptics{PairCorrelation}, u₀, tspan::Tuple; kargs...)

    G = interaction_matrix(problem)
    G_c = (2 / Γ).*G ## I NEED IT
	G_c = im.*G_c ## TO VERIFY IF I NEED IT
    solution = time_evolution(problem, u₀, tspan, G_c; kargs...)

    return solution
end
function time_evolution(problem::NonLinearOptics{PairCorrelation}, u₀, tspan::Tuple, G::AbstractMatrix; kargs...)
    N = problem.atoms.N
    Gconj = conj.(G)
    Γⱼₘ = real.(G)

    # `laser_field = (-im/2)*Ω`, but I need only `Ω`
    Ω⁻ = laser_field(problem.laser, problem.atoms) / LASER_FACTOR
    Ω⁺ = conj.(Ω⁻)

    Δ = problem.laser.Δ

    parameters = N, view(G,:,:), view(Gconj,:,:), view(Γⱼₘ,:,:), view(Ω⁻,:), view(Ω⁺,:), Δ

    ### calls for solver
    problemFunction = get_evolution_function(problem)
    prob = ODEProblem(problemFunction, u₀, tspan, parameters)
    solution = OrdinaryDiffEq.solve(prob, RDPK3SpFSAL35(); kargs...)


    return solution
end
get_evolution_function(problem::NonLinearOptics{PairCorrelation}) = PairCorrelation!

function PairCorrelation!(du, u, params, t; Γ=1)
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

    # σ⁻σᶻ = transpose(σᶻσ⁻) # original
	# σ⁺σᶻ = σᶻσ⁻' # original

    σ⁺σᶻ = transpose(σᶻσ⁻) # from Nicolas Morazotti
    # from Nicolas Morazotti
    σ⁻σᶻ = @tullio σ[j,m] := begin
		if m==j
			-σᶻσ⁻[j,m]
		else
			conj(σᶻσ⁻[j,m])
		end
	end
    σᶻσ⁺ = transpose(σ⁻σᶻ) ## aqui

	σ⁺σ⁻σ⁻ = truncation(σ⁺, σ⁻, σ⁻, σ⁺σ⁻, σ⁺σ⁻, σ⁻σ⁻)
	σᶻσᶻσ⁻ = truncation(σᶻ, σᶻ, σ⁻, σᶻσᶻ, σᶻσ⁻, σᶻσ⁻)
	σ⁺σ⁻σᶻ = truncation(σ⁺, σ⁻, σᶻ, σ⁺σ⁻, σ⁺σᶻ, σ⁻σᶻ)
	σᶻσ⁻σ⁻ = truncation(σᶻ, σ⁻, σ⁻, σᶻσ⁻, σᶻσ⁻, σ⁻σ⁻)
	σ⁺σᶻσ⁻ = truncation(σ⁺, σᶻ, σ⁻, σ⁺σᶻ, σ⁺σ⁻, σᶻσ⁻)
	σᶻσ⁺σ⁻ = truncation(σᶻ, σ⁺, σ⁻, σᶻσ⁺, σᶻσ⁻, σ⁺σ⁻) # computing the missing term
    # σᶻσ⁺σ⁻ = σ⁺σ⁻σᶻ # original

	dₜ_σ⁻ = @tullio σ[j] := begin
		if m≠j
			(Γ/2) * G[j,m] .* σᶻσ⁻[j,m]
		else
			(im*Δ - Γ/2)*σ⁻[j] + im*(Ω⁺[j]/2)*σᶻ[j]
		end
	end
	dₜ_σᶻ = @tullio σ[j] := begin
		if m≠j
			-Γ*(G[j,m]*σ⁺σ⁻[j,m] + conj(G[j,m] * σ⁺σ⁻[j,m]))
		else
			im*(Ω⁻[j]*σ⁻[j] - conj(Ω⁻[j]*σ⁻[j])) - Γ*(1 + σᶻ[j])
		end
	end

	dₜ_σᶻσ⁻_p1 = @tullio σ[j,m] := begin
			(im*Δ - 3Γ/2)*σᶻσ⁻[j,m] - Γ*σ⁻[m] + im*(Ω⁻[j]*σ⁻σ⁻[j,m] - Ω⁺[j]*σ⁺σ⁻[j,m] + 0.5*Ω⁺[m]*σᶻσᶻ[j,m]) - Γ*Γⱼₘ[j,m]*σ⁻σᶻ[j,m]  - (Γ/2)* Gconj[j,m]*σ⁻[j]
		end
    dₜ_σᶻσ⁻_p2 = @tullio σ[j,m] := begin
			if (k≠j) && (k≠m)
				- Γ*( G[j,k]*σ⁺σ⁻σ⁻[j,m,k] + Gconj[j,k]*σ⁺σ⁻σ⁻[k,m,j] ) +  0.5Γ*( G[m,k]*σᶻσᶻσ⁻[m,j,k]  )
			else
                zero(eltype(u))
			end
		end
    dₜ_σᶻσ⁻ = dₜ_σᶻσ⁻_p1 + dₜ_σᶻσ⁻_p2

    dₜ_σ⁺σ⁻_p1 = @tullio σ[j,m] := begin
        -Γ*σ⁺σ⁻[j,m] - 0.5im*(Ω⁻[j]*σᶻσ⁻[j,m] - Ω⁺[m]*σ⁺σᶻ[j,m])+(Γ/4)*(G[j,m]*σᶻ[m] + Gconj[j,m]*σᶻ[j]) + (Γ/2)*Γⱼₘ[j,m]*σᶻσᶻ[j,m]
    end
    dₜ_σ⁺σ⁻_p2 = @tullio σ[j,m] := begin
        if (k≠j) && (k≠m)
            +(Γ/2)*(Gconj[j,k]*σ⁺σ⁻σᶻ[k,m,j] + G[m,k]*σ⁺σ⁻σᶻ[j,k,m])
        else
            zero(ComplexF64)
        end
    end
    dₜ_σ⁺σ⁻ = dₜ_σ⁺σ⁻_p1 + dₜ_σ⁺σ⁻_p2


    dₜ_σ⁻σ⁻_p1 = @tullio σ[j,m] := begin
        (2im*Δ - Γ)*σ⁻σ⁻[j,m] + 0.5im*(Ω⁺[j]*σᶻσ⁻[j,m] + Ω⁺[m]*σᶻσ⁻[m,j])
    end
    dₜ_σ⁻σ⁻_p2 = @tullio σ[j,m] := begin
        if (k≠j) && (k≠m)
            +(Γ/2)*(G[j,k]*σᶻσ⁻σ⁻[j,m,k] + G[m,k]*σᶻσ⁻σ⁻[m,j,k])
        else
            zero(ComplexF64)
        end
    end
    dₜ_σ⁻σ⁻ = dₜ_σ⁻σ⁻_p1 + dₜ_σ⁻σ⁻_p2


    dₜ_σᶻσᶻ_p1 = @tullio σ[j,m] := begin
        -Γ * (σᶻ[j] + σᶻ[m] + 2σᶻσᶻ[j,m]) + im*(Ω⁻[j]*σᶻσ⁻[m,j] + Ω⁻[m]*σᶻσ⁻[j,m] - conj(Ω⁻[j]*σᶻσ⁻[m,j] + Ω⁻[m]*σᶻσ⁻[j,m])) + 2Γ * Γⱼₘ[j,m]*(σ⁺σ⁻[j,m] + conj(σ⁺σ⁻[j,m]))
    end
    dₜ_σᶻσᶻ_p2 = @tullio σ[j,m] := begin
        if (k≠j) && (k≠m)
            -Γ*(G[j,k]*σ⁺σᶻσ⁻[j,m,k] + G[m,k]*σᶻσ⁺σ⁻[j,m,k] + conj(G[j,k]*σ⁺σᶻσ⁻[j,m,k] + G[m,k]*σᶻσ⁺σ⁻[j,m,k]))
            # -2Γ*real(G[j,k]*σ⁺σᶻσ⁻[j,m,k] + G[m,k]*σ⁺σ⁻σᶻ[k,j,m])
        else
            zero(ComplexF64)
        end
    end
    dₜ_σᶻσᶻ = dₜ_σᶻσᶻ_p1 + dₜ_σᶻσᶻ_p2

    ## Diagonal must be zero.
    ## Take for example, σ⁻σ⁻: once on ground,
    ## the application of the down operator leads to zero
    dₜ_σᶻσ⁻[diagind(dₜ_σᶻσ⁻)] .= 0.0 + 0.0im
    dₜ_σ⁺σ⁻[diagind(dₜ_σ⁺σ⁻)] .= 0.0 + 0.0im
    dₜ_σ⁻σ⁻[diagind(dₜ_σ⁻σ⁻)] .= 0.0 + 0.0im
    dₜ_σᶻσᶻ[diagind(dₜ_σᶻσᶻ)] .= 0.0 + 0.0im

	du[:] .= vcat(dₜ_σ⁻, dₜ_σᶻ,
            reshape(dₜ_σᶻσ⁻,(N^2)),
            reshape(dₜ_σ⁺σ⁻,(N^2)),
            reshape(dₜ_σ⁻σ⁻,(N^2)),
            reshape(dₜ_σᶻσᶻ,(N^2))
            )

    σ⁺σ⁻σ⁻ = σᶻσᶻσ⁻ = σ⁺σ⁻σᶻ = σᶻσ⁻σ⁻ = σ⁺σᶻσ⁻ = 1
	nothing
end
function truncation(A, B, C, AB, AC, BC)
	@tullio σ[j, m, k] := begin
		AB[j, m]*C[k] + AC[j, k]*B[m] + BC[m, k]*A[j] - 2A[j]*B[m]*C[k]
	end
end