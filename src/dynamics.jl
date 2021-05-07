function get_steady_state(problem::SimulationScalar)
    G = get_interaction_matrix(problem)
    Ωₙ = 0.5im * laser_over_atoms(problem.laser, problem.atoms)
    βₛ = G \ Ωₙ
    return βₛ
end

function get_steady_state(problem::SimulationMeanField)
    full_evolution = time_evolution(problem, time_max=50)
    return vcat(full_evolution.βₜ[end], full_evolution.zₜ[end])
end

### --------------- SCALAR ---------------
function time_evolution(problem::SimulationScalar;time_min = 0.0, time_max = 10, dt=1e-10, abstol=1e-10, reltol=1e-10)
    ############################################
    ### parameters == constant vector and matrices
    H = get_interaction_matrix(problem)
    Ωₙ = -0.5im*laser_over_atoms(problem.laser, problem.atoms)

    parameters = view(H,:,:), view(Ωₙ,:)
    
    ############################################
    ### initial conditions
    u₀ = get_initial_conditions(problem)
    tspan = (time_min, time_max)

    ############################################
    ### calls for solver
    prob = ODEProblem(Scalar!, u₀, tspan, parameters)
    solution = DifferentialEquations.solve(prob, VCABM(), adaptive=true, dt=dt, reltol=reltol, abstol=abstol)
    
    ############################################
    time_array, βₜ = extract_solution_from_Scalar_Problem(solution)
            
    return (time_array=time_array, βₜ=βₜ)
end
function get_initial_conditions(problem::SimulationScalar)
    λ, ψ = get_spectrum(problem)
    β₀ = ψ[:, end] # the last mode is the most subradiant == more loccalized
    return β₀
end
function Scalar!(du, u, p, t)
    G, Ωₙ = p

    # du[:] = G*u + Ωₙ
    BLAS.gemv!('N', ComplexF64(1.0), G, u, ComplexF64(0.0), du) 
    du[:] += Ωₙ
    return nothing
end
function extract_solution_from_Scalar_Problem(solution)
    nSteps = length(solution.u)
    time_array = [solution.t[i]  for i in 1:nSteps]
    βₜ = [solution.u[i]  for i in 1:nSteps]
    return time_array, βₜ
end


### --------------- MEAN FIELD ---------------
function time_evolution(problem::SimulationMeanField;time_min = 0.0, time_max = 10, dt=1e-10, abstol=1e-10, reltol=1e-10)
    ############################################
    ### parameters
    H = get_interaction_matrix(problem)

    ## -> I don't sum over diagonal elements during time evolution
    ## to avoid an IF statement, I put a zero on diagonal
    H[diagind(H)] .= zero(eltype(H))

    ## -> Also, the definition for Mean Field evolution needs
    ## +(Γ/2), where scalar kernel return -(Γ/2).
    ## Just multiply by -1 fixes
    H = -H

    Ωₙ = laser_over_atoms(problem.laser, problem.atoms)
    Wₙ = similar(Ωₙ)
    G_βₙ = similar(Ωₙ)
    parameters = view(H,:,:), view(diag(H),:), view(Ωₙ,:), Wₙ, problem.laser.Δ, problem.atoms.N, G_βₙ

    ############################################
    ### initial conditions
    u₀ = get_initial_conditions(problem)
    tspan = (time_min, time_max)
    
    ############################################
    ### calls for solver
    prob = ODEProblem(MeanField!, u₀, tspan, parameters)
    solution = DifferentialEquations.solve(prob, VCABM(), adaptive=true, dt=dt, reltol=reltol, abstol=abstol)
    
    ############################################
    time_array, βₜ, zₜ = extract_solution_from_MeanField_Problem(solution)
        
    return (time_array=time_array, βₜ=βₜ, zₜ=zₜ)
end
function get_initial_conditions(problem::SimulationMeanField)
    λ, ψ = get_spectrum(problem)
    β₀ = ψ[:, end] # the last mode is the most subradiant == more loccalized
    z₀ = 2β₀.*conj.(β₀) .- 1
    u₀ = vcat(β₀, z₀)
    return u₀
end
function MeanField!(du, u, p, t)
    # parameters
    G, diagG, Ωₙ, Wₙ, Δ, N, G_βₙ = p
    
    βₙ = @view u[1:N]
    zₙ = @view u[N+1:end]
    
    ### Code below is equivalento to
    # Wₙ[:] .= Ωₙ/2 - im*(G*βₙ - diagG.*βₙ) # don't forget the element wise multiplication of "diagG.*βₙ"
    # @. du[1:N] = (im*Δ - Γ/2)*βₙ + im*Wₙ*zₙ
    # @. du[N+1:end] = -Γ*(1 + zₙ) - 4*imag(βₙ*conj(Wₙ))
    
    BLAS.gemv!('N', ComplexF64(1.0), G, βₙ, ComplexF64(0.0), G_βₙ) # == G_βₙ = G*βₙ
    @simd for i = 1:N
        @inbounds Wₙ[i] = Ωₙ[i]/2 - im*(G_βₙ[i] - diagG[i]*βₙ[i])
    end
    @simd for i = 1:N
        @inbounds du[i] = (im*Δ - Γ/2)*βₙ[i] + im*Wₙ[i]*zₙ[i]
    end
    @simd for i = 1:N
        @inbounds du[i+N] = -Γ*(1 + zₙ[i]) - 4*imag(βₙ[i]*conj(Wₙ[i]))
    end

    return nothing
end
function extract_solution_from_MeanField_Problem(solution)
    nTimeSteps = length(solution.u)
    time_array = [solution.t[i]  for i in 1:nTimeSteps   ]
    
    N = length(solution.u[1])÷2
    βₜ = [solution.u[i][1:N]  for i in 1:nTimeSteps   ]
    zₜ = [solution.u[i][(N+1):end]  for i in 1:nTimeSteps   ]
    return time_array, βₜ, zₜ
end