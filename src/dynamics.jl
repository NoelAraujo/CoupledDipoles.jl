function get_steady_state(problem::SimulationScalar)
    G = get_interaction_matrix(problem)
    Ωₙ = 0.5im * laser_over_atoms(problem.laser, problem.atoms)
    βₛ = G \ Ωₙ
    return βₛ
end

# function get_SteadyState(cloud_type, cloud_data; s=1e-6, time_max = 100, dt=1e-10)
#     solution_full = nonlinear_evolution(cloud_type, cloud_data; s=s, time_max = time_max, dt=dt)
#     solution = vcat(solution_full.βₜ, solution_full.σᶻₜ)
#     return solution
# end

### --------------- SCALAR ---------------
function time_evolution(problem::SimulationScalar;time_min = 0.0, time_max = 10, dt=1e-10, abstol=1e-10, reltol=1e-10)
    ############################################
    ### parameters == constant vector and matrices
    H = get_interaction_matrix(problem)
    Ωₙ = -0.5im*laser_over_atoms(problem.laser, problem.atoms)

    parameters = H, Ωₙ
    
    ############################################
    ### initial conditions
    u₀ = get_initial_conditions(problem)
    tspan = (time_min, time_max)

    ############################################
    ### calls for solver
    prob = ODEProblem(Scalar!, u₀, tspan, parameters)
    solution = DifferentialEquations.solve(prob, Tsit5(), adaptive=true, dt=dt, reltol=reltol, abstol=abstol)
    
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
     
    du[:] = G*u + Ωₙ
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

    parameters = H, Ωₙ, problem.laser.Δ, problem.atoms.N # simulation_parameters = (H=H, Ωₙ=Ωₙ, Δ=cloud_type.Δ, Γ=cloud_type.Γ)

    ############################################
    ### initial conditions
    u₀ = get_initial_conditions(problem)
    tspan = (time_min, time_max)
    
    ############################################
    ### calls for solver
    prob = ODEProblem(MeanField!, u₀, tspan, parameters)
    solution = DifferentialEquations.solve(prob, Tsit5(), adaptive=true, dt=dt, reltol=reltol, abstol=abstol)
    # solution = solve(prob, KenCarp4(autodiff=false), adaptive=true, dt=dt, reltol=reltol, abstol=abstol)
    
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
    G, Ωₙ, Δ, N = p
    
    βₙ = @view u[1:N]
    zₙ = @view u[N+1:end]
    
    Wₙ = Ωₙ/2 .- im*(G*βₙ - diag(G).*βₙ)
    du_temp1 = (im*Δ - Γ/2).*βₙ .+ im*Wₙ.*zₙ
    du_temp2 = -Γ.*(1 .+ zₙ) .- 4*imag(βₙ.*conj(Wₙ))
    
    du_temp = vcat(du_temp1, du_temp2)
    du[:] .= du_temp
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