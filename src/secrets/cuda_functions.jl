# ENV["JULIA_CUDA_MEMORY_POOL"] = "none"
ENV["JULIA_CUDA_USE_COMPAT"] = false
ENV["JULIA_CUDA_USE_BINARYBUILDER"] = true
using CUDA
# using OrdinaryDiffEq
# using LinearAlgebra

# import Pkg
# Pkg.add(name="CUDA", version="3.8")

function steady_state(problem::LinearOptics{Scalar})
    if ENV["COUPLED_DIPOLES_USE_GPU"]
        G = interaction_matrix(problem) |> CuArray
    else
        G = interaction_matrix(problem)
    end

    if problem.atoms.N > 1
        if ENV["COUPLED_DIPOLES_USE_GPU"]
            Ωₙ = vec(laser_field(problem, problem.atoms.r)) |> CuArray
        else
            Ωₙ = vec(laser_field(problem, problem.atoms.r))
        end

        βₛ = -(G \ Ωₙ)
    else
        # the conversion from matrix to vector, avoids an extra function to deal single atom
        Ωₙ = laser_field(problem, vec(problem.atoms.r))
        βₛ =   [-(Ωₙ / G[1])]   # I need a vector for single element
    end
    return βₛ |> Array
end
function steady_state(problem::LinearOptics{Vectorial})
    G = interaction_matrix(problem) |> CuArray
    Ωₙ = laser_field(problem, problem.atoms.r)
    ## Ωₙ_eff  = [all X - all Y - all Z]
    Ωₙ_eff = vcat(view(Ωₙ, 1, :), view(Ωₙ, 2, :), view(Ωₙ, 3, :))  |> CuArray
    if problem.atoms.N > 1
        βₛ = Array((G \ Ωₙ_eff))
    else
        βₛ = -(Ωₙ / G[1])
    end
    # transpose and NOT transpose conjugated
    # because i am just changing the array format to create
    # an effetive result
    N = problem.atoms.N
    βₛ_x = transpose(βₛ[1:N])
    βₛ_y = transpose(βₛ[N+1:2N])
    βₛ_z = transpose(βₛ[2N+1:3N])
    βₛ_eff = vcat(βₛ_x, βₛ_y, βₛ_z)
    return βₛ_eff
end

function CoupledDipoles.steady_state(problem::NonLinearOptics{MeanField})
    @debug "start: stedy state - NonLinearOptics"

    u₀ = default_initial_condition(problem) |> CuArray
    tspan = (0.0, 500)
    steady_state = CoupledDipoles.time_evolution(problem, u₀, tspan; progress=true, save_on=false)

    @debug "end  : stedy state - NonLinearOptics"
    return Array(steady_state.u[end])
end

function CoupledDipoles.time_evolution(problem::NonLinearOptics{MeanField}, u₀, tspan::Tuple; kargs...)
    @debug "start: time evolution - NonLinearOptics"
    G = interaction_matrix(problem)

    #=
        I don't sum over diagonal elements during time evolution
     thus, to avoid an IF statement, I put a zero on diagonal
    =#
    saveDiag = diagind(G)
    G[diagind(G)] .= zero(eltype(G))
    Ggpu = CuArray(G)
    # laser_field returns `(-im/2)*Ω`, but I need only `Ω`
    Ωₙ = aser_field(problem.laser, problem.atoms) / (-im / 2) |> CuArray
    Wₙ = similar(Ωₙ)
    G_βₙ = similar(Ωₙ)
    temp1 = similar(Ωₙ)
    temp2 = similar(Ωₙ)

    parameters = view(Ggpu, :, :), view(Ωₙ, :), Wₙ, problem.laser.Δ, problem.atoms.N, G_βₙ, temp1, temp2

    ### calls for solver
    problemFunction = get_evolution_function(problem)
    prob = ODEProblem(problemFunction, u₀, tspan, parameters)
    solution = DifferentialEquations.solve(prob, OwrenZen3(); reltol=1e-8, kargs...)

    # !!!! restore diagonal !!!!
    Ggpu = 1 # clear memory

    @debug "end  : time evolution - NonLinearOptics"
    return solution
end

# function my_pairwise(a, b=a)
#     na = length(a)
#     r = Array{Float64}(undef, na, na)

#     @inbounds for (j, bj) in enumerate(b), (i, ai) in enumerate(a)
#         # After debuggingm, I discovered the correct expression is 'bj - ai'.
#         # And NOT 'ai-bj' as one would expect
#         r[i, j] = bj - ai
#     end
#     return r
# end
# @views function CoupledDipoles.scattering_intensity(problem::NonLinearOptics{MeanField}, atomic_states, measurement_positions, scattering_func::Function)
#     _laser = problem.laser
#     _r = problem.atoms.r
#     _physics = problem.physic

#     _sensors = measurement_positions
#     _states = atomic_states
#     _func = scattering_func

#     Xₙₘ = my_pairwise(_r[1, :])
#     Yₙₘ = my_pairwise(_r[2, :])
#     Zₙₘ = my_pairwise(_r[3, :])
#     Xₙₘ[diagind(Xₙₘ)] .= 0
#     Yₙₘ[diagind(Yₙₘ)] .= 0
#     Zₙₘ[diagind(Zₙₘ)] .= 0

#     Xₙₘ_gpu, Yₙₘ_gpu, Zₙₘ_gpu = CuArray(Xₙₘ), CuArray(Yₙₘ), CuArray(Zₙₘ)
#     ii = ones(N, N)
#     ii[diagind(ii)] .= 0
#     ii_gpu = CuArray(ii)

#     cis_rₙₘ_gpu = CuArray(Array{ComplexF64}(undef, N, N))
#     βₙₘ_gpu = similar(cis_rₙₘ_gpu)

#     n_sensors = CoupledDipoles._get_number_elements(_sensors)
#     if n_sensors == 1
#         return CoupledDipoles._OnePoint_Intensity(_physics, _laser, _r, _sensors, _states, _func, Xₙₘ_gpu, Yₙₘ_gpu, Zₙₘ_gpu, ii_gpu)
#     else
#         scat_int = zeros(n_sensors)

#         for i in 1:n_sensors # REMOVE MULTI THREADING
#             scat_int[i] = CoupledDipoles._OnePoint_Intensity(_physics, _laser, _r, _sensors[:, i], _states, _func, Xₙₘ_gpu, Yₙₘ_gpu, Zₙₘ_gpu, ii_gpu, cis_rₙₘ_gpu, βₙₘ_gpu)
#         end
#         return scat_int
#     end
# end
# function CoupledDipoles._OnePoint_Intensity(physic::MeanField, laser, R⃗, sensor, β, scattering_func, Xₙₘ_gpu, Yₙₘ_gpu, Zₙₘ_gpu, ii_gpu, cis_rₙₘ_gpu, βₙₘ_gpu; k₀=1, Γ=1)
#     Ω = laser_field(laser, sensor) / (-0.5im)

#     r = norm(sensor)
#     n̂ = sensor / r
#     N = size(R⃗, 2)

#     σ⁻ = β[1:N]
#     σ⁺ = conj.(σ⁻)
#     σᶻ = β[(N + 1):end]

#     σ⁻_gpu = CuArray(σ⁻)
#     σ⁺_gpu = CuArray(σ⁺)

#     term1 = abs2(Ω) / 4
#     term2 = real(-im * Ω * (exp(-im * k₀ * r) / (im * k₀ * r)) * CoupledDipoles.ThreadsX.sum(σ⁺[j] * cis(+k₀ * (n̂ ⋅ R⃗[:, j])) for j in 1:N))
#     term3 = CoupledDipoles._term3(σ⁻_gpu, σ⁺_gpu, n̂, R⃗, N, Xₙₘ_gpu, Yₙₘ_gpu, Zₙₘ_gpu, ii_gpu, cis_rₙₘ_gpu, βₙₘ_gpu)
#     term4 = CoupledDipoles.ThreadsX.sum((1 + σᶻ[j]) / 2 for j in 1:N)
#     intensity_oneSensor = term1 + (Γ / 2) * term2 + (Γ / (2k₀ * r))^2 * (term3 + term4)

#     return real(intensity_oneSensor)
# end
# function CoupledDipoles._term3(σ⁻, σ⁺, n̂, R⃗, N, Xₙₘ_gpu, Yₙₘ_gpu, Zₙₘ_gpu, ii_gpu, cis_rₙₘ_gpu, βₙₘ_gpu; k₀=1)
#     cis_rₙₘ_gpu .= cis.(-k₀ .* (n̂[1] .* Xₙₘ_gpu + n̂[2] .* Yₙₘ_gpu + n̂[3] .* Zₙₘ_gpu))
#     βₙₘ_gpu .= (σ⁺ * transpose(σ⁻))
#     βₙₘ_gpu .= βₙₘ_gpu .* ii_gpu

#     intensity = sum(βₙₘ_gpu .* cis_rₙₘ_gpu)

#     return real(intensity)
# end
