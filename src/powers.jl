"""
    scattered_power(problem, atomic_states; part=:total)

driver function to compute scattered power

- part = [:total] | :coherent | :incoherent
"""
function scattered_power(problem, atomic_states; part=:total)
    if part==:total
        P_coh = _coherent_power(problem, atomic_states)
        P_inc = _incoherent_power(problem, atomic_states)
        return P_coh + P_inc
    elseif part==:coherent
        P_coh = _coherent_power(problem, atomic_states)
        return P_coh
    elseif part==:incoherent
        P_inc = _incoherent_power(problem, atomic_states)
        return P_inc
    end
end


sum_upper_diagonal(beta, r) = @tullio s := begin
    if m > j
        rjm = sqrt((r[j, 1] - r[m, 1])^2 + (r[j, 2] - r[m, 2])^2 + (r[j, 3] - r[m, 3])^2)
        2*real((sin(rjm) / rjm) * beta[j] * conj(beta[m]))
    else
        0.0
    end
end


## SCALAR
function _coherent_power(problem::LinearOptics{Scalar}, atomic_states)
    N = problem.atoms.N
    r = transpose(problem.atoms.r)
    σ⁻ = view(atomic_states, 1:N)
    R = how_far_is_farField(problem.atoms)

    result = sum_upper_diagonal(σ⁻, r) # upper diagonal

    if N == 1 # problems with matrix format and not arrays
        result += abs2(σ⁻[1][1]) # diaognal
    else
        result += sum(abs2,  σ⁻) # diaognal
    end
    result *= 4π*(Γ^2)/(2k₀*R)

    return result
end
function _incoherent_power(problem::LinearOptics{Scalar}, atomic_states)
    return 0.0
end


## VECTORIAL
function _coherent_power(problem::LinearOptics{Vectorial}, atomic_states)
    # TO DO
    return nothing
end
function _incoherent_power(problem::LinearOptics{Vectorial}, atomic_states)
    # TO DO
    return nothing
end


## MEANFIELD
function _coherent_power(problem::NonLinearOptics{MeanField}, atomic_states)
    N = problem.atoms.N
    r = transpose(problem.atoms.r)
    σ⁻ = view(atomic_states, 1:N)
    R = how_far_is_farField(problem)

    result = sum_upper_diagonal(σ⁻, r) # upper diagonal
    result += sum(abs2, σ⁻) # diagonal
    result *= 4π*(Γ^2)/(2k₀*R)

    return result
end
function _incoherent_power(problem::NonLinearOptics{MeanField}, atomic_states)
	N = problem.atoms.N
    R = how_far_is_farField(problem)

    σ⁻ = view(atomic_states, 1:N)
    σᶻ = real.( view(atomic_states, N+1:2N) )
    result = sum(-abs2(σ⁻[j]) + (1 + σᶻ[j]) / 2 for j ∈ 1:N)
    result *= 4π*(Γ^2)/(2k₀*R)

    return result
end
