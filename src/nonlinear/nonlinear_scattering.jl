@views function _OnePoint_Intensity(physic::MeanField, laser, R⃗, sensor, β, scattering_func)
    N = size(R⃗, 2)

    r = norm(sensor)
    n = sensor ./ r

    σ⁻ = β[1:N]
    σᶻ = β[(N + 1):end]
    R = R⃗

    Ω = laser_field(laser, sensor) / (-0.5im)
    # y = -Γ / 2 * (cis(k₀ * r) / (im * k₀ * r)) * sum(σ⁻[j] * cis(-k₀ * (n ⋅ R[:, j])) for j in 1:N)

    y = zero(ComplexF64)
    dot_nR = zero(Float64)
    for j in 1:N
        dot_nR *= 0.0
        oneAtom = R[:, j]
        for i in eachindex(oneAtom)
        @inbounds    dot_nR += n[i]*R[i, j]
        end
        @inbounds y += σ⁻[j] * cis(-k₀ * (dot_nR))
    end
    y *= -Γ / 2 * (cis(k₀ * r) / (im * k₀ * r))

    E = Ω + y
    # intensity_oneSensor = conj(E) * E + (Γ^2) * (1 / (2 * k₀ * r)^2) * (-sum(σ⁻ .* conj.(σ⁻)) + sum(1 .+ σᶻ) / 2)

    intensity_oneSensor = conj(E) * E
    for j in 1:N
        intensity_oneSensor += (1 / (2 * k₀ * r)^2) * (-(σ⁻[j] * conj(σ⁻[j])) + (1 + σᶻ[j]) / 2)
    end
    return real(intensity_oneSensor)
end

function _OnePoint_Intensity_legacy(physic::MeanField, laser, R⃗, sensor, β, scattering_func)
    Ω = laser_field(laser, sensor) / (-0.5im)

    r = norm(sensor)
    n̂ = sensor / r
    N = size(R⃗, 2)

    σ⁻ = β[1:N]
    σ⁺ = conj.(σ⁻)
    σᶻ = β[(N + 1):end]

    term1 = abs2(Ω) / 4
    term2 = real(-im * Ω * (cis(-k₀ * r) / (im * k₀ * r)) * ThreadsX.sum(σ⁺[j] * cis(+k₀ * (n̂ ⋅ R⃗[:, j])) for j in 1:N))
    term3 = _term3(σ⁻, σ⁺, n̂, R⃗)
    term4 = ThreadsX.sum((1 + σᶻ[j]) / 2 for j in 1:N)

    intensity_oneSensor = term1 + (Γ / 2) * term2 + (Γ / (2k₀ * r))^2 * (term3 + term4)
    return real(intensity_oneSensor)
end
"""
    Hard core optimizations for term3. Check benchmarks folder for details.
"""
@views function _term3(σ⁻, σ⁺, n̂, R⃗)
    N = length(σ⁻)
    number_configurations = ((N^2) ÷ 2 - N ÷ 2)

    βₙₘ = Array{eltype(σ⁻)}(undef, number_configurations)
    cont = 1
    for n in 1:N
        for m in (n + 1):N
            βₙₘ[cont] = σ⁻[n] * σ⁺[m]
            cont += 1
        end
    end

    rₙₘ = Array{eltype(R⃗)}(undef, 3, number_configurations)
    cont = 1
    for n in 1:N
        r_n = R⃗[:, n]
        for m in (n + 1):N
            rₙₘ[1, cont] = r_n[1] - R⃗[1, m]
            rₙₘ[2, cont] = r_n[2] - R⃗[2, m]
            rₙₘ[3, cont] = r_n[3] - R⃗[3, m]
            cont += 1
        end
    end

    intensity = ThreadsX.mapreduce(+, 1:number_configurations) do cont
        (
            begin
                @inbounds dot_n_r = n̂[1] * rₙₘ[1, cont] + n̂[2] * rₙₘ[2, cont] + n̂[3] * rₙₘ[3, cont]
                @inbounds βₙₘ[cont] * cis(-k₀ * dot_n_r)
            end
        )
    end

    return 2real(intensity)
end

function get_intensity_over_an_angle(problem::NonLinearOptics{MeanField}, atoms_states::Vector, θ::Float64)
    @debug "start : get intensity over an angle - NonLinearOptics{MeanField}"

    if is_integration_const_term_available(problem)
        Gₙₘ = problem.data[:Gₙₘ]
    else
        Gₙₘ = _get_meanfield_constant_term(problem.atoms, θ)
        problem.data[:Gₙₘ] = Gₙₘ
    end

    N = problem.atoms.N
    β = view(atoms_states, 1:N)
    z = view(atoms_states, (N + 1):(2N))

    βₙₘ = transpose(β * β') # I have to do "transpose" and NOT "adjoint = complex+tranpose"
    βₙₘ[diagind(βₙₘ)] .= (1 .+ z) ./ 2

    # IMPORTANT FOR NEXT LINE: you should use ELEMENT WISE multiplication.
    # Also, you can reduce memory allocation with inplace multiplication
    βₙₘ .*= Gₙₘ
    intensity = real(sum(βₙₘ))

    @debug "end  : get intensity over an angle - NonLinearOptics{MeanField}"
    return intensity
end

function is_integration_const_term_available(problem)
    if haskey(problem.data, :Gₙₘ)
        return true
    else
        return false
    end
end

function _get_meanfield_constant_term(atoms, Θ)
    N, r = atoms.N, atoms.r

    xₙₘ, yₙₘ, zₙₘ = get_xyz_distances(r)
    k₀sinΘ = k₀ * sin(Θ)
    cos_Θ = cos(Θ)

    Gₙₘ = Array{ComplexF64,2}(undef, N, N)
    @sync for n in 1:N
        Threads.@spawn for m in 1:N # we had to compute all terms, and not the upper part
            @inbounds Gₙₘ[n, m] = _constant_term_core_computation(xₙₘ, yₙₘ, zₙₘ, n, m, cos_Θ, k₀sinΘ)
        end
    end

    #=
        IF I want to compute only the upper part,
        I have to multiply all terms by π/2:   Gₙₘ = (π/2)*Gₙₘ

        I decided to don't make this, to don't appear with factors
        not mentioned on theory.
    =#
    #= 
        Before returning, we HAVE to do some memory cleaning,      
        EVEN losing some performance. 

        Without this cleaning, the garbage collector gets lost 
        outside this function - when many simulation occurs at the same time.
    =#
    xₙₘ = yₙₘ = zₙₘ = 1
    GC.gc()  # DO NOT DELETE
    return Gₙₘ
end
function _constant_term_core_computation(xₙₘ, yₙₘ, zₙₘ, n, m, cos_Θ, k₀sinΘ)
    a = zero(ComplexF64)
    a = cis(k₀ * zₙₘ[n, m] * cos_Θ) * besselj(0, k₀sinΘ * sqrt(xₙₘ[n, m]^2 + yₙₘ[n, m]^2))
    return a
end
function get_xyz_distances(r) # @memoize  --> creating warning, let's ignore it right now
    dimensions = size(r, 1)
    N = size(r, 2)

    xₙₘ = Array{Float64,2}(undef, N, N)
    yₙₘ = Array{Float64,2}(undef, N, N)
    zₙₘ = Array{Float64,2}(undef, N, N)

    r_shared = Array{Float64,2}(undef, dimensions, N)
    r_shared .= r
    @sync for n in 1:N
        r_n = view(r_shared, :, n)
        Threads.@spawn for m in 1:N
            xₙₘ[n, m] = r_n[1] - r_shared[1, m]
            yₙₘ[n, m] = r_n[2] - r_shared[2, m]
            zₙₘ[n, m] = r_n[3] - r_shared[3, m]
        end
    end

    return xₙₘ, yₙₘ, zₙₘ
end
