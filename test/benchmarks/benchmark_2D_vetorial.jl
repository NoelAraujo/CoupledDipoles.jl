#=
    For this particular positions
        r = [1 2; 3 5]
    the expression below are exact solutions.

    Check sagemath notebook for the proof.
=#
N = 2
r = [1 2; 3 5]

using SpecialFunctions, LinearAlgebra
H_0 = [1 besselh.(0, sqrt(13)); besselh.(0, sqrt(13)) 1]
H_m = [
    0 exp(-2im * atan(3, 2))*besselh.(2, sqrt(13))
    exp(-2im * atan(3, 2))*besselh.(2, sqrt(13)) 0
]

H_p = [
    0 exp(+2im * atan(3, 2))*besselh.(2, sqrt(13))
    exp(+2im * atan(3, 2))*besselh.(2, sqrt(13)) 0
]

K_exact = [H_0 H_p; H_m H_0]

#=
        Creating matrices individually
=#
let
    N = 2
    r = [1 2; 3 5]
    Dx = [r[j, 1] - r[l, 1] for j in 1:N, l in 1:N]
    Dy = [r[j, 2] - r[l, 2] for j in 1:N, l in 1:N]
    r_jl = sqrt.(Dx .^ 2 + Dy .^ 2)
    φ_jl = atan.(Dy, Dx)

    #Hankel functions
    H0 = besselh.(0, r_jl)
    H0[findall(isnan.(H0))] .= 0.0
    H0[diagind(H0)] .= 1

    Hp = exp.(+2 * 1im * φ_jl) .* besselh.(2, r_jl)
    Hp[diagind(Hp)] .= 0
    Hm = exp.(-2 * 1im * φ_jl) .* besselh.(2, r_jl)
    Hm[diagind(Hm)] .= 0

    K = [H0 Hp; Hm H0]
    K .≈ K_exact
end

#=
        Creating matrices generically
=#
let
    N = 2
    r = [1 2; 3 5]
    δ(x, y) = Int(==(x, y))
    Dx = [r[j, 1] - r[l, 1] for j in 1:N, l in 1:N]
    Dy = [r[j, 2] - r[l, 2] for j in 1:N, l in 1:N]
    r_jl = sqrt.(Dx .^ 2 + Dy .^ 2)
    φ_jl = atan.(Dy, Dx)

    H0 = besselh.(0, r_jl)
    H2 = besselh.(2, r_jl)

    A = []
    α_range = β_range = [-1, +1]
    for (α_idx, α) in enumerate(α_range)
        B = []
        for (β_idx, β) in enumerate(β_range)
            term1 = I(N) * δ(α, β)
            term2 = ones(N, N) - I(N)

            term3 = δ(α, β) * H0
            term3[findall(isnan.(term3))] .= 0
            term3[findall(isinf.(term3))] .= 0

            term4 = (1 - δ(α, β)) * exp.(2 * α * im * φ_jl) .* H2
            term4[findall(isnan.(term4))] .= 0

            temp = term1 + term2 .* (term3 + term4)

            push!(B, temp)
        end
        push!(A, vcat(B[1:length(β_range)]...))
    end

    K = hcat(A[1:length(α_range)]...)
    K .≈ K_exact
end
