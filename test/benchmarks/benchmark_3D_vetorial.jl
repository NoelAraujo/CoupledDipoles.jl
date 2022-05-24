#=
    For this particular positions
        r = [1 2 5; 3 5 -1/2]
    the expression below are exact solutions.

    Check sagemath notebook for the proof.
=#
xx = [1im -(3 * √(173) / 5177717)*(-250 * √(173) * im-26661)*exp(0.5 * im * √(173)); -(3 * √(173) / 5177717)*(-250 * √(173) * im-26661)*exp(0.5 * im * √(173)) 1im]
xy = [0 -(3 * √(173) / 5177717)*(144 * √(173) * im+3864)*exp(0.5 * im * √(173)); -(3 * √(173) / 5177717)*(144 * √(173) * im+3864)*exp(0.5 * im * √(173)) 0]
xz = [0 -(3 * √(173) / 5177717)*(-264 * √(173) * im-7084)*exp(0.5 * im * √(173)); -(3 * √(173) / 5177717)*(-264 * √(173) * im-7084)*exp(0.5 * im * √(173)) 0]

yx = [0 -(3 * √(173) / 5177717)*(144 * √(173) * im+3864)*exp(0.5 * im * √(173)); -(3 * √(173) / 5177717)*(144 * √(173) * im+3864)*exp(0.5 * im * √(173)) 0]
yy = [1im -(3 * √(173) / 5177717)*(-130 * √(173) * im-23441)*exp(0.5 * im * √(173)); -(3 * √(173) / 5177717)*(-130 * √(173) * im-23441)*exp(0.5 * im * √(173)) 1im]
yz = [0 -(3 * √(173) / 5177717)*(-396 * √(173) * im-10626)*exp(0.5 * im * √(173)); -(3 * √(173) / 5177717)*(-396 * √(173) * im-10626)*exp(0.5 * im * √(173)) 0]

zx = [0 -(3 * √(173) / 5177717)*(-264 * √(173) * im-7084)*exp(0.5 * im * √(173)); -(3 * √(173) / 5177717)*(-264 * √(173) * im-7084)*exp(0.5 * im * √(173)) 0]
zy = [0 -(3 * √(173) / 5177717)*(-396 * √(173) * im-10626)*exp(0.5 * im * √(173)); -(3 * √(173) / 5177717)*(-396 * √(173) * im-10626)*exp(0.5 * im * √(173)) 0]
zz = [1im -(3 * √(173) / 5177717)*(380 * √(173) * im-9756)*exp(0.5 * im * √(173)); -(3 * √(173) / 5177717)*(380 * √(173) * im-9756)*exp(0.5 * im * √(173)) 1im]

K_exact = [xx xy xz; yx yy yz; zx zy zz]

#=
        Creating individual elements
=#
let
    using LinearAlgebra
    N = 2
    r = [1 2 5; 3 5 -1/2]
    δ(x, y) = Int(==(x, y))

    P(x) = 1 - 1 / x + 1 / x^2
    Q(x) = -1 + 3 / x - 3 / x^2

    α_range = β_range = [-1, +1, 0]

    for (α_idx, α) in enumerate(α_range)
        for (β_idx, β) in enumerate(β_range)
            for j in 1:N
                for n in 1:N
                    if j ≠ n
                        r_jn_vec = r[j, :] - r[n, :]
                        r_jn = norm(r_jn_vec)

                        term1 = im .* δ(j, n) * δ(α, β)
                        term2 = (1 - δ(j, n)) * (3 / 2) * (exp(im * r_jn) / r_jn)

                        term3 = P(im * r_jn) * δ(α, β)
                        term4 = Q(im * r_jn) * r_jn_vec[α_idx] * r_jn_vec[β_idx] / r_jn^2
                    else
                        term1 = im .* δ(j, n) * δ(α, β)
                        term2 = term3 = term4 = 0
                    end
                    temp = term1 + term2 * (term3 + term4)
                    println("α:$(α) ... β:$(β) ... j:$(j) ... n:$(n) ... temp:$(temp)")
                end
            end
            println("---------")
        end
    end
end

#=
        Creating sub matrices but not storing them
=#
let
    using LinearAlgebra
    N = 2
    r = [1 2 5; 3 5 -1/2]
    δ(x, y) = Int(==(x, y))
    Xt, Yt, Zt = r[:, 1], r[:, 2], r[:, 3]

    Xjn = Xt * ones(1, N) - ones(N, 1) * Xt'
    Yjn = Yt * ones(1, N) - ones(N, 1) * Yt'
    Zjn = Zt * ones(1, N) - ones(N, 1) * Zt'

    array_XYZ_jn = [Xjn, Yjn, Zjn]
    Rjn = sqrt.(Xjn .^ 2 + Zjn .^ 2 + Yjn .^ 2)

    P(x) = 1 - 1 / x + 1 / x^2
    Q(x) = -1 + 3 / x - 3 / x^2

    α_range = β_range = [-1, +1, 0]

    for (α_idx, α) in enumerate(α_range)
        for (β_idx, β) in enumerate(β_range)
            term1 = im * I(N) .* δ(α, β)
            term2 = (3 / 2) * exp.(im * Rjn) ./ Rjn
            term3 = P.(im * Rjn) .* δ(α, β)
            term4 = (array_XYZ_jn[α_idx] .* array_XYZ_jn[β_idx]) .* Q.(im * Rjn) ./ Rjn .^ 2

            term2[findall(isnan.(term2))] .= 0
            term2[findall(isinf.(term2))] .= 0

            term3[findall(isnan.(term3))] .= 0
            term3[findall(isinf.(term3))] .= 0

            term4[findall(isnan.(term4))] .= 0
            term4[findall(isinf.(term4))] .= 0

            K = term1 + term2 .* (term3 .+ term4)
            println("------")
            display(K)
        end
    end
end

#=
        Creating matrices generically
=#
let
    using LinearAlgebra
    N = 2
    r = [1 2 5; 3 5 -1/2]
    δ(x, y) = Int(==(x, y))

    Xt, Yt, Zt = r[:, 1], r[:, 2], r[:, 3]

    Xjn = Xt * ones(1, N) - ones(N, 1) * Xt'
    Yjn = Yt * ones(1, N) - ones(N, 1) * Yt'
    Zjn = Zt * ones(1, N) - ones(N, 1) * Zt'

    array_XYZ_jn = [Xjn, Yjn, Zjn]
    Rjn = sqrt.(Xjn .^ 2 + Zjn .^ 2 + Yjn .^ 2)

    P(x) = 1 - 1 / x + 1 / x^2
    Q(x) = -1 + 3 / x - 3 / x^2

    α_range = β_range = [-1, +1, 0]

    term2 = (3 / 2) * exp.(im * Rjn) ./ Rjn
    term2[findall(isnan.(term2))] .= 0
    term2[findall(isinf.(term2))] .= 0

    P_Rjn = P.(im * Rjn)
    P_Rjn[findall(isnan.(P_Rjn))] .= 0

    Q_Rjn_over_Rjn_squared = Q.(im * Rjn) ./ Rjn .^ 2
    Q_Rjn_over_Rjn_squared[findall(isnan.(Q_Rjn_over_Rjn_squared))] .= 0

    A = []
    for (α_idx, α) in enumerate(α_range)
        B = []
        for (β_idx, β) in enumerate(β_range)
            term1 = im * I(N) .* δ(α, β)
            # term2 = (3/2)*exp.(im*Rjn)./Rjn  ## Defined outside      
            term3 = P_Rjn .* δ(α, β)
            term4 = Q_Rjn_over_Rjn_squared .* (array_XYZ_jn[α_idx] .* array_XYZ_jn[β_idx])

            K = term1 + term2 .* (term3 .+ term4)
            push!(B, K)
        end
        push!(A, vcat(B[1:length(β_range)]...))
    end

    K = hcat(A[1:length(α_range)]...)
    all(K .≈ K_exact)
end