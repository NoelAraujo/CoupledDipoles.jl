using LinearAlgebra

function scalar_3D(r_jm::Vector; k₀=1)
    r = norm(r_jm)
    G = exp(+im*k₀*r)/(im*k₀*r)
    return G
end


δ(μ, η) = float(μ == η)
function vectorial_3D(r_jm::Vector; k₀=1)
    G = zeros(ComplexF64, 3, 3)
    r = k₀*norm(r_jm)
    r2 = r^2
    r̂ = r_jm./norm(r_jm)
    for μ in 1:3, η in 1:3
        term1 = (3/2)cis(r)/r
        term2 = (δ(μ, η) - r̂[μ]*r̂[η])
        term3 = (δ(μ, η) - 3r̂[μ]*r̂[η])*(im/r - 1/r2)
        G[μ, η] = term1*( term2 + term3 )
    end
    return G
end

# examples
r_jm = [1,3,-4]
scalar_3D(r_jm)
vectorial_3D(r_jm)