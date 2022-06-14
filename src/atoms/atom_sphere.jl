"""
    Atom(T <: ThreeD, ::Int64, ::Union{Real, Integer})
"""
function Atom(geometry::Sphere, N::Int64, kR::Union{Real,Integer})
    @debug "start: Shape - Sphere"

    dimensions = 3
    ρ = 3N / (4π * kR^3)
    rₘᵢₙ = get_rₘᵢₙ(ρ)
    if rₘᵢₙ ≥ kR / 10
        rₘᵢₙ = kR / 100
    end

    if geometry.isGaussian
        createFunction = ftn_AtomsOnSphere_Gaussian
    else
        createFunction = ftn_AtomsOnSphere
    end

    r = get_atoms(dimensions, N, rₘᵢₙ; createFunction, kR)

    @debug "end  : Shape - Sphere"
    return Atom(Sphere(; gaussian=geometry.isGaussian), r, N, Float64(kR))
end

function Atom(geometry::Sphere, r::Matrix, kR::Union{Real,Integer})
    N = size(r, 2) # remember to use each collum as a atom position
    return Atom(Sphere(), r, N, Float64(kR))
end

"""
    expected 'kwargs': kR
based upon: https://datagenetics.com/blog/january32020/index.html
"""
function ftn_AtomsOnSphere(; kwargs...)
    kR = kwargs[:kR]
    new_atom = zeros(3)
    U = rand()^(1 ./ 3.0)
    x = (2rand() - 1)
    y = (2rand() - 1)
    z = (2rand() - 1)
    mag = sqrt(x^2 + y^2 + z^2)
    
    new_atom[1] = kR*U*x/mag
    new_atom[2] = kR*U*y/mag
    new_atom[3] = kR*U*z/mag
    
    return new_atom
end
function ftn_AtomsOnSphere_Gaussian(; kwargs...)
    μ = [0.0, 0.0, 0.0]
    σ = get(kwargs, :kR, 1.0)
    new_atom = rand(MvNormal(μ, σ))

    return new_atom
end
