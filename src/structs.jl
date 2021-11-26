#######################################################
### ------------------ DEFINITIONS -------------------
#######################################################

abstract type Dimension end
abstract type TwoD <: Dimension end
abstract type ThreeD <: Dimension end

struct Square <: TwoD end
struct Circle <: TwoD end
struct Cube <: ThreeD end
struct Sphere <: ThreeD end

struct Atom{T<:Dimension}
    shape::T
    r::AbstractMatrix
    N::Int64
    sizes::Any
end

get_dimension(atom::Atom{T}) where {T<:TwoD} = 2
get_dimension(atom::Atom{T}) where {T<:ThreeD} = 3


abstract type Pump end
abstract type PlaneWave <: Pump end
abstract type Gaussian <: Pump end

struct PlaneWave2D <: PlaneWave
    direction::Vector
end
struct PlaneWave3D <: PlaneWave
    direction::Vector
end
struct Gaussian2D <: Gaussian
    w₀::Float64
end
struct Gaussian3D <: Gaussian
    w₀::Float64
    function Gaussian3D(w₀)
        λ = 2π/k₀
        if w₀ < 2λ
            @warn "Waist is Too Small, it should be larger than 2λ (~ $(round(2λ, digits=3))). 
            \n If the input value is correct, at least make sure that `w₀` is larger than `system size`.
            \n This advice garantee that Transmission curves will have values below 1." maxlog=50
        end
        return new(w₀)
    end
end


struct Laser{T<:Pump}
    pump::T
    s::Float64
    Δ::Float64
end
struct Lasers{T<:Pump}
    pump::Vector{T}
    s::Vector{Float64}
    Δ::Vector{Float64}
end





abstract type Physics end
abstract type Linear <: Physics end
abstract type NonLinear <: Physics end

struct Scalar <: Linear end
struct Vectorial <: Linear end
struct MeanField <: NonLinear end
struct BBGKY <: NonLinear end

struct LinearOptics{T<:Linear}
    physic::T
    atoms::Atom
    laser::Laser
    kernelFunction::Function
    spectrum::Dict
    data::Dict
end

struct NonLinearOptics{T<:NonLinear}
    physic::T
    atoms::Atom
    laser::Laser
    excitations::Dict
    data::Dict
end



#######################################################
### ------------------ CONSTRUCTORS -------------------
#######################################################



# """
#     b₀_of(atoms::Cube)

# Formula : `ρ^2*N/( (4π/3 * (3/(16π))^3))^(1/3)  )`
# """
# b₀_of(atoms::Cube) = (ρ_of(atoms) .^ 2 * atoms.N / ((4π / 3) * (3 / (16π))^3))^(1 / 3)
# """
#     ρ_of(atoms::Cube)
# Formula : `N/( kL^3 )`
# """
# ρ_of(atoms::Cube) = atoms.N / atoms.kL^3


# """
#     b₀_of(atoms::Cube)

# Formula : `ρ^2*N/( (4π/3 * (3/(16π))^3))^(1/3)  )`
# """
# b₀_of(atoms::Sphere) = (ρ_of(atoms) .^ 2 * atoms.N / ((4π / 3) * (3 / (16π))^3))^(1 / 3)
# ρ_of(atoms::Sphere) = atoms.N / VolumeSphere(atoms.kR)



# ### --------------- LASER ---------------
# """
#     PlaneWave_3D(direction=:z, s=1e-6, Δ=0)

# `direction` can be [:x, :y, :z]
# """
# function PlaneWave_3D(direction=:z, s=1e-6, Δ=0)
#     if direction == :x
#         matrix_slace = 1
#     elseif direction == :y
#         matrix_slace = 2
#     elseif direction == :z
#         matrix_slace = 3
#     else
#         @error("No support for this direction yet")
#     end
#     return PlaneWave(matrix_slace, s, Δ)
# end
# """
#     PlaneWave_2D(direction=:z, s=1e-6, Δ=0)

# `direction` can be [:x, :y]
# """
# function PlaneWave_2D(direction=:x, s=1e-6, Δ=0)
#     if direction == :x
#         matrix_slace = 1
#     elseif direction == :y
#         matrix_slace = 2
#     else
#         @error("No support for this direction yet")
#     end
#     return PlaneWave(matrix_slace, s, Δ)
# end
# """
