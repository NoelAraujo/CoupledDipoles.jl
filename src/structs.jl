#######################################################
### ------------------ DEFINITIONS -------------------
#######################################################

abstract type Dimension end
abstract type TwoD <: Dimension end
abstract type ThreeD <: Dimension end

struct Square <: TwoD end
struct Circle <: TwoD end
struct Cube <: ThreeD end
struct Sphere <: ThreeD
    isGaussian::Bool

    function Sphere(args...; kwargs...)
        isGaussian = get(kwargs, :gaussian, false)
        if isGaussian
            return new(true)
        else
            return new(false)
        end
    end
end
struct Cylinder <: ThreeD end

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
end
struct PlaneWave3D <: PlaneWave
end
struct Gaussian2D <: Gaussian
    w₀::Float64
end
struct Gaussian3D <: Gaussian
    w₀::Float64
    function Gaussian3D(w₀)
        λ = 2π / k₀
        if w₀ < 2λ
            @warn "Waist is Too Small, it should be larger than 2λ (~ $(round(2λ, digits=3))).
            \n If the input value is correct, at least make sure that `w₀` is larger than `system size`.
            \n This advice garantee that Transmission curves will have values below 1." maxlog = 50
        end
        return new(w₀)
    end
end

mutable struct Laser{T<:Pump}
    pump::T
    s::Float64
    Δ::Float64
    direction::AbstractArray
    polarization::AbstractArray
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