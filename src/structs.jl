#######################################################
### ------------------ DEFINITIONS -------------------
#######################################################

abstract type Dimension end
abstract type TwoD <: Dimension end
abstract type ThreeD <: Dimension end

struct Square <: TwoD end
struct Circle <: TwoD end

"""
        Cube()
"""
struct Cube <: ThreeD end
"""
        Sphere(; gaussian=false)

If `gaussian=true`, produces a Gaussian Sphere with μ = 0 (mean = 0) and variance = kR (σ^2 = kR)
"""
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
"""
        Cylinder()
"""
struct Cylinder <: ThreeD end

struct Atom{T<:Dimension}
    shape::T            # Cube, Sphere, Cylinder
    r::Matrix{Float64}  # atomic matrix (column-major)
    N::Int64            # number of atoms
    sizes::Any          # characteristic of the shape (e.g, 'Sphere' is the 'radius')
end

get_dimension(atom::Atom{T}) where {T<:TwoD} = 2
get_dimension(atom::Atom{T}) where {T<:ThreeD} = 3

abstract type Pump end
abstract type PlaneWave <: Pump end
abstract type Gaussian <: Pump end

struct PlaneWave2D <: PlaneWave
end
"""
    PlaneWave3D()
"""
struct PlaneWave3D <: PlaneWave
end
struct Gaussian2D <: Gaussian
    w₀::Float64
end

"""
    Gaussian3D(w₀)

- w₀ is the Beam Waist.

Please, avoid w₀ < 2λ (λ = 2π / k₀) 
"""
struct Gaussian3D <: Gaussian
    w₀::Float64
    function Gaussian3D(w₀)
        λ = 2π / k₀
        if w₀ < 2λ
            @warn "Waist is Too Small, it should be larger than 2λ (~ $(round(2λ, digits=3))).
            \n If the input value is correct, at least make sure that `w₀` is larger than `system size`.
            \n This advice garantee that Transmission curves will have values below 1." maxlog = 1
        end
        return new(w₀)
    end
end
"""
    Laser(pump, s, Δ, direction, polarization)

- pump is e.g PlanweWave3D or Gaussian3D
- s: saturation on ressonance (used for raby_frequency = Γ √(s / 2))
- Δ: the detunning (laser - atomic frequency) 
- direction: Array is the propagation direction (default is [0,0,1])
- polarization: Array is the polarization direction (default is [0,0,0])

Note that `direction` and `polarization` have to be orthogonal.
"""
mutable struct Laser{T}         # 'mutable' because 'Δ' usually is is altered
    pump::T                     # PlanweWave3D, Gaussian3D
    s::Float64                  # saturation at ressonance
    Δ::Float64                  # atom-laser detunning
    direction::AbstractArray    # propagation direction
    polarization::AbstractArray # must be: orthogonal to 'direction'
end


abstract type Physics end
abstract type Linear <: Physics end
abstract type NonLinear <: Physics end

struct Scalar <: Linear end
struct Vectorial <: Linear end
struct MeanField <: NonLinear end
struct PairCorrelation <: NonLinear end

struct LinearOptics{T<:Linear}
    physic::T                   # Scalar, Vectorial
    atoms::Atom                 # result from Atom struct
    laser::Laser                # result from Laser struct
    kernelFunction!::Function   # computes the atomic interaction
    spectrum::Dict              # stores eigenvalues/eigenvectors
    data::Dict                  # empty dict for general use
end

struct NonLinearOptics{T<:NonLinear}
    physic::T
    atoms::Atom
    laser::Laser
    data::Dict
end