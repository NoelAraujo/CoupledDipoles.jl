#######################################################
### ------------------ DEFINITIONS -------------------
#######################################################

abstract type Dimension end
abstract type TwoD   <: Dimension end
abstract type ThreeD <: Dimension end

struct Square <: TwoD end
struct Circle <: TwoD end
struct Cube   <: ThreeD end
struct Sphere <: ThreeD end

struct Atom{T <: Dimension}
    shape::T
    r::Matrix
    N::Int64
    sizes::Any
end

get_dimension(atom::Atom{T}) where T <: TwoD = 2
get_dimension(atom::Atom{T}) where T <: ThreeD = 3


abstract type Pump end
abstract type PlaneWave <: Pump end
abstract type Gaussian  <: Pump end

struct PlaneWave2D <: PlaneWave 
    direction::Vector
end
struct PlaneWave3D <: PlaneWave 
    direction::Vector
end
struct Gaussian2D  <: Gaussian 
    w₀::Float64
end
struct Gaussian3D  <: Gaussian 
    w₀::Float64
end


struct Laser{T <: Pump}
    pump::T
    s::Float64
    Δ::Float64
end
struct Lasers{T <: Pump}
    pump::Vector{T}
    s::Vector{Float64}
    Δ::Vector{Float64}
end





abstract type Physics end
abstract type Linear    <: Physics end
abstract type NonLinear <: Physics end

struct Scalar    <: Linear end
struct Vectorial <: Linear end
struct MeanField <: NonLinear end
struct BBGKY     <: NonLinear end

struct LinearOptics{T <: Linear}
    physic::T
    atoms::Atom
    laser::Laser
    kernelFunction::Function
    spectrum::Dict
    data::Dict
end

struct NonLinearOptics{T <: NonLinear}
    physic::T
    atoms::Atom
    laser::Laser
    excitations::Dict
    data::Dict
end



#######################################################
### ------------------ CONSTRUCTORS -------------------
#######################################################

### --------------- SHAPES ---------------
## To create the Constructor, you need:
## (1) Create the function to create atoms on `atoms.jl` file
## (2) b₀_of(atoms::Shape)
## (3) ρ_of(atoms::Shape)
## (4) Base.size(x::Shape)
## (5) Base.show(io::IO, cloud::Shape)
# """
#     Cube(N::Integer, kL::Real)

# The cube goes from [-kL/2, kL/2] (with homogeneous distribution)
# """
# function Cube(N::Integer, kL::Real, createFunction=ftn_AtomsOnCube!::Function)
#     dimensions = 3
#     ρ = N / kL^3
#     rₘᵢₙ = get_rₘᵢₙ(ρ)
#     rₘᵢₙ = 0.1
#     # if rₘᵢₙ ≥ kL / 10
#     #     rₘᵢₙ = kL / 100
#     # end
#     r = get_atoms(dimensions, N, rₘᵢₙ; createFunction, kL)

#     return Cube(r, N, kL)
# end
# """
#     Cube(r::Matrix{Float64}, kL)

# I expect that: size(r) = (N atoms, N dimensions)
# """
# function Cube(r::Matrix{Float64}, kL)
#     N = size(r, 1)
#     dimensions = 3

#     r_Static = SArray[]
#     for n in 1:N
#         push!(r_Static, SVector{dimensions}(r[n, :]))
#     end
#     return Cube(r_Static, N, kL)
# end
# """
#     Cube(r::Vector{StaticArrays.SArray})

# I expect that: length(r) = N atoms; length(r[1]) =  N dimensions
# """
# function Cube(r::Vector{StaticArrays.SArray}, kL)
#     N = length(r)
#     return Cube(r, N, kL)
# end
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
#     Sphere(N::Integer, kR::Real)

# Sphere at the origin with radius kR (with homogeneous distribution)
# """
# function Sphere(N::Integer, kR::Real, createFunction=ftn_AtomsOnSphere!::Function)
#     dimensions = 3
#     ρ = N / VolumeSphere(kR)
#     rₘᵢₙ = get_rₘᵢₙ(ρ)
#     if rₘᵢₙ ≥ kR / 200
#         rₘᵢₙ = kR / 200
#     end
#     r = get_atoms(dimensions, N, rₘᵢₙ; createFunction, kR)

#     return Sphere(r, N, kR)
# end
# VolumeSphere(kR::Real) = (4π / 3) * kR^3
# """
#     Sphere(r::Matrix{Float64}, kR)

# Insert the Sphere Matrix `r`, where I assume the distributation centered at the origin, and radius kR.
# """
# function Sphere(r::Matrix{Float64}, kR)
#     N = size(r, 1)
#     dimensions = 3

#     r_Static = SArray[]
#     for n in 1:N
#         push!(r_Static, SVector{dimensions}(r[n, :]))
#     end

#     return Sphere(r_Static, N, kR)
# end
# """
#     Sphere(r::Vector{StaticArrays.SArray}, kR)

# Insert the Sphere Matrix `r`, where I assume the distributation centered at the origin, and radius kR.
# """
# function Sphere(r::Vector{StaticArrays.SArray}, kR)
#     N = size(r, 1)
#     return Sphere(r, N, kR)
# end
# """
#     b₀_of(atoms::Cube)

# Formula : `ρ^2*N/( (4π/3 * (3/(16π))^3))^(1/3)  )`
# """
# b₀_of(atoms::Sphere) = (ρ_of(atoms) .^ 2 * atoms.N / ((4π / 3) * (3 / (16π))^3))^(1 / 3)
# ρ_of(atoms::Sphere) = atoms.N / VolumeSphere(atoms.kR)

# ### --------------- PHYSICS ---------------
# """
#     ScalarProblem(atoms, laser)
# """
# function ScalarProblem(atoms, laser)
#     KernelFunction = get_ScalarKernelFunction(atoms)

#     return SimulationScalar(atoms, laser, KernelFunction, nothing, nothing, nothing, nothing)
# end
# function get_ScalarKernelFunction(atoms::T) where {T<:ThreeD}
#     return green_scalar!
# end
# """
#     MeanFieldProblem(atoms, laser)
# """
# function MeanFieldProblem(atoms, laser)
#     return SimulationMeanField(atoms, laser, nothing, nothing, nothing, nothing)
# end



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
#     Gaussian_3D(ω₀=1, s=1e-6, Δ=0)
# """
# function Gaussian_3D(ω₀=1, s=1e-6, Δ=0)
#     return Gaussian3D(ω₀, s, Δ)
# end

