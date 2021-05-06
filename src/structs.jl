#######################################################
### ------------------ DEFINITIONS -------------------
#######################################################

### --------------- SHAPES ---------------
abstract type TwoD end
abstract type ThreeD end

struct Square <: TwoD
    r::Any
    N::Integer
    kL::Number
end
struct Circle <: TwoD
    r::Any
    N::Integer
    kR::Number
end

struct Cube <: ThreeD
    r::Any
    N::Integer
    kL::Number
end
struct Sphere <: ThreeD
    r::Any
    N::Integer
    kR::Number
end

### --------------- PHYSICS ---------------
abstract type Scalar end
# abstract type VectorialProblem end
abstract type MeanFieldProblem end

mutable struct SimulationScalar <: Scalar
    atoms::Any
    laser::Any
    KernelFunction!::Function
    λ::Any
    ψ::Any
    ξ
    R1
end

mutable struct SimulationMeanField <: MeanFieldProblem
    atoms
    laser
    λ
    ψ
    β
    σᶻ
end



### --------------- LASER ---------------
abstract type Laser end
abstract type PlaneWaveL <: Laser end
abstract type GaussianL <: Laser end

mutable struct PlaneWave <: PlaneWaveL
    direction::Any
    s::Real
    Δ::Real
end
mutable struct Gaussian3D <: GaussianL
    # implement the direction when gaussian beam points along x/y direction
    # also, i had to change `apply_gaussian3D_field` for this case
    # direction

    ω₀::Real
    s::Real
    Δ::Real
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
"""
    Cube(N::Integer, kL::Real)

The cube goes from [-kL/2, kL/2] (with homogeneous distribution)
"""
function Cube(N::Integer, kL::Real, createFunction=ftn_AtomsOnCube!::Function)
    dimensions = 3
    ρ = N / kL^3
    rₘᵢₙ = get_rₘᵢₙ(ρ)
    rₘᵢₙ = 0.1
    # if rₘᵢₙ ≥ kL / 10
    #     rₘᵢₙ = kL / 100
    # end
    r = get_atoms(dimensions, N, rₘᵢₙ; createFunction, kL)

    return Cube(r, N, kL)
end
"""
    Cube(r::Matrix{Float64}, kL)

I expect that: size(r) = (N atoms, N dimensions)
"""
function Cube(r::Matrix{Float64}, kL)
    N = size(r, 1)
    dimensions = 3

    r_Static = SArray[]
    for n in 1:N
        push!(r_Static, SVector{dimensions}(r[n, :]))
    end
    return Cube(r_Static, N, kL)
end
"""
    Cube(r::Vector{StaticArrays.SArray})

I expect that: length(r) = N atoms; length(r[1]) =  N dimensions
"""
function Cube(r::Vector{StaticArrays.SArray}, kL)
    N = length(r)
    return Cube(r, N, kL)
end

b₀_of(atoms::Cube) = (ρ_of(atoms) .^ 2 * atoms.N / ((4π / 3) * (3 / (16π))^3))^(1 / 3)
ρ_of(atoms::Cube) = atoms.N / atoms.kL^3

"""
    Sphere(N::Integer, kR::Real)

Sphere at the origin with radius kR (with homogeneous distribution)
"""
function Sphere(N::Integer, kR::Real, createFunction=ftn_AtomsOnSphere!::Function)
    dimensions = 3
    ρ = N / VolumeSphere(kR)
    rₘᵢₙ = get_rₘᵢₙ(ρ)
    if rₘᵢₙ ≥ kR / 200
        rₘᵢₙ = kR / 200
    end
    r = get_atoms(dimensions, N, rₘᵢₙ; createFunction, kR)

    return Sphere(r, N, kR)
end
VolumeSphere(kR::Real) = (4π / 3) * kR^3
"""
    Sphere(r::Matrix{Float64}, kR)

Insert the Sphere Matrix `r`, where I assume the distributation centered at the origin, and radius kR.
"""
function Sphere(r::Matrix{Float64}, kR)
    N = size(r, 1)
    dimensions = 3

    r_Static = SArray[]
    for n in 1:N
        push!(r_Static, SVector{dimensions}(r[n, :]))
    end

    return Sphere(r_Static, N, kR)
end
"""
    Sphere(r::Vector{StaticArrays.SArray}, kR)

Insert the Sphere Matrix `r`, where I assume the distributation centered at the origin, and radius kR.
"""
function Sphere(r::Vector{StaticArrays.SArray}, kR)
    N = size(r, 1)
    return Sphere(r, N, kR)
end

b₀_of(atoms::Sphere) = (ρ_of(atoms) .^ 2 * atoms.N / ((4π / 3) * (3 / (16π))^3))^(1 / 3)
ρ_of(atoms::Sphere) = atoms.N / VolumeSphere(atoms.kR)

### --------------- PHYSICS ---------------
function ScalarProblem(atoms, laser)
    KernelFunction = get_ScalarKernelFunction(atoms)

    return SimulationScalar(atoms, laser, KernelFunction, nothing, nothing, nothing, nothing)
end
function get_ScalarKernelFunction(atoms::T) where {T<:ThreeD}
    return green_scalar!
end

function MeanFieldProblem(atoms, laser)
    return SimulationMeanField(atoms, laser, nothing, nothing, nothing, nothing)
end



### --------------- LASER ---------------
function PlaneWave_3D(direction=:z, s=1e-6, Δ=0)
    if direction == :x
        matrix_slace = 1
    elseif direction == :y
        matrix_slace = 2
    elseif direction == :z
        matrix_slace = 3
    else
        @error("No support for this direction yet")
    end
    return PlaneWave(matrix_slace, s, Δ)
end
function PlaneWave_2D(direction=:x, s=1e-6, Δ=0)
    if direction == :x
        matrix_slace = 1
    elseif direction == :y
        matrix_slace = 2
    else
        @error("No support for this direction yet")
    end
    return PlaneWave(matrix_slace, s, Δ)
end

function Gaussian_3D(ω₀=1, s=1e-6, Δ=0)
    return Gaussian3D(ω₀, s, Δ)
end

#######################################################
### -------------------- IO RELATED -------------------
#######################################################
Base.size(atoms::Cube) = atoms.kL
Base.size(atoms::Sphere) = atoms.kR

function make_short(x)
    if (abs(x) ≈ 0)
        return x
    elseif (abs(x) ≥ 1e3) || (abs(x) ≤ 1e-3)
        return @sprintf("%.2E", x)
    else
        return round(x; digits=2)
    end
end

function Base.show(io::IO, atoms::Cube)
    return print(
        io,
        Crayon(; foreground=:blue),
        "Cube with ",
        Crayon(; foreground=:red),
        "N=$(atoms.N) ",
        Crayon(; foreground=:blue),
        "atoms and size ",
        Crayon(; foreground=:red),
        "kL=$(make_short(size(atoms)))",
    )
end

function Base.show(io::IO, atoms::Sphere)
    return print(
        io,
        Crayon(; foreground=:blue),
        "Sphere with ",
        Crayon(; foreground=:red),
        "N=$(atoms.N) ",
        Crayon(; foreground=:blue),
        "atoms and radius ",
        Crayon(; foreground=:red),
        "kR=$(make_short(size(atoms)))",
    )
end

function Base.show(io::IO, laser::PlaneWaveL)
    return print(
        io,
        Crayon(; foreground=:blue),
        "Plane Wave with ",
        Crayon(; foreground=:red),
        "s=$(make_short(laser.s)) ",
        Crayon(; foreground=:blue),
        "and ",
        Crayon(; foreground=:red),
        "Δ=$(make_short(laser.Δ))",
    )
end

function Base.show(io::IO, laser::Gaussian3D)
    return print(
        io,
        Crayon(; foreground=:blue),
        "Gaussian Wave with ",
        Crayon(; foreground=:red),
        "ω₀=$(make_short(laser.ω₀))",
        Crayon(; foreground=:blue),
        ", ",
        Crayon(; foreground=:red),
        "s=$(make_short(laser.s)) ",
        Crayon(; foreground=:blue),
        "and ",
        Crayon(; foreground=:red),
        "Δ=$(make_short(laser.Δ))",
    )
end

function Base.show(io::IO, simulation::SimulationScalar)
    return print(
        io,
        Crayon(; foreground=:red),
        "Scalar ",
        Crayon(; foreground=:blue),
        "Problem with a ",
        Crayon(; foreground=:red),
        "$(typeof(simulation.atoms)) ",
        Crayon(; foreground=:blue),
        "of ",
        Crayon(; foreground=:red),
        "N=$(make_short(simulation.atoms.N)) ",
        Crayon(; foreground=:blue),
        "atoms ",
    )
end

function Base.show(io::IO, simulation::SimulationMeanField)
    return print(
        io,
        Crayon(; foreground=:red),
        "Mean Field ",
        Crayon(; foreground=:blue),
        "Problem with a ",
        Crayon(; foreground=:red),
        "$(typeof(simulation.atoms)) ",
        Crayon(; foreground=:blue),
        "of ",
        Crayon(; foreground=:red),
        "N=$(make_short(simulation.atoms.N)) ",
        Crayon(; foreground=:blue),
        "atoms ",
    )
end

