abstract type Dimension end
abstract type TwoD   <: Dimension end
abstract type ThreeD <: Dimension end

struct Square <: TwoD end
struct Circle <: TwoD end
struct Cube   <: ThreeD end
struct Sphere <: ThreeD end

struct Shape{T <: Dimension}
    dimension::T
    r::Any
    N::Integer
    sizes::Any
end
struct Shapes{T <: Dimension}
    dimension::Vector{T}
    r::Vector{<:Any}
    N::Vector{<:Integer}
    sizes::Vector{Any}
end




abstract type Pump end
abstract type PlaneWave <: Pump end
abstract type Gaussian  <: Pump end

struct PlaneWave2D <: PlaneWave 
    direction::Any
end
struct PlaneWave3D <: PlaneWave 
    direction::Any
end
struct Gaussian2D  <: Gaussian 
    w₀::Number
end
struct Gaussian3D  <: Gaussian 
    w₀::Number
end


struct Laser{T <: Pump}
    pump::T
    s::Real
    Δ::Real
end
struct Lasers{T <: Pump}
    pump::Vector{T}
    s::Vector{<:Real}
    Δ::Vector{<:Real}
end





abstract type Physics end
abstract type Linear    <: Physics end
abstract type NonLinear <: Physics end

struct Scalar    <: Linear end
struct Vectorial <: Linear end
struct MeanField <: NonLinear end
struct BBGKY     <: NonLinear end

mutable struct LinearProblem{T <: Linear}
    physic::T
    atoms::Shape
    laser::Laser
    kernelFunction!::Function
    spectrum::Tuple
    extra::Any
end
mutable struct LinearProblems{T <: Linear}
    physic::Vector{T}
    atoms::Shapes
    laser::Lasers
    kernelFunction!::Vector{Function}
    spectrum::Vector{Tuple}
    extra::Vector{<:Any}
end


mutable struct NonLinearProblem{T <: NonLinear}
    physic::T
    atoms::Shape
    laser::Laser
    excitations::Any
    extra::Any
end
mutable struct NonLinearProblems{T <: NonLinear}
    physic::Vector{T}
    atoms::Shapes
    laser::Lasers
    excitations::Vector{<:Any}
    extra::Vector{<:Any}
end

const u_input = Union{Real, Integer}

get_dimension(shape::Shape{T}) where T <: TwoD = 2
get_dimension(shape::Shape{T}) where T <: ThreeD = 3

get_kernelFunction(physics::Scalar,    dimension::Shape{<:TwoD}) = sin
get_kernelFunction(physics::Vectorial, dimension::Shape{<:TwoD}) = cos

get_kernelFunction(physics::Scalar,    dimension::Shape{<:ThreeD}) = sin
get_kernelFunction(physics::Vectorial, dimension::Shape{<:ThreeD}) = cos




cloud = Shape(Cube(), rand(30), 3,3)
ensemble = Shapes([Cube() for i=1:2],[rand(3) for i=1:2], [3 for i=1:2], [3.0 for i=1:2])
laser = Laser(Gaussian2D(1), 3, 2)

simulation = LinearProblem(Scalar(), cloud, laser, sin, (2,3), (2,3,4))