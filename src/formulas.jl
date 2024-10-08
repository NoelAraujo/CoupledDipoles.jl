b₀_of(problem) = σ₀_of(problem.physic) * ρ_of(problem.atoms) * geometric_integral(problem.atoms)

σ₀_of(problem::Scalar) = 4π / k₀^2
σ₀_of(problem::Vectorial) = 6π / k₀^2

"""
    ρ_of(atoms::Shape{T})
Computes Volumetric Density
"""
ρ_of(atoms::Atom{T}) where {T} = atoms.N / Volume_of(atoms)

"""
    Volume_of(atoms::Shape{Circle})
"""
Volume_of(atoms::Atom{Circle}) = π * atoms.sizes^2 # sizes == circle radius

"""
    Volume_of(atoms::Shape{Cube})
"""
Volume_of(atoms::Atom{Cube}) = atoms.sizes^3 # sizes == length of Cube's Sides

"""
    Volume_of(atoms::Shape{Cylinder})
"""
Volume_of(atoms::Atom{Cylinder}) = atoms.sizes[:h] * π * atoms.sizes[:R]^2

"""
    Volume_of(atoms::Shape{Sphere})
"""
function Volume_of(atoms::Atom{Sphere})
    if atoms.shape.isGaussian
        return (2π)^(3 / 2) * atoms.sizes^3 # Check tutorials folder
    else
        return (4π / 3) * atoms.sizes^3 # sizes == radius of Sphere
    end
end
#=
These "geometric_integral" values holds only for homogeneous distributions
- They represent the light path acroos the cloud -
For example, to cross a sphere, we need the diameter of the sphere (diameter = 2*radius)
=#
function geometric_integral(atoms::Atom{Sphere})
    if atoms.shape.isGaussian
        return √(2π) * atoms.sizes # "√2π R" comes from an integral. Check tutorials folder
    else
        return 2 * atoms.sizes # diameter = 2*radius
    end
end
geometric_integral(atoms::Atom{Cube}) = atoms.sizes # diameter = cube side = L

#=
(usually) we want to know the physics in the laser direction, that is, at z-direction,
thus, we want to know the light path across the height.
=#
function geometric_integral(atoms::Atom{Cylinder})
    @info "using the cylinder height as light path"
    return atoms.sizes[:h]
end

"""
    sph2cart([θ, ϕ, r])

    θ = "elevation" or "Polar" = projection on Z-axis (in radians)
    ϕ = azimuth = projection on XY-plane (in radians)
    r = radius

	ref: https://doc.sagemath.org/html/en/thematic_tutorials/vector_calculus/vector_calc_change.html
"""
function sph2cart(spherical_coordinate)
    θ = spherical_coordinate[1]
    ϕ = spherical_coordinate[2]
    r = spherical_coordinate[3]
    x = r * cos(ϕ) * sin(θ)
    y = r * sin(ϕ) * sin(θ)
    z = r * cos(θ)
    return [x, y, z]
end
function cart2sph(cartesian_coordinate)
    x = cartesian_coordinate[1]
    y = cartesian_coordinate[2]
    z = cartesian_coordinate[3]
    ϕ = atan(y, x)
    θ = atan(sqrt(x^2 + y^2), z)
    r = sqrt(x^2 + y^2 + z^2)
    return [θ, ϕ, r]
end


function laser_angles(direction::AbstractVector)
    θ = acos( direction[3]/norm(direction) )
    if (direction[2]^2 + direction[1]^2) ≠ 0
        ϕ = acos(direction[1]/sqrt(direction[2]^2 + direction[1]^2))
    else
        ϕ = 0.0
    end
    return θ, ϕ
end


function how_far_is_farField(problem)
    return how_far_is_farField(problem.atoms)
end
"""
    how_far_is_farField(atoms) = 50*size(atoms)^2
"""
function how_far_is_farField(atoms::Atom{T}) where T<: Dimension
    R = 50*size(atoms)^2
    return R
end

@views function inverseMatrix2x2(A)
    constDivision = 1/det(A)
    [A[4]*constDivision -A[3]*constDivision; -A[2]*constDivision A[1]*constDivision]
end
@views function inverseMatrix3x3(A)
    constDivision = 1/det(A)
    a1 = (A[2,2]*A[3,3] -A[2,3]*A[3,2])*constDivision
    a2 = (A[1,3]*A[3,2] -A[1,2]*A[3,3])*constDivision
    a3 = (A[1,2]*A[2,3] -A[1,3]*A[2,2])*constDivision
    a4 = (A[2,3]*A[3,1] -A[2,1]*A[3,3])*constDivision
    a5 = (A[1,1]*A[3,3] -A[1,3]*A[3,1])*constDivision
    a6 = (A[1,3]*A[2,1] -A[1,1]*A[2,3])*constDivision
    a7 = (A[2,1]*A[3,2] -A[2,2]*A[3,1])*constDivision
    a8 = (A[1,2]*A[3,1] -A[1,1]*A[3,2])*constDivision
    a9 = (A[1,1]*A[2,2] -A[1,2]*A[2,1])*constDivision
    return [a1 a2 a3; a4 a5 a6; a7 a8 a9]
end

single_sensor(x) = reshape(x, length(x), 1)