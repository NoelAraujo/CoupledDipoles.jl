b₀_of(problem) =
    σ₀_of(problem.physic) * ρ_of(problem.atoms) * geometric_integral(problem.atoms)

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
    Volume_of(atoms::Shape{Sphere})
"""
Volume_of(atoms::Atom{Sphere}) = (4π / 3) * atoms.sizes^3 # sizes == radius of Sphere


# these "geometric_integral" values holds only for homogeneous distributions
geometric_integral(atoms::Atom{Sphere}) = 2 * atoms.sizes # diameter = 2*Radius
geometric_integral(atoms::Atom{Cube}) = atoms.size # diameter = cube side = L


"""
    estimate_E₀(laser)
"""
estimate_E₀(laser) = √(laser.s * (1 + 4(laser.Δ / Γ)^2) / 2)


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