b₀_of(problem) = σ₀_of(problem.physic)*ρ_of(problem.atoms)*geometric_integral(problem.atoms)

σ₀_of(problem::Scalar) = 4π/k₀^2
σ₀_of(problem::Vectorial) = 6π/k₀^2


"""
    ρ_of(atoms::Shape{T})
Computes Volumetric Density
"""
ρ_of(atoms::Atom{T}) where T = atoms.N / Volume_of(atoms)


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
Volume_of(atoms::Atom{Sphere}) = (4π/3)*atoms.sizes^3 # sizes == radius of Sphere


# these "geometric_integral" values holds only for homogeneous distributions
geometric_integral(atoms::Atom{Sphere}) = 2*atoms.sizes # diameter = 2*Radius
geometric_integral(atoms::Atom{Cube}) = atoms.size # diameter = cube side = L


"""
    estimate_E₀(laser)
"""
estimate_E₀(laser) = √(laser.s * (1 + 4(laser.Δ / Γ)^2) / 2)