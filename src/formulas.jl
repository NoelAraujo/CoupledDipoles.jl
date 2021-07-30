"""
    ρ_of(atoms::Shape{T})
Computes Volumetric Density
"""
ρ_of(atoms::Atom{T}) where T = atoms.N / Volume_of(atoms)


"""
    Volume_of(atoms::Shape{Cube})
"""
Volume_of(atoms::Atom{Cube}) = π * atoms.sizes^2 # sizes == length of Cube's Sides

"""
    Volume_of(atoms::Shape{Cube})
"""
Volume_of(atoms::Atom{Circle}) = atoms.sizes^3 # sizes == length of Cube's Sides


"""
    b₀_of(atoms::Shape{Cube})
"""
b₀_of(atoms::Atom{Cube}) = (ρ_of(atoms) .^ 2 * atoms.N / ((4π / 3) * (3 / (16π))^3))^(1 / 3)

"""
    estimate_E₀(laser)
"""
estimate_E₀(laser) = √(laser.s * (1 + 4(laser.Δ / Γ)^2) / 2)