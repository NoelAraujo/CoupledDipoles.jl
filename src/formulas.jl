"""
    ρ_of(atoms::Shape{T})
Computes Volumetric Density
"""
ρ_of(atoms::Shape{T}) = atoms.N / Volume_of(atoms)


"""
    Volume_of(atoms::Shape{Cube})
"""
Volume_of(atoms::Shape{Cube}) = atoms.sizes^3 # sizes == length of Cube's Sides

"""
    b₀_of(atoms::Shape{Cube})
"""
b₀_of(atoms::Shape{Cube}) = (ρ_of(atoms) .^ 2 * atoms.N / ((4π / 3) * (3 / (16π))^3))^(1 / 3)