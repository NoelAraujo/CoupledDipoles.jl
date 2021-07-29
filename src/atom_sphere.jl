"""
    Atom(T <: ThreeD, ::Int64, ::Union{Real, Integer})
"""
function Atom(geometry::Sphere, N::Int64, kR::Union{Real, Integer}; createFunction=ftn_AtomsOnSphere::Function)
    @debug "start: Shape - Sphere"
    
    dimensions = 3
    ρ = 3N / (4π*kR^3)
    rₘᵢₙ = get_rₘᵢₙ(ρ)
    if rₘᵢₙ ≥ kR / 10
        rₘᵢₙ = kR / 100
    end
    
    r = get_atoms(dimensions, N, rₘᵢₙ; createFunction, kR)
    
    @debug "end  : Shape - Sphere"
    return Atom(Sphere(), r, N, Float64(kR))
end


function ftn_AtomsOnSphere(; kwargs...)
    kR = kwargs[:kR]
    new_atom = zeros(3)
    new_atom[1] = 2π * rand() # azimuth
    new_atom[2] = asin(2 * rand() - 1) # elevation
    new_atom[3] = kR * (rand()^(1 ./ 3.0)) # radii
    new_atom[:] = sph2cart(new_atom)
    return new_atom
end
"""
    spherical_coordinate=[azimuth, elevation, r]

    azimuth = projection on XY-plane (in radians)
    "elevation" or "Polar" = projection on Z-axis (in radians)
    r = radius

	ref: https://www.mathworks.com/help/matlab/ref/sph2cart.html
"""
function sph2cart(spherical_coordinate)
    azimuth = spherical_coordinate[1]
    elevation = spherical_coordinate[2]
    radius = spherical_coordinate[3]
    x = radius * cos(elevation) * cos(azimuth)
    y = radius * cos(elevation) * sin(azimuth)
    z = radius * sin(elevation)
    return [x, y, z]
end
