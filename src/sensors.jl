function get_sensors_ring(; num_pts=180, φ_inicial=0, φ_final=2π, kR=1, θ=π / 2)
    ϕ_range = [(φ_final-φ_inicial)*k/num_pts for k=1:num_pts]
    θ_range = ones(length(ϕ_range)) .* θ

    x = kR .* cos.(ϕ_range) .* sin.(θ_range)
    y = kR .* sin.(ϕ_range) .* sin.(θ_range)
    z = kR .* cos.(θ_range)
    return Array(transpose([x y z]))
end

"""
    get_sensors_sphereSurface_Fibonacci(; num_pts::Integer=100, R=5k₀, θ_lims=(0,π))

Function uses the `Fibonacci Sphere Algorithm` to create `num_pts` points over a sphere shell of radius `R`. This is a Legacy code, prefer the Lebedev points.


Returns each axis (in cartesian coordinates) in row: [x y z]

size([x y z]) = (3,num_pts)
"""
function get_sensors_sphereSurface_Fibonacci(; num_pts::Integer=100, R=5k₀, θ_lims=(0, π))
    indices = 1:num_pts
    θ_min, θ_max = θ_lims

    θ = acos.(range(cos(θ_min); stop=cos(θ_max), length=num_pts))
    ϕ = π * (1 + sqrt(5)) .* indices

    x, y, z = R .* cos.(ϕ) .* sin.(θ), R .* sin.(ϕ) .* sin.(θ), R .* cos.(θ)
    return Array(transpose([x y z]))
end

"""
    get_sensors_sphereSurface(; R=100k₀, leb_order=125)

return xyz, w

---

Returns points in a sphere of radius R and their corresponding weights, following Lebedev quadrature.

The number of points depends on the order, check Lebedev.jl docs.  

We used `leb_order = 125` to procude N=5294 points. That was enough for all inhouse tests.  
Sensors will be positioned in Far Field, therefore, `R` should be "large" - which dependent on the cloud size. 

**Note**: Increase this value if your system size is already ~ R=100k₀.
"""
function get_sensors_sphereSurface(; R=100k₀, leb_order=125)
    x, y, z, w = lebedev_by_order(leb_order)

    xyz = zeros(3, length(x))
    for i in eachindex(x)
        _inSpherical = cart2sph([x[i], y[i], z[i]])

        _inSpherical[3] = R * _inSpherical[3] # change sensor radius

        xyz[:, i] = sph2cart(_inSpherical)
    end

    return xyz, w
end
