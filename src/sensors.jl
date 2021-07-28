function get_sensors_ring(; num_pts=180, φ_inicial=0, φ_final=2π, kR=1, θ=π/2)
    ϕ_range = range(φ_inicial, φ_final; length=num_pts)
    θ_range = ones(length(ϕ_range)) .* θ

    x = kR .* cos.(ϕ_range) .* sin.(θ_range)
    y = kR .* sin.(ϕ_range) .* sin.(θ_range)
    z = kR .* cos.(θ_range)
    return [x y z]'
end
