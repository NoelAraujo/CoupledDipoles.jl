"""
    ring_on_space(;num_pts=180, φ_inicial=0, φ_final=2π, kR=1, θ=π/2)

"""
function ring_on_space(; num_pts=180, φ_inicial=0, φ_final=2π, kR=1, θ=π / 2)
    ϕ_range = range(φ_inicial, φ_final; length=num_pts)
    θ_range = ones(length(ϕ_range)) .* θ

    x, y, z = kR .* cos.(ϕ_range) .* sin.(θ_range),
    kR .* sin.(ϕ_range) .* sin.(θ_range),
    kR .* cos.(θ_range)

    sensors = SArray[]
    for n in 1:num_pts
        push!(sensors, SVector{3}([x[n] y[n] z[n]]))
    end
    return sensors
end

"""
    ring_on_plane(;num_pts=180, φ_inicial=0, φ_final=2π, kR=1)

"""
function ring_on_plane(; num_pts=180, φ_inicial=0, φ_final=2π, kR=1)
    xyz = ring_on_space(;
        num_pts=num_pts, φ_inicial=φ_inicial, φ_final=φ_final, kR=kR, θ=π / 2
    )
    sensors = SArray[]
    for n in 1:num_pts
        push!(sensors, SVector{2}([xyz[n][1] xyz[n][2]]))
    end
    return sensors
end

function get_number_sensors(sensors)
    if typeof(sensors)==Vector{StaticArrays.SArray}
        return length(sensors)
    elseif typeof(sensors)==StaticArrays.SVector{2, Float64} || typeof(sensors)==StaticArrays.SVector{3, Float64}
        return 1 # single sensor case
    end
end
get_one_sensor(sensors, i) = sensors[i]

select_sensor_axes(sensors, direction) = select_matrix_axes(sensors, direction)
