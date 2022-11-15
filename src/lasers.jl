"""
    laser_field(laser, points)

    Compute -(im/2)*Ω.(points), where Ω is the laser (spatial distribution + parameters)
"""
function laser_field(laser::Laser{T}, sensor::AbstractVector) where T <: Pump
    Ω₀ = raby_frequency(laser)
    if laser.polarization == [0,0,0]
        return (-im/2)*Ω₀*_scalar_laser_field(laser, sensor)
    else
        return (-im/2)*Ω₀*_vectorial_laser_field(laser, laser.polarization, sensor)
    end
end
function laser_field(laser::Laser{T}, sensors::AbstractMatrix) where T <: Pump
    Ω₀ = raby_frequency(laser)
    if laser.polarization == [0,0,0]
        return map(eachcol(sensors)) do sensor
                    (-im/2)*Ω₀*_scalar_laser_field(laser, sensor)
                end
    else
        return map(eachcol(sensors)) do sensor
                    (-im/2)*Ω₀*_vectorial_laser_field(laser, laser.polarization, sensor)
                end
    end
end
function laser_field(laser::Laser{T}, atoms::Atom{D}) where {T <: Pump, D<:Dimension}
   return laser_field(laser, atoms.r)
end
