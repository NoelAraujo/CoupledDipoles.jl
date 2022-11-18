function laser_field(problem::NonLinearProblem{MeanField}, sensor)
    Ω₀ = raby_frequency(problem.laser)
    return (-im/2)*Ω₀*_meanfield_laser_field(problem.laser, sensor)
end


## PUMP FIELD
function _meanfield_laser_field(laser::Union{Laser{PlaneWave3D},Laser{Gaussian3D}}, sensor)
    return _scalar_laser_field(laser, sensor) ## it will choose take care of PlaneWave3D ou Gaussian3D
end


## SCATTERING FIELD
function scattering_far_field(problem::NonLinearProblem{MeanField}, β, sensor)
    ## TO DO
    return nothing
end
function scattering_near_field(problem::NonLinearProblem{MeanField}, β, sensor)
    ## TO DO
    return nothing
end
