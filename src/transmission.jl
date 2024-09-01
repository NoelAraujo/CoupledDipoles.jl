"""
    transmission(problem, detunning_range::AbstractArray; regime=:near_field)

One point is created at `far_field_distance = 100*size(problem.atoms)` in any direction.
Then, transmission (Intensity Scattered over Intensity Laser) is evaluated with on steady state (computed automatically) for different detunnings.
"""
function transmission(problem, detunning_range::AbstractArray; regime=:near_field)
    sensor = getSensor_at_front(problem)

    # I need a dummy problem just to compute the electric field intensity
    # necessary for normalization
    I_laser = let         
        _N, _R, _L = 3, 3, 3
        _atoms = Atom(Cylinder(), _N, _R, _L)
        
        _Δ = 0.0 # the detunning does not matter for the laser field
        # BUT, the saturation DOES matter
        _laser = similar_laser(problem.laser, problem.laser.s, _Δ)
        
        ## dummy problem
        _problem = similar_problem(problem, _atoms, _laser)
        
        laser_intensity(_problem, sensor)[1]
    end

    Transmissions = map(detunning_range) do Δ       
        _laser = similar_laser(problem.laser, problem.laser.s, Δ)        
        _problem = similar_problem(problem, problem.atoms, _laser)
        beta = steady_state(_problem)

        _temp = laser_and_scattered_intensity(_problem, beta, sensor; regime=regime)[1]
        _temp / I_laser
    end
    Transmissions
end

@inline function getSensor_at_front(problem)
    new_R = 100*size(problem.atoms) # near field values have less errors
    sensor = reshape(new_R*problem.laser.direction, :, 1)
    return sensor
end


function similar_laser(laser::Laser{PlaneWave3D}, s, Δ)
    direction = laser.direction

    Laser(PlaneWave3D(), s, Δ;         
            direction = direction) 
end
function similar_laser(laser::Laser{Gaussian3D}, s, Δ)
    w₀ = laser.pump.w₀
    polarization = laser.polarization
    direction = laser.direction

    Laser(Gaussian3D(w₀), s, Δ; 
            polarization = polarization,
            direction = direction) 
end

function similar_problem(problem::LinearOptics{Scalar}, atoms, laser)
    LinearOptics(Scalar(), atoms, laser)
end
function similar_problem(problem::LinearOptics{Vectorial}, atoms, laser)
    LinearOptics(Vectorial(), atoms, laser)
end
function similar_problem(problem::NonLinearOptics{MeanField}, atoms, laser)
    NonLinearOptics(MeanField(), atoms, laser)
end
function similar_problem(problem::NonLinearOptics{PairCorrelation}, atoms, laser)
    NonLinearOptics(PairCorrelation(), atoms, laser)
end
