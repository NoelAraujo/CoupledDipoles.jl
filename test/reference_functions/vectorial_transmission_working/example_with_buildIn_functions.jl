using LinearAlgebra, Random, CoupledDipoles

Δ_range = range(-10, 10; length=16)
Transmissions = zeros(length(Δ_range))

nRepetitions = 2

N = 500
R = 17.71830997686217301634315
L = 8.177681527782540982229875
s, Δ = 1e-3, -10.0
w₀ =  R/3
new_R = 100*R

sensor, E_laser, I_laser = let 
    atoms = Atom(Cylinder(), N, R, L)
    
    _Δ = 0.0 # the detunning does not matter for the laser field
    # BUT, the saturation DOES matter
    laser = Laser(Gaussian3D(w₀), s, _Δ; polarization=[1,0,0]) 

    problem = LinearOptics(Vectorial(), atoms, laser)
    sensor = reshape(new_R*problem.laser.direction, :, 1)
    
    E_laser = laser_field(problem, sensor)
    I_laser = laser_intensity(problem, sensor)[1]
    sensor, E_laser, I_laser
end

for ib in 1:length(Δ_range)
    println(ib)
    _Δ_ = Δ_range[ib]
    
    _temp = mapreduce(+, 1:nRepetitions) do j
        Random.seed!(j)
                   
        atoms = Atom(Cylinder(), N, R, L)
        laser = Laser(Gaussian3D(w₀), s, _Δ_; polarization=[1,0,0])
        problem = LinearOptics(Vectorial(), atoms, laser)
        betass = steady_state(problem)

        laser_and_scattered_intensity(problem, betass, sensor; regime=:near_field)[1]
    end

    Transmissions[ib] = (_temp/nRepetitions) / I_laser
end

extrema(Transmissions)

