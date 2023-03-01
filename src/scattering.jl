"""
    scattered_electric_field(problem, atomic_states, sensor_positions; regime=:far_field)

Returns a Matrix{ComplexF64} with value of the `Eletric Laser` + `Electric Scattered` from atoms

- problem: `LinearOptics` or `NonLinearOptics`
- atomic_states: β for `Scalar`/`Vectorial` Model, or [β,z] for `Mean Field` Model
- sensor_positions: matrix with measurement points

Note:
- `Scalar` problem returns a Matrix and not a Vector, to maintain consistency
    with `Vectorial` problem that necessary returns a Matrix,
    where each column has the [Ex, Ey, Ez] components of the field.

- Also, even for single sensor, returns a Matrix of one element.

# Example

```julia
using CoupledDipoles, Random
Random.seed!(111)
N = 5
kR, kh = 1.0, 1.0
atoms = Atom(Cylinder(), N, kR, kh)

s, Δ = 1e-5, 1.0
laser = Laser(PlaneWave3D(), s, Δ; polarization=[1,0,0])

problem_scalar = LinearOptics(Scalar(), atoms, laser)
problem_vectorial = LinearOptics(Vectorial(), atoms, laser)

atomic_states_scalar = steady_state(problem_scalar)
atomic_states_vectorial = steady_state(problem_vectorial)

## 1 sensor
Random.seed!(222)
nSensors = 1
sensor = Matrix(rand(3, nSensors)) # '3' == sensor in position in 3D space
scattered_electric_field(problem_scalar, atomic_states_scalar, sensor)
scattered_electric_field(problem_vectorial, atomic_states_vectorial, sensor)


## 10 sensors
Random.seed!(333)
nSensors = 10
sensor = rand(3, nSensors) # '3' == sensor in position in 3D space
scattered_electric_field(problem_scalar, atomic_states_scalar, sensor)
scattered_electric_field(problem_vectorial, atomic_states_vectorial, sensor)
```
"""
function scattered_electric_field(problem, atomic_states, sensor_positions; regime=:far_field)
    states = prepare_input_state(problem, atomic_states)
    measurement_positions = prepare_input_position(problem, sensor_positions)
    if regime == :far_field
        single_point_field = scattering_far_field
    elseif regime == :near_field
        single_point_field = scattering_near_field
    else
        @error "the regime $(regime) does not exist. The options are ':far_field' or ':near_field'"
    end
    if problem.atoms.N ≤ 25 # SEQUENCIAL map worths for atoms below ~25
        _electric_fields = map(eachcol(measurement_positions)) do sensor
            single_point_field(problem, states, sensor)
        end
    else
        _electric_fields = ThreadsX.map(eachcol(measurement_positions)) do sensor
            single_point_field(problem, states, sensor)
        end
    end
    electric_fields::Matrix{ComplexF64} = hcat(_electric_fields...)

    return electric_fields

end

function laser_and_scattered_electric_field(problem, atomic_states, sensor_positions; regime=:far_field)
    laserField = laser_field(problem, sensor_positions)
    scatteredField = scattered_electric_field(problem, atomic_states, sensor_positions; regime = regime)
    total_field = laserField + scatteredField
    return total_field
end






function prepare_input_state(problem::LinearOptics{Scalar}, atomic_states::AbstractVecOrMat)
    n_states = length(atomic_states)
    if     n_states == problem.atoms.N# if is a vector, return it as a view to reduce memory access
        return view(atomic_states, :, :)
    else
        @error "Atomic States are Invalid. Use an Array with the same length as the Number of Atoms."
    end
end
function prepare_input_state(problem::LinearOptics{Vectorial}, atomic_states::AbstractVecOrMat)
    ## TO DO
    return atomic_states
end
function prepare_input_state(problem::NonLinearOptics{MeanField}, atomic_states::AbstractVector)
    ## TO DO
    return atomic_states
end



function prepare_input_position(problem, sensor_positions::AbstractVector)
    n_axes_components = length(sensor_positions)
    system_dimension = get_dimension(problem.atoms)
    if     n_axes_components == system_dimension
        return view( Array(Matrix(sensor_positions')'), :, :) ## convert array into a matrix of single coluumn
    else
        @error "Your `Measurement Position Vector`  does not match the dimensionality of the system."
    end
end
function prepare_input_position(problem, sensor_positions::AbstractMatrix)
    dimensions, n_sensors = size(sensor_positions)
    system_dimension = get_dimension(problem.atoms)
    if     dimensions == system_dimension
        return view(sensor_positions, :, :)
    else
        @error "Your `Measurement Position Matrix`  does not match the dimensionality of the system."
    end
end




## INTENSITY
"""
    scattered_intensity(problem, atomic_states, sensor_positions; regime=:far_field)

Returns a Vector{Float64} with value of the `|Electric Scattered|^2` from atoms

- problem: `LinearOptics` or `NonLinearOptics`
- atomic_states: β for `Scalar`/`Vectorial` Model, or [β,z] for `Mean Field` Model
- sensor_positions: matrix with measurement points

# Example

```julia
using CoupledDipoles, Random
Random.seed!(111)
N = 5
kR, kh = 1.0, 1.0
atoms = Atom(Cylinder(), N, kR, kh)

s, Δ = 1e-5, 1.0
laser = Laser(PlaneWave3D(), s, Δ; polarization=[1,0,0])

problem_scalar = LinearOptics(Scalar(), atoms, laser)
problem_vectorial = LinearOptics(Vectorial(), atoms, laser)

atomic_states_scalar = steady_state(problem_scalar)
atomic_states_vectorial = steady_state(problem_vectorial)

## 1 sensor
Random.seed!(222)
nSensors = 1
sensor = Matrix(rand(3, nSensors)) # '3' == sensor in position in 3D space
scattered_intensity(problem_scalar, atomic_states_scalar, sensor)
scattered_intensity(problem_vectorial, atomic_states_vectorial, sensor)


## 10 sensors
Random.seed!(333)
nSensors = 10
sensor = rand(3, nSensors) # '3' == sensor in position in 3D space
scattered_intensity(problem_scalar, atomic_states_scalar, sensor)
scattered_intensity(problem_vectorial, atomic_states_vectorial, sensor)
```
"""
function scattered_intensity(problem, atomic_states, sensor_positions; regime=:far_field)
    fields = scattered_electric_field(problem, atomic_states, sensor_positions; regime = regime)
    intesities = _get_intensity(problem, fields)
    return intesities
end

"""
    laser_and_scattered_intensity(problem, atomic_states, sensor_positions; regime=:far_field)

Returns a Vector{Float64} with value of the `|Electric Laser + Electric Scattered|^2` from atoms

- problem: `LinearOptics` or `NonLinearOptics`
- atomic_states: β for `Scalar`/`Vectorial` Model, or [β,z] for `Mean Field` Model
- sensor_positions: matrix with measurement points
"""
function laser_and_scattered_intensity(problem, atomic_states, sensor_positions; regime=:far_field)
    total_field = laser_and_scattered_electric_field(problem, atomic_states, sensor_positions; regime = regime)
    intesities = _get_intensity(problem, total_field)
    return intesities
end



function scattered_intensity(problem::NonLinearOptics{MeanField}, atomic_states, sensor_positions; regime=:far_field)
    fields = scattered_electric_field(problem, atomic_states, sensor_positions; regime = regime)

    if regime == :far_field
        intensities = ThreadsX.map(pairs(eachcol(sensor_positions))) do sensor_field
            field, sensor = fields[sensor_field[1]], sensor_field[2]
            _get_intensity_far_field(problem, field, atomic_states, norm(sensor))
        end
    elseif regime == :near_field
        r = problem.atoms.r
        intensities = ThreadsX.map(pairs(eachcol(sensor_positions))) do sensor_field
            field, sensor = fields[sensor_field[1]], sensor_field[2]
            _get_intensity_near_field(problem, field, atomic_states, sensor, r)
        end
    else
        @error "the regime $(regime) does not exist. The options are ':far_field' or ':near_field'"
    end

    return intensities
end

function laser_and_scattered_intensity(problem::NonLinearOptics{MeanField}, atomic_states, sensor_positions; regime=:far_field)
    total_field = laser_and_scattered_electric_field(problem, atomic_states, sensor_positions; regime = regime)

    if regime == :far_field
        intensities = ThreadsX.map(pairs(eachcol(sensor_positions))) do sensor_field
            field, sensor = total_field[sensor_field[1]], sensor_field[2]
            _get_intensity_far_field(problem, field, atomic_states, norm(sensor))
        end
    elseif regime == :near_field
        r = problem.atoms.r
        intensities = ThreadsX.map(pairs(eachcol(sensor_positions))) do sensor_field
            field, sensor = total_field[sensor_field[1]], sensor_field[2]
            _get_intensity_near_field(problem, field, atomic_states, sensor, r)
        end
    else
        @error "the regime $(regime) does not exist. The options are ':far_field' or ':near_field'"
    end
    return intensities
end