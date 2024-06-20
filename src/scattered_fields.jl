# --------------------------------- GENERAL FUNCTIONS ---------------------------------
"""
    scattered_electric_field(problem, atomic_states, sensors; regime=:near_field, use_sequencial=false)

Returns a Matrix{ComplexF64} with value of the `Eletric Laser` + `Electric Scattered` from atoms

- problem: `LinearOptics` or `NonLinearOptics`
- atomic_states: A vector for `Scalar`, `MeanField` and `PairCorrelation`. A matrix `Vectorial` Model
- sensors: matrix with measurement points
- `use_sequencial` turn on/off the internal parallelism

Note:
- `Scalar` (and other models) problem returns a **Matrix and NOT a Vector**, to maintain consistency
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
sensor = Matrix(rand(3, nSensors)) # '3' == sensor lives in a 3D space
scattered_electric_field(problem_scalar, atomic_states_scalar, sensor)
scattered_electric_field(problem_vectorial, atomic_states_vectorial, sensor)


## 10 sensors
Random.seed!(333)
nSensors = 10
sensor = rand(3, nSensors) # '3' == sensor lives in a 3D space
scattered_electric_field(problem_scalar, atomic_states_scalar, sensor)
scattered_electric_field(problem_vectorial, atomic_states_vectorial, sensor)
```
"""
function scattered_electric_field(problem, atomic_states, sensors; regime=:near_field, use_sequencial=false)
    ## validate inputs
    states = prepare_states(problem, atomic_states)
    measurement_positions = prepare_sensors(problem, sensors)

    ## define the 'function' to be used
    if regime == :far_field
        single_point_field = scattering_far_field
    elseif regime == :near_field
        single_point_field = scattering_near_field
    else
        @error "the regime $(regime) does not exist. The options are ':far_field' or ':near_field'"
    end

    ## compute fields
    if problem.atoms.N ≤ 25 || use_sequencial==true # SEQUENCIAL map is faster for atoms below ~25
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



function prepare_sensors(problem, sensors::AbstractMatrix)
    data_dimension, n_sensors = size(sensors)
    system_dimension = get_dimension(problem.atoms)
    if  data_dimension == system_dimension
        return view(sensors, :, :)
    else
        @error "Your `Measurement Position Matrix`  does not match the dimensionality of the system."
    end
end
function prepare_sensors(problem, sensors::AbstractVector)
    @error "Your `Measurement Positions` are not Valid. You must specify a Matrix and not a Vector. \n Maybe `Array(Matrix(sensors')')` is a solution."
end


function laser_and_scattered_electric_field(problem, atomic_states, sensors; regime = :near_field)
    E_laser = laser_field(problem, sensors)
    E_scattered = scattered_electric_field(problem, atomic_states, sensors; regime=regime)
    E_total = E_laser + E_scattered
    return E_total
end




# --------------------------------- MODEL SPECIFIC ---------------------------------
function prepare_states(problem::LinearOptics{Scalar}, atomic_states)
    N = length(atomic_states)
    if N == problem.atoms.N# if is a vector, return it as a view to reduce memory access
        σ⁻ = view(atomic_states, :)
        return σ⁻
    else
        @error "Atomic States are Invalid. Use an Array with the same length as the Number of Atoms."
    end
end
function prepare_states(problem::LinearOptics{Vectorial}, atomic_states)
    N = size(atomic_states, 2)
    if N == problem.atoms.N# if is a matrix, return it as a view to reduce memory access
        σ⁻ = view(atomic_states, :, :)
        return σ⁻
    else
        @error "Atomic States are Invalid. Use a Matrix with polarization components at each line."
    end
end
function prepare_states(problem::NonLinearOptics{MeanField}, atomic_states)
    expected_size = 2*problem.atoms.N
    if length(atomic_states) == expected_size
        σ⁻ = view(atomic_states, 1:problem.atoms.N)
        return σ⁻
    else
        @error "Atomic States are Invalid. Use an Array with the same length as the Number of Atoms."
    end
end
function prepare_states(problem::NonLinearOptics{PairCorrelation}, atomic_states)
    expected_size = 2*problem.atoms.N + 4*(problem.atoms.N^2)
    if length(atomic_states) == expected_size
        σ⁻ = view(atomic_states, 1:problem.atoms.N)
        return σ⁻
    else
        @error "Atomic States are Invalid. Use an Array with the same length as the Number of Atoms."
    end
end






#=
            SCATTERED FIELD: :near_field
=#
function scattering_near_field(problem::LinearOptics{Scalar}, β, sensor)
    _scalar_scattering_near_field(problem.atoms, β, sensor)
end
function scattering_near_field(problem::LinearOptics{Vectorial}, β, sensor)
    _vectorial_scattering_near_field(problem.atoms, β, sensor)
end

## The code below is NOT a mistake (so far). The formulas for electric field depends only on sigma^-,
## and have the same formula as the Scalar case
function scattering_near_field(problem::NonLinearOptics{MeanField}, β, sensor)
    _scalar_scattering_near_field(problem.atoms, β, sensor)
end
function scattering_near_field(problem::NonLinearOptics{PairCorrelation}, β, sensor)
    _scalar_scattering_near_field(problem.atoms, β, sensor)
end




#=
            SCATTERED FIELD: :far_field
=#
function scattering_far_field(problem::LinearOptics{Scalar}, β, sensor)
    _scalar_scattering_far_field(problem.atoms, β, sensor)
end
function scattering_far_field(problem::LinearOptics{Vectorial}, β, sensor)
    _vectorial_scattering_far_field(problem.atoms, β, sensor)
end
function scattering_far_field(problem::NonLinearOptics{MeanField}, β, sensor)
    _scalar_scattering_near_field(problem.atoms, β, sensor)
end
function scattering_far_field(problem::NonLinearOptics{PairCorrelation}, β, sensor)
    _scalar_scattering_near_field(problem.atoms, β, sensor)
end