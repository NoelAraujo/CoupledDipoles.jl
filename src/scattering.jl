
## INTENSITY
"""
    scattered_intensity(problem, atomic_states, sensor_positions; regime=:far_field, use_sequencial=false)

Returns a Vector{Float64} with value of the `|Electric Scattered|^2` from atoms

- problem: `LinearOptics` or `NonLinearOptics`
- atomic_states: β for `Scalar`/`Vectorial` Model, or [β,z] for `Mean Field` Model
- sensor_positions: matrix with measurement points
- `use_sequencial` turn on/off the internal parallelism

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
function scattered_intensity(problem::LinearOptics{T}, atomic_states, sensor_positions; regime=:far_field, use_sequencial=false) where {T<:Union{Scalar,Vectorial}}
    fields = scattered_electric_field(problem, atomic_states, sensor_positions; regime=regime, use_sequencial=use_sequencial)

    intesities = mapreduce(vcat, eachcol(fields)) do field
        _get_intensity(problem, field)
    end
    return intesities
end

function _get_intensity(problem::LinearOptics{Scalar}, field::AbstractVector)
    return abs2.(field)
end

function _get_intensity(problem::LinearOptics{Vectorial}, field::AbstractVector)
    return vec([sum(abs2, field)])
end

function _get_intensity(problem::NonLinearOptics{MeanField}, field::AbstractVector)
    return abs2.(field)
end
function _get_intensity(problem::NonLinearOptics{PairCorrelation}, field::AbstractVector)
    return abs2.(field)
end

function scattered_intensity(problem::NonLinearOptics{T}, atomic_states, sensors; regime=:far_field, inelasticPart=true, use_sequencial=false) where {T<:Union{MeanField,PairCorrelation}}

    ## define the 'function' to be used
    if regime == :far_field
        get_intensity = _get_intensity_far_field
    elseif regime == :near_field
        get_intensity = _get_intensity_near_field
    else
        @error "the regime $(regime) does not exist. The options are ':far_field' or ':near_field'"
    end

    ## it is more efficient to compute the Scalar intensity here, and pass as arguments
    ## right now: this is usefull for MeanField, but not on PairCorrelation

    _dummy = LinearOptics(Scalar(), problem.atoms, problem.laser)
    _atomic_states = view(atomic_states, 1:problem.atoms.N)
    intensities_scalar = scattered_intensity(_dummy, _atomic_states, sensors)

    intensities = get_intensity(problem, atomic_states, sensors, intensities_scalar; inelasticPart=inelasticPart)
    return intensities
end





## mean field
function _get_intensity_far_field(problem::NonLinearOptics{MeanField}, atomic_states::AbstractVector, sensors, _intensities; inelasticPart=true)
    N = problem.atoms.N
    σ⁻ = view(atomic_states, 1:N)
    σᶻ = view(atomic_states, (N+1):2N)
    r = how_far_is_farField(problem)

    I_meanField = _intensities
    if inelasticPart
        inelast_factor = (Γ^2 / (4*(k₀^2) * r^2)) * sum(-abs2.(σ⁻) .+ (1 .+ σᶻ) ./ 2)
        I_meanField = I_meanField .+ inelast_factor
    end
    return real.(I_meanField)
end
function _get_intensity_near_field(problem::NonLinearOptics{MeanField}, atomic_states::AbstractVector, sensors, _intensities; inelasticPart=true)
    N = problem.atoms.N
    σ⁻ = view(atomic_states, 1:N)
    σᶻ = view(atomic_states, (N+1):2N)
    r = view(problem.atoms.r, :, :)

    I_meanField = _intensities
    if inelasticPart
        inelast_factor = mapreduce(vcat, eachcol(sensors)) do sensor
            sum(abs2(σ⁻[j]) / norm(sensor - r[:, j])^2 + (1 + σᶻ[j]) / (2 * norm(sensor - r[:, j])^2) for j = 1:N)
        end
        I_meanField = I_meanField .+ (Γ^2 / (4*(k₀^2))) * inelast_factor
    end
    return real.(I_meanField)
end

## pair correlation
function _get_intensity_far_field(problem::NonLinearOptics{PairCorrelation}, atomic_states::AbstractVector, sensors, _intensities; inelasticPart=true)
    N = problem.atoms.N
    r = view(problem.atoms.r, :, :)
    σᶻ = view(atomic_states, (N+1):2N)
    σ⁻σ⁺ = reshape(atomic_states[2*N+1+N^2:2*N+2*N^2], (N, N))
    σ⁺σ⁻ = transpose(σ⁻σ⁺)
    r_far_field = how_far_is_farField(problem)

    I_pairCorrelation = mapreduce(vcat, eachcol(sensors)) do sensor
        sensor_norm = k₀ * norm(sensor)
        n̂ = sensor / sensor_norm
        _I_pairCorrelation = sum(σ⁺σ⁻[j, m] * cis(dot(n̂, r[:, j] - r[:, m])) for m = 1:N, j = 1:N if j ≠ m)
        if inelasticPart
            inelast_factor = sum((1 .+ σᶻ) ./ 2)
            _I_pairCorrelation = _I_pairCorrelation .+ inelast_factor
        end
        _I_pairCorrelation = (Γ^2 / (4 * (k₀^2) * r_far_field^2)) * _I_pairCorrelation
    end

    nSensors = size(sensors,2)
    if nSensors == 1
        return [real.(I_pairCorrelation)]
    else
        return real.(I_pairCorrelation)
    end
end
function _get_intensity_near_field(problem::NonLinearOptics{PairCorrelation}, atomic_states::AbstractVector, sensors, _intensities; inelasticPart=true)
    N = problem.atoms.N
    r = view(problem.atoms.r, :, :)

    σᶻ = view(atomic_states, (N+1):2N)
    σ⁺σ⁻ = reshape(atomic_states[2*N+1+N^2:2*N+2*N^2], (N, N))

    # atomsPhases =[ j≠m ? cis( dot(sensor,r[:,j]-r[:,m]) )/( norm(sensor-r[:,j])*norm(sensor-r[:,m])  ) : 0im for j=1:N, m=1:N]
    # I_pairCorrelation = sum(atomsPhases.*σ⁺σ⁻) # NOT A MATRIX MULTIPLICATION, IT IS A ELEMENT-WISE MULTIPLICATION

    I_pairCorrelation = mapreduce(vcat, eachcol(sensors)) do sensor
        _I_pairCorrelation = sum(
            σ⁺σ⁻[j,m]*
            cis( norm(sensor - r[:,j]) - norm(sensor - r[:,m])  ) / (norm(sensor-r[:,j])*norm(sensor-r[:,m]))
            for j=1:N, m=1:N  if j ≠ m)
        if inelasticPart
            _I_pairCorrelation = _I_pairCorrelation + sum(0.5 * (1 + σᶻ[j]) / norm(sensor - r[:, j])^2 for j = 1:N)
        end
        real(_I_pairCorrelation)
    end

    nSensors = size(sensors,2)
    if nSensors == 1
        return [(Γ^2 / (4 * k₀^2)) * I_pairCorrelation]
    else
        return (Γ^2 / (4 * k₀^2)) * I_pairCorrelation
    end
end









# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
"""
    laser_and_scattered_intensity(problem, atomic_states, sensor_positions; regime=:far_field)

Returns a Vector{Float64} with value of the `|Electric Laser + Electric Scattered|^2` from atoms

- problem: `LinearOptics` or `NonLinearOptics`
- atomic_states: β for `Scalar`/`Vectorial` Model, or [β,z] for `Mean Field` Model
- sensors: matrix with measurement points
"""
function laser_and_scattered_intensity(problem, atomic_states, sensors; regime=:far_field)
    E_laser = laser_field(problem, sensors)
    E_scattered = scattered_electric_field(problem, atomic_states, sensors;regime=regime)

    I_laser = laser_intensity(problem, sensors)
    I_scattered = scattered_intensity(problem, atomic_states, sensors; regime=regime)

    total_field = I_laser + I_scattered + 2*real.(_mul_fields(problem, E_laser, E_scattered))

    return total_field
end

function _mul_fields(problem::LinearOptics{Scalar}, E_laser, E_scattered)
    vec(E_laser.*E_scattered)
end
function _mul_fields(problem::LinearOptics{Vectorial}, E_laser, E_scattered)
    vec(sum(E_laser.*E_scattered, dims=1))
end
function _mul_fields(problem::NonLinearOptics{MeanField}, E_laser, E_scattered)
    vec(E_laser.*E_scattered)
end
function _mul_fields(problem::NonLinearOptics{PairCorrelation}, E_laser, E_scattered)
    vec(E_laser.*E_scattered)
end