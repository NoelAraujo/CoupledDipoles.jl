# function scattering_fuction(distance::Symbol, dimension::Symbol)
#     if dimension == :ThreeD
#         if distance == :nearField
#             sf = default_farField3D_Field
#         elseif distance == :farField
#             sf = default_farField3D_Field
#         else
#             @error "Only options are `:nearField` and `:farField`"
#         end
#     else
#         @error "Only :ThreeD is supported."
#     end
#     return sf
# end

# """
#     default_nearField3D_Field(atoms::AbstractMatrix, β::AbstractArray, sensor::AbstractArray)

# - atoms: each column is an atom
# - atomic_states: β for Scalar Model, ou [β,z] for Mean Field Model
# - sensor: measurement position

# Returns a Complex Value of the Eletric Field: +(Γ/2) * ∑ⱼ exp(-i*k₀* n̂⋅R⃗ⱼ)/(k₀* sensor⋅R⃗ⱼ)
# - R⃗ⱼ : atom j
# """
# function default_nearField3D_Field(atoms::AbstractMatrix, β::AbstractArray, sensor::AbstractArray)
#     E_scatt = zero(eltype(β))

#     j = 1
#     @inbounds for atom in eachcol(atoms)
#         d_SensorAtom = sqrt((sensor[1] - atom[1])^2 + (sensor[2] - atom[2])^2 + (sensor[3] - atom[3])^2)
#         E_scatt += cis(k₀ * d_SensorAtom) * (β[j] / d_SensorAtom)
#         j += 1
#     end
#     E_scatt = +(Γ / 2) * im * E_scatt
#     return E_scatt
# end

# """
#     default_farField3D_Field(atoms::AbstractMatrix, β::AbstractArray, sensor::AbstractArray)

# - atoms: each column is an atom
# - atomic_states: β for Scalar Model, ou [β,z] for Mean Field Model
# - sensor: measurement position

# Returns a Complex Value of the Eletric Field: -(Γ/2) * (exp(ikr) / ikr) * ∑ⱼ exp(-i*k₀* n̂⋅R⃗ⱼ)

# - r : distance sensor to origin (r = norm(sensor))
# - n̂ : norm of sensor ( n̂ = sensor / norm(sensor) )
# - R⃗ⱼ : atom j
# """
# function default_farField3D_Field(atoms::AbstractMatrix, β::AbstractArray, sensor::AbstractArray)
#     E_scatt = zero(eltype(β))
#     n̂ = sensor / norm(sensor)

#     dot_n_r = zero(eltype(β))
#     j = 1
#     @inbounds for atom in eachcol(atoms)
#         dot_n_r = n̂[1] * atom[1] + n̂[2] * atom[2] + n̂[3] * atom[3]
#         dot_n_r = cis(-k₀ * dot_n_r)
#         E_scatt += dot_n_r * β[j]
#         j += 1
#     end

#     ikr = im * k₀ * norm(sensor)
#     E_scatt = -(Γ / 2) * E_scatt * exp(ikr) / ikr
#     return E_scatt
# end

# """
#     scattering_intensity(problem, atomic_states, measurement_positions, scattering_func::Function)

#     Computes the Total Electric Field (Laser Pump + Scattering), then returns the Intensity
# """
# @views function scattering_intensity(problem, atomic_states, measurement_positions, scattering_func::Function; useThreads=true)
#     _laser = problem.laser
#     _r = problem.atoms.r
#     _physics = problem.physic

#     _sensors = measurement_positions
#     _states = atomic_states
#     _func = scattering_func

#     n_sensors = _get_number_elements(_sensors)

#     if n_sensors == 1
#         return _OnePoint_Intensity(_physics, _laser, _r, _sensors, _states, _func)
#     else
#         scat_int = Array{Float64}(undef, n_sensors)
#         if useThreads
#             Threads.@threads for i in 1:n_sensors
#                 scat_int[i] = _OnePoint_Intensity(_physics, _laser, _r, view(_sensors, :, i), _states, _func)
#             end
#         else
#             for i in 1:n_sensors
#                 scat_int[i] = _OnePoint_Intensity(_physics, _laser, _r, view(_sensors, :, i), _states, _func)
#             end
#         end
#         return scat_int
#     end
# end





"""
    scattered_electric_field(problem, atomic_states, sensor_positions; regime=:far_field)

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

    _electric_fields = ThreadsX.map(eachcol(measurement_positions)) do sensor
        single_point_field(problem, states, sensor)
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
function scattered_intensity(problem, atomic_states, sensor_positions; regime=:far_field)
    fields = scattered_electric_field(problem, atomic_states, sensor_positions; regime = regime)
    intesities = _get_intensity(problem, fields)
    return intesities
end

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