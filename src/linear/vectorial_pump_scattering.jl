# function _get_intensity(problem::LinearOptics{Vectorial}, field::AbstractVector)
#     return sum(abs2, field)
# end
# function _get_intensity(problem::LinearOptics{Vectorial}, fields::AbstractMatrix)
#     return vec(sum(abs2, fields; dims=1))
# end



## INTENSITY
# function scattered_intensity(problem::LinearOptics{Vectorial}, atomic_states, sensor_positions; regime=:far_field)
#     fields = scattered_electric_field(problem, atomic_states, sensor_positions; regime = regime)
#     intesities = map(eachcol(fields)) do field
#                     # sum(abs2, field) # == mapreduce(abs2, +, field) == norm(field)^2
#                     _get_intensity(problem, field)
#                  end
#     return intesities
# end