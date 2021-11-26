using LinearAlgebra

"""
    Far Fiel Condition: r_emitter^2 / 2r_observation ≪ 1 
    
Condition comes from is the Eq (4.4) in  
    https://engineering.purdue.edu/wcchew/ece604f18/latex%20lecture%20notes/LectureNotes20.pdf  
(A copy of the file is inside the `benchmarks/pdf` folder)  

In practice, we verify if the condition is smaller than `1e-1`, 
because we tested that the estimate error is around 1%.
"""
function test_farField_validity(r_observation, r_emitter)
    if abs( norm(r_emitter)^2/(2*norm(r_observation))) ≤ 1e-1
        return true
    else
        return false
    end
    nothing
end

r_atom = 10*rand(3)
scalings = 10.0.^range(1, 7, length=100)

all_error_real = []
all_error_imag = []
conditions = []
for rep=1:200
    for scaling_sensor = scalings
        r_sensor = scaling_sensor.*rand()*ones(3)

        vector_distance_SensorAtom = r_sensor - r_atom
        distance_SensorAtom = norm(vector_distance_SensorAtom)

        distance_sensor = norm(r_sensor)
        direction_sensor = r_sensor / distance_sensor


        E_exact = -0.5im * cis(distance_SensorAtom) / (im * distance_SensorAtom)
        E_farfield = -0.5im *cis(distance_sensor) * cis(-dot(direction_sensor, r_atom)) / (im * distance_sensor)

        diffE = E_exact - E_farfield
        error_real = 100*abs(real(diffE)/real(E_exact))
        error_imag = 100*abs(imag(diffE)/imag(E_exact))

        push!(all_error_real, error_real)
        push!(all_error_imag, error_imag)

        condition = abs( norm(r_atom)^2/(2*norm(r_sensor)))
        push!(conditions, condition)
    end
end

using Plots

scatter(conditions, all_error_real, label="error real",markershape=:circle,markersize=6, xflip=true)
scatter!(conditions, all_error_imag, label="error imag", markershape=:rect, markersize=5,scale=:log10)
xlabel!("Validity Condition")
ylabel!("Percentage Error")
hline!([1], c=:black, label="1%", linestyle=:dash, lw=2)
vline!([1e-1], c=:red, label="Aceptable Valid Condition", lw=2)

#=
    I state that an error aceptable is around 1% - the black dashed line.
    Therefore, the `Validity Condition` required is around 10^-1.
    Which is exact reason of the `10^-1` used in test_farField_validity() function.

    Now that I have a defined expression, I can figure out the Far Field distance:

    |r_emitter|^2/(2*|r_observation|) ≤ 1e-1
    |r_emitter|^2/(2*1e-1) ≤ |r_observation|
    |r_emitter|^2/(0.2) ≤ |r_observation|
    |r_emitter|^2/(1/5) ≤ |r_observation|
    
    |r_observation| ≥ 5|r_emitter|^2

    Just to be a little bit above this condition, I will define
    a function with "5.01" and not "5.0" 
=#
function estimate_far_Field(r_emitter)
    return 5.01*(norm(r_emitter)^2)
end
