function _get_intensity_near_field(problem::NonLinearOptics{PairCorrelation}, atomic_states, sensor; inelasticPart=true)
    N = problem.atoms.N
    r = problem.atoms.r

    σᶻ = atomic_states[N+1:2N]
    σ⁻σ⁺ = reshape(atomic_states[2*N+1+N^2:2*N+2*N^2], (N, N))
    σ⁺σ⁻ = transpose(σ⁻σ⁺)

    # atomsPhases =[ j≠m ? cis( dot(sensor,r[:,j]-r[:,m]) )/( norm(sensor-r[:,j])*norm(sensor-r[:,m])  ) : 0im for j=1:N, m=1:N]
    # I_pairCorrelation = sum(atomsPhases.*σ⁺σ⁻) # NOT A MATRIX MULTIPLICATION, IT IS A ELEMENT-WISE MULTIPLICATION

    I_pairCorrelation = sum(σ⁺σ⁻[j,m]*cis( dot(sensor,r[:,j]-r[:,m]) )/( norm(sensor-r[:,j])*norm(sensor-r[:,m])  ) for m=1:N, j=1:N if j≠m);
    if inelasticPart
        I_pairCorrelation = I_pairCorrelation + sum(  0.5*(1 + σᶻ[j])/norm(sensor-r[:,j])^2 for j=1:N)
    end
    return (Γ^2/4)*I_pairCorrelation
end

function scattered_intensity(problem::NonLinearOptics{PairCorrelation}, atomic_states, sensor_positions; regime=:far_field, inelasticPart=true)
    if regime == :far_field
        intensities = ThreadsX.map(eachcol(sensor_positions)) do sensor
            # _get_intensity_far_field(problem,  atomic_states, norm(sensor); inelasticPart=inelasticPart)
            rand() # TO DO
        end
    elseif regime == :near_field
        r = problem.atoms.r
        intensities = ThreadsX.map(eachcol(sensor_positions)) do sensor
            _get_intensity_near_field(problem, atomic_states, sensor; inelasticPart=inelasticPart)
        end
    else
        @error "the regime $(regime) does not exist. The options are ':far_field' or ':near_field'"
    end

    return intensities
end