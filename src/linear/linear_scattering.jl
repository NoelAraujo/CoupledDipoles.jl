"""
    get_intensity_over_an_angle(problem::LinearOptics{Scalar}, atoms_states::Vector{ComplexF64}, θ::Number; tol=exp10(-7.4))

Used for the single angle and single single state (most probably user case).

# Example:

```julia
using CoupledDipoles, Random
Random.seed!(111)
N = 5
kR, kh = 1.0, 1.0
atoms = Atom(Cylinder(), N, kR, kh)

s, Δ = 1e-5, 1.0
laser = Laser(PlaneWave3D(), s, Δ; polarization=[1,0,0])

problem_scalar = LinearOptics(Scalar(), atoms, laser)
atomic_states_scalar = steady_state(problem_scalar)

θ = deg2rad(48)
get_intensity_over_an_angle(problem_scalar, atomic_states_scalar, θ)
```
"""
function get_intensity_over_an_angle(problem::LinearOptics{Scalar}, atoms_states::Vector{ComplexF64}, angle::Number; tol=exp10(-7.4), exact_solution=false)
    if exact_solution
        xⱼₘ2, yⱼₘ2, zⱼₘ = _rⱼₘ_distances(problem)
        return _intensity_angle_exact_parallel(problem, atoms_states, angle, xⱼₘ2, yⱼₘ2, zⱼₘ)
    else
        return _intensity_angle_approx_quadradure(problem, atoms_states, angle; tol=tol)
    end
end



"""
    get_intensity_over_an_angle(problem::LinearOptics{Scalar}, atoms_states::Vector{Vector{ComplexF64}}, θ::Number; tol=exp10(-7.4), exact_solution=false)

Used for the single angle and different states (for example, the output of `time_evolution`).

# Example:
```julia
using CoupledDipoles, Random
Random.seed!(111)
N = 5
kR, kh = 1.0, 1.0
atoms = Atom(Cylinder(), N, kR, kh)

s, Δ = 1e-5, 1.0
laser = Laser(PlaneWave3D(), s, Δ; polarization=[1,0,0])

problem_scalar = LinearOptics(Scalar(), atoms, laser)
u0 = default_initial_condition(problem_scalar)
tspan = (0.0, 10.0)
solutions = time_evolution(problem_scalar, u0, tspan)
states = solutions.u

θ = deg2rad(48)
get_intensity_over_an_angle(problem_scalar, states, θ)
```
"""
function get_intensity_over_an_angle(problem::LinearOptics{Scalar}, atoms_states::Vector{Vector{ComplexF64}}, angle::Number; tol=exp10(-7.4), exact_solution=false)
    if exact_solution
        # i avoid creating these matrices many times creating them here, and making the same program
        # as the 'single angle and single single state', defined above
        xⱼₘ2, yⱼₘ2, zⱼₘ = _rⱼₘ_distances(problem)
        return map(atoms_states) do β
            _intensity_angle_exact_parallel(problem, β, angle, xⱼₘ2, yⱼₘ2, zⱼₘ)
        end
    end

    return map(atoms_states) do β # O(N)
        _intensity_angle_approx_quadradure(problem, β, angle; tol=tol)
    end
end


function get_intensity_over_an_angle(problem::LinearOptics{Scalar}, atoms_states::Vector{Vector{ComplexF64}}, angles::AbstractArray; tol=exp10(-7.4), exact_solution=false)
    # Many States + Many Angles
    return reduce(hcat, map(atoms_states) do β # O(N)
        get_intensity_over_an_angle(problem, β, angles; tol = tol, exact_solution=exact_solution)
    end)

end
function get_intensity_over_an_angle(problem::LinearOptics{Scalar}, atoms_state::Vector{ComplexF64}, angles::AbstractArray; tol=exp10(-7.4), exact_solution=false)
    # Single States + Many Angles
    return map(angles) do θ
        get_intensity_over_an_angle(problem, atoms_state, θ; tol=tol, exact_solution=exact_solution)
    end
end



@views function _rⱼₘ_distances(problem)
    N = problem.atoms.N
    r = problem.atoms.r

    xⱼ, yⱼ, zⱼ = r[1,:], r[2,:], r[3,:]

    number_configurations = ((N^2) ÷ 2 - N ÷ 2)

    xⱼₘ2 = Array{Float64}(undef, number_configurations)
    yⱼₘ2, zⱼₘ = similar(xⱼₘ2), similar(xⱼₘ2)
    count = 1
    for m in 1:N
        for n in (m+1):N
            xⱼₘ2[count] = (xⱼ[n] - xⱼ[m])^2
            yⱼₘ2[count] = (yⱼ[n] - yⱼ[m])^2
            zⱼₘ[count] = zⱼ[n] - zⱼ[m]
            count += 1
        end
    end

    return xⱼₘ2, yⱼₘ2, zⱼₘ
end



@views function _intensity_angle_exact_parallel(problem::LinearOptics{Scalar}, atoms_states::AbstractVector, θ::Number, xⱼₘ2, yⱼₘ2, zⱼₘ)
    N = problem.atoms.N
    βₙ = view(atoms_states, 1:N)
    βₘ = conj.(βₙ)

    number_configurations = ((N^2) ÷ 2 - N ÷ 2)

    βₙₘ = Array{ComplexF64}(undef, number_configurations)
    count = 1
    for m in 1:N
        b_m = βₘ[m]
        for n in (m+1):N
            b_n = βₙ[n]
            βₙₘ[count] = b_n*b_m
            count += 1
        end
    end

    k₀sinθ = abs(k₀*sin(θ))
    k₀cosθ = k₀*cos(θ)
    total_intensity = ThreadsX.mapreduce(+, 1:number_configurations) do ii
        βₙₘ[ii]*cis(-zⱼₘ[ii]*k₀cosθ)*Bessels.besselj0(k₀sinθ*sqrt(xⱼₘ2[ii] + yⱼₘ2[ii]))
    end
    total_intensity += sum(abs2, βₙ)/2

    k₀R = k₀*how_far_is_farField(problem)

    return (2π)*Γ^2/(4*k₀R^2)*2real(total_intensity)

end





@views function _intensity_angle_approx_quadradure(problem::LinearOptics{Scalar}, atoms_states::AbstractVector, θ::Number; tol=exp10(-7.4))
    R = how_far_is_farField(problem)

    sensor = Array{Float64}(undef,3)
    (intensity, _e) = hcubature([0.0], [2π], rtol = tol, atol=tol) do ϕ
        sensor[1] = R*sin(θ)*cos(ϕ[1])
        sensor[2] = R*sin(θ)*sin(ϕ[1])
        sensor[3] = R*cos(θ)
        scattered_intensity(problem, atoms_states, sensor; regime=:far_field)
    end
    return  intensity[1]
end