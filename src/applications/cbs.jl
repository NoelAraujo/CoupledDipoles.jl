"""
    CBS_scalar(settings, angles, saveat; r_min=0.0)

Coherent Backscattering Matrix (eachrow = angle; eachcol = time)

# Arguments
- `settings`: named tuple (atomsSettings, laserSettings), used as `Atom(atomsSettings...)`, `Laser(laserSettings...)`
- `angles`: Vector or Number
- `saveat`: the times which to save the atomic states, in switch-off condition

# Example
```julia
using CoupledDipoles, Random
b₀, N = 5, 50
kR = √(2N / b₀)

w₀, s, Δ = 4π, 1e-5, 0.0

saveat = range(0, stop=250, length=25)
θ = range(0, 2π; length=180)

atomsSettings = CoupledDipoles.Sphere(; gaussian=true), N, kR
laserSettings = PlaneWave3D(), s, Δ
settings = (atoms=atomsSettings, laser=laserSettings)


Random.seed!(444)
CBS_scalar(settings,  θ, saveat);
CBS_scalar(settings,  θ, [0.0]);
CBS_scalar(settings,  0.0, saveat);
```
"""
function CBS_scalar(settings, angles, saveat; r_min=0.0)
    atomsSettings, laserSettings = settings.atoms, settings.laser
    atoms = Atom(atomsSettings...; r_min=r_min)
    laser = Laser(laserSettings...)
    problem = LinearOptics(Scalar(), atoms, laser) # need far field regime

    ## time evolution
    u₀ = steady_state(problem)
    turn_laser_off!(problem)
    tspan = extrema(saveat)
    β_over_time = time_evolution(problem, u₀, tspan; saveat = saveat)

    ## intensities
    xₙₘ, yₙₘ, zₙₘ = CoupledDipoles._xyz_distances(problem)
    eachrow_angle_eachcol_time = mapreduce(hcat, β_over_time.u) do βₜ
        CoupledDipoles.get_intensity_over_an_angle_fast(problem, βₜ, angles, xₙₘ, yₙₘ, zₙₘ)
    end

    return eachrow_angle_eachcol_time
end


"""
    Evaluate CBS only at steady state
"""
function CBS_scalar(settings, angles; r_min=0.0)
    atomsSettings, laserSettings = settings.atoms, settings.laser
    atoms = Atom(atomsSettings...; r_min=r_min)
    laser = Laser(laserSettings...)
    problem = LinearOptics(Scalar(), atoms, laser) # need far field regime
    u₀ = steady_state(problem)

    xₙₘ, yₙₘ, zₙₘ = CoupledDipoles._xyz_distances(problem)
    intensity_each_angle =  CoupledDipoles.get_intensity_over_an_angle_fast(problem, u₀, angles, xₙₘ, yₙₘ, zₙₘ)

    return intensity_each_angle
end






@views function get_intensity_over_an_angle_fast(problem::LinearOptics{Scalar}, atoms_states::AbstractVector, θ::Number, xₙₘ, yₙₘ, zₙₘ)
    N = problem.atoms.N
    r = problem.atoms.r
    βₙ = view(atoms_states, 1:N)
    βₘ = conj.(βₙ)

    k₀sinθ = abs(k₀*sin(θ))
    k₀cosθ = k₀*cos(θ)
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

    total_intensity = ThreadsX.mapreduce(+, 1:number_configurations) do ii
        βₙₘ[ii]*exp(-im*zₙₘ[ii]*k₀cosθ)*Bessels.besselj0(k₀sinθ*sqrt(xₙₘ[ii]^2+yₙₘ[ii]^2))
    end
    total_intensity += sum(abs2, βₙ)/2
    return 2real(total_intensity)

end
@views function get_intensity_over_an_angle_fast(problem::LinearOptics{Scalar}, atoms_states::AbstractVector, θ_range::AbstractVector, xₙₘ, yₙₘ, zₙₘ)
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

    intensities = zeros(length(θ_range))
    for (idx_θ, θ) in enumerate(θ_range)
        k₀sinθ = abs(k₀*sin(θ))
        k₀cosθ = k₀*cos(θ)

        _intensity = ThreadsX.mapreduce(+, 1:number_configurations) do ii
            βₙₘ[ii]*cis(-zₙₘ[ii]*k₀cosθ)*Bessels.besselj0(k₀sinθ*sqrt(xₙₘ[ii]^2+yₙₘ[ii]^2))
        end
        _intensity += sum(abs2, βₙ)/2

        intensities[idx_θ] = 2real(_intensity)
    end

    return intensities
end



@views function _xyz_distances(problem)
    N = problem.atoms.N
    r = problem.atoms.r
    xₙ, yₙ, zₙ = r[1,:], r[2,:], r[3,:]

    number_configurations = ((N^2) ÷ 2 - N ÷ 2)

    xₙₘ = Array{Float64}(undef, number_configurations)
    yₙₘ, zₙₘ = similar(xₙₘ), similar(xₙₘ)
    count = 1
    for m in 1:N
        for n in (m+1):N
            xₙₘ[count] = xₙ[n] - xₙ[m]
            yₙₘ[count] = yₙ[n] - yₙ[m]
            zₙₘ[count] = zₙ[n] - zₙ[m]
            count += 1
        end
    end
    return xₙₘ, yₙₘ, zₙₘ
end