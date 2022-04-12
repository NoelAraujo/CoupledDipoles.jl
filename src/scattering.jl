function scattering_fuction(distance::Symbol, dimension::Symbol)
    if dimension == :ThreeD
        if distance == :nearField
            sf = default_farField3D_Field
        elseif distance == :farField
            sf = default_farField3D_Field
        else
            @error "Only options are `:nearField` and `:farField`"
        end
    else
        @error "Only :ThreeD is supported."
    end
    return sf
end

"""
    default_nearField3D_Field(atoms::AbstractMatrix, β::AbstractArray, sensor::AbstractArray)

- atoms: each column is an atom
- atomic_states: β for Scalar Model, ou [β,z] for Mean Field Model
- sensor: measurement position

Returns a Complex Value of the Eletric Field: +(Γ/2) * ∑ⱼ exp(-i*k₀* n̂⋅R⃗ⱼ)/(k₀* sensor⋅R⃗ⱼ)
- R⃗ⱼ : atom j
"""
function default_nearField3D_Field(atoms::AbstractMatrix, β::AbstractArray, sensor::AbstractArray)
    E_scatt = zero(eltype(β))

    j = 1
    @inbounds for atom in eachcol(atoms)
        d_SensorAtom = sqrt(
            (sensor[1] - atom[1])^2 + (sensor[2] - atom[2])^2 + (sensor[3] - atom[3])^2,
        )
        E_scatt += cis(k₀ * d_SensorAtom) * (β[j] / d_SensorAtom)
        j += 1
    end
    E_scatt = +(Γ / 2) * im * E_scatt
    return E_scatt
end

"""
    default_farField3D_Field(atoms::AbstractMatrix, β::AbstractArray, sensor::AbstractArray)

- atoms: each column is an atom
- atomic_states: β for Scalar Model, ou [β,z] for Mean Field Model
- sensor: measurement position

Returns a Complex Value of the Eletric Field: -(Γ/2) * (exp(ikr) / ikr) * ∑ⱼ exp(-i*k₀* n̂⋅R⃗ⱼ) 
     
- r : distance sensor to origin (r = norm(sensor))
- n̂ : norm of sensor ( n̂ = sensor / norm(sensor) )
- R⃗ⱼ : atom j
"""
function default_farField3D_Field(atoms::AbstractMatrix, β::AbstractArray, sensor::AbstractArray)
    E_scatt = zero(eltype(β))
    n̂ = sensor / norm(sensor)

    dot_n_r = zero(eltype(β))
    j = 1
    @inbounds for atom in eachcol(atoms)
        dot_n_r = n̂[1] * atom[1] + n̂[2] * atom[2] + n̂[3] * atom[3]
        dot_n_r = cis(-k₀ * dot_n_r)
        E_scatt += dot_n_r * β[j]
        j += 1
    end

    ikr = im * k₀ * norm(sensor)
    E_scatt = -(Γ / 2) * E_scatt * exp(ikr) / ikr
    return E_scatt
end

"""
    scattering_intensity(SP::ScatteringProblem)

    Computes the Total Electric Field (Laser Pump + Scattering), then returns the Intensity
"""
@views function scattering_intensity(
    problem,
    atomic_states,
    measurement_positions,
    scattering_func::Function,
)
    _laser = problem.laser
    _r = problem.atoms.r
    _physics = problem.physic
    
    _sensors = measurement_positions
    _states = atomic_states
    _func = scattering_func

    n_sensors = _get_number_elements(_sensors)
    
    # N = size(_r, 2)
    # number_configurations = ((N^2)÷2 - N÷2)
    # βₙₘ = Array{eltype(_states[1])}(undef, number_configurations)

    if n_sensors == 1
        return _OnePoint_Intensity(_physics, _laser, _r, _sensors, _states, _func)
    else
        scat_int = zeros(n_sensors)
        for i = 1:n_sensors #Threads.@threads 
            scat_int[i] = _OnePoint_Intensity(_physics, _laser, _r, _sensors[:, i], _states, _func)
        end
        return scat_int
    end
end

function _OnePoint_Intensity(physic::Union{Scalar,Vectorial}, laser, atoms, sensor, β, scattering_func)
    E_L = laser_field(laser, sensor)
    E_scatt = scattering_func(atoms, β, sensor)
    return abs2(E_L + E_scatt)
end






function get_intensity_over_an_angle(
    problem::LinearOptics{Scalar},
    atoms_states::AbstractVector,
    θ::Number,
)
    N = problem.atoms.N
    β = view(atoms_states, 1:N)

    ϕ_range = range(0, 2π, length = 30)
    vr = view(problem.atoms.r, :, :)

    complex_intensity = zeros(ComplexF64, N)
    total_intensity = 0.0

    for k = 1:length(ϕ_range)
        ϕ = ϕ_range[k]
        Threads.@threads for j = 1:N
            complex_intensity[j] =
                cis(
                    -k₀ * (
                        vr[1, j] * sin(θ) * cos(ϕ) +
                        vr[2, j] * sin(θ) * sin(ϕ) +
                        vr[3, j] * cos(θ)
                    ),
                ) * β[j]
        end
        total_intensity += abs2(sum(complex_intensity))
    end
    return total_intensity / length(ϕ_range)
end

function get_intensity_over_an_angle_shared(
    problem::LinearOptics{Scalar},
    β::AbstractVector,
    θ::Number,
    r_shared,
)
    N = problem.atoms.N

    ϕ_range = range(0, 2π, length = 30)

    complex_intensity = zeros(ComplexF64, N)
    total_intensity = 0.0

    for k = 1:length(ϕ_range)
        ϕ = ϕ_range[k]
        for j = 1:N
            complex_intensity[j] =
                cis(
                    -k₀ * (
                        r_shared[1, j] * sin(θ) * cos(ϕ) +
                        r_shared[2, j] * sin(θ) * sin(ϕ) +
                        r_shared[3, j] * cos(θ)
                    ),
                ) * β[j]
        end
        total_intensity += abs2(sum(complex_intensity))
    end
    return total_intensity / length(ϕ_range)
end

function get_intensity_over_an_angle(
    problem::LinearOptics{Scalar},
    atoms_states::Matrix,
    θ::Number,
)
    timeSteps = size(atoms_states, 2)
    intensities = zeros(timeSteps)

    # r_shared = SharedArray{Float64,2}(3, problem.atoms.N)
    r_shared = problem.atoms.r

    Threads.@threads for i = 1:timeSteps
        oneState = view(atoms_states, :, i)
        intensities[i] = get_intensity_over_an_angle_shared(problem, oneState, θ, r_shared)
    end

    return intensities
end


# ### --------------- MEAN FIELD ---------------
function _OnePoint_Intensity(physic::MeanField, laser, R⃗, sensor, β, scattering_func)
    
    Ω = laser_field(laser, sensor)
    
    r = norm(sensor)
    n̂ = sensor/r
    N = size(R⃗, 2)
    
    σ⁻ = β[1:N]
    σ⁺ = conj.(σ⁻)
    σᶻ = β[ (N+1) : end]

    term1 = abs2(-im*Ω)
    term2 = real(2Ω*(exp(-im*k₀*r)/(im*k₀*r))*ThreadsX.sum(σ⁺[j]*cis(+k₀ * (n̂⋅R⃗[:,j])) for j = 1:N), )
    term3 = _term3(σ⁻, σ⁺, n̂, R⃗)
    term4 = ThreadsX.sum((1 + σᶻ[j]) / 2 for j = 1:N)

    intensity_oneSensor = term1 + (Γ / 2)*term2 + (Γ/(2k₀*r))^2 * (term3 + term4)

    return real(intensity_oneSensor)
end
"""
    Hard core optimizations for term3. Check benchmarks folder for details.
"""
@views function _term3(σ⁻, σ⁺, n̂, R⃗)
    N = length(σ⁻)
    number_configurations = ((N^2)÷2 - N÷2)

    βₙₘ = Array{eltype(σ⁻)}(undef, number_configurations)
    cont = 1
    for n=1:N
        for m=(n+1):N            
            βₙₘ[cont] = σ⁺[n]*σ⁻[m]
            cont += 1
        end
    end
    

    rₙₘ = Array{eltype(R⃗)}(undef, 3, number_configurations)
    cont = 1
    for n=1:N
        r_n = R⃗[:,n]
        for m=(n+1):N
            rₙₘ[1,cont] = r_n[1] - R⃗[1,m]
            rₙₘ[2,cont] = r_n[2] - R⃗[2,m]
            rₙₘ[3,cont] = r_n[3] - R⃗[3,m]
            cont += 1
        end
    end

    intensity = ThreadsX.mapreduce(+, 1:number_configurations) do cont
        (  
           begin 
                @inbounds dot_n_r = n̂[1]*rₙₘ[1, cont] + n̂[2]*rₙₘ[2, cont] + n̂[3]*rₙₘ[3, cont]
                @inbounds βₙₘ[cont]*cis( -k₀*dot_n_r  )
           end 
        )
    end
    
    return 2real(intensity)
end


function get_intensity_over_an_angle(
    problem::NonLinearOptics{MeanField},
    atoms_states::Vector,
    θ::Float64,
)
    @debug "start : get intensity over an angle - NonLinearOptics{MeanField}"

    if is_integration_const_term_available(problem)
        Gₙₘ = problem.data[:Gₙₘ]
    else
        Gₙₘ = _get_meanfield_constant_term(problem.atoms, θ)
        problem.data[:Gₙₘ] = Gₙₘ
    end

    N = problem.atoms.N
    β = view(atoms_states, 1:N)
    z = view(atoms_states, (N+1):2N)

    βₙₘ = transpose(β * β') # I have to do "transpose" and NOT "adjoint = complex+tranpose"
    βₙₘ[diagind(βₙₘ)] .= (1 .+ z) ./ 2

    # IMPORTANT FOR NEXT LINE: you should use ELEMENT WISE multiplication.
    # Also, you can reduce memory allocation with inplace multiplication
    βₙₘ .*= Gₙₘ
    intensity = real(sum(βₙₘ))

    @debug "end  : get intensity over an angle - NonLinearOptics{MeanField}"
    return intensity
end

function is_integration_const_term_available(problem)
    if haskey(problem.data, :Gₙₘ)
        return true
    else
        return false
    end
end

function _get_meanfield_constant_term(atoms, Θ)
    N, r = atoms.N, atoms.r

    xₙₘ, yₙₘ, zₙₘ = get_xyz_distances(r)
    k₀sinΘ = k₀ * sin(Θ)
    cos_Θ = cos(Θ)

    Gₙₘ = Array{ComplexF64,2}(undef, N, N)
    @sync for n = 1:N
        Threads.@spawn for m = 1:N # we had to compute all terms, and not the upper part
            @inbounds Gₙₘ[n, m] =
                _constant_term_core_computation(xₙₘ, yₙₘ, zₙₘ, n, m, cos_Θ, k₀sinΘ)
        end
    end

    #=
        IF I want to compute only the upper part,
        I have to multiply all terms by π/2:   Gₙₘ = (π/2)*Gₙₘ

        I decided to don't make this, to don't appear with factors
        not mentioned on theory.
    =#
    #= 
        Before returning, we HAVE to do some memory cleaning,      
        EVEN losing some performance. 

        Without this cleaning, the garbage collector gets lost 
        outside this function - when many simulation occurs at the same time.
    =#
    xₙₘ = yₙₘ = zₙₘ = 1
    GC.gc()  # DO NOT DELETE
    return Gₙₘ
end
function _constant_term_core_computation(xₙₘ, yₙₘ, zₙₘ, n, m, cos_Θ, k₀sinΘ)
    a = zero(ComplexF64)
    a = cis(k₀ * zₙₘ[n, m] * cos_Θ) * besselj(0, k₀sinΘ * sqrt(xₙₘ[n, m]^2 + yₙₘ[n, m]^2))
    return a
end
function get_xyz_distances(r) # @memoize  --> creating warning, let's ignore it right now
    dimensions = size(r, 1)
    N = size(r, 2)

    xₙₘ = Array{Float64,2}(undef, N, N)
    yₙₘ = Array{Float64,2}(undef, N, N)
    zₙₘ = Array{Float64,2}(undef, N, N)

    r_shared = Array{Float64,2}(undef, dimensions, N)
    r_shared .= r
    @sync for n = 1:N
        r_n = view(r_shared, :, n)
        Threads.@spawn for m = 1:N
            xₙₘ[n, m] = r_n[1] - r_shared[1, m]
            yₙₘ[n, m] = r_n[2] - r_shared[2, m]
            zₙₘ[n, m] = r_n[3] - r_shared[3, m]
        end
    end

    return xₙₘ, yₙₘ, zₙₘ
end
