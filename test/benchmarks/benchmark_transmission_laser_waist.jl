using CoupledDipoles, Revise
using ProgressMeter
using Plots
using Random

function find_min_transmission(N, ρ, laser_increase_range, Δ_range, scatt_func)
    Random.seed!(2243)
    new_cloud = Atom(CoupledDipoles.Cylinder(), cylinder_inputs(N, ρ)...)
    s = 1e-5

    Ts = zeros(length(laser_increase_range), length(Δ_range))
    p = Progress(length(laser_increase_range) * length(Δ_range); showspeed=true)

    #=
        Multithreading did not returned good enough perfomance here.
        If you try, make sure to use spawn, because the process do not
        takes the exact amount of time.
    =#
    for l_idx in 1:length(laser_increase_range)
        laser_increase = laser_increase_range[l_idx]
        w₀ = size(new_cloud) * laser_increase

        for idx in 1:length(Δ_range)
            Δ = Δ_range[idx]

            _laser = Laser(Gaussian3D(w₀), s, Δ)
            _problem = LinearOptics(Scalar(), new_cloud, _laser)
            _βₙ = steady_state(_problem)

            Ts[l_idx, idx] = transmission(_problem, _βₙ, scatt_func)
            ProgressMeter.next!(p)
        end
    end

    return Ts
end

N = 500
ρ = 0.4
laser_increase_range = [2.0, 2.25, 2.5, 3.0]
Δ_range = range(-100, 100; length=20)

scatt_func = scattering_fuction(:nearField, :ThreeD)

Ts = find_min_transmission(N, ρ, laser_increase_range, Δ_range, scatt_func);

chosen = :Spectral_10
color_pallete = cgrad(chosen; rev=false).colors
n_of_colors = length(color_pallete)
my_colors = cgrad(color_pallete[round.(Int, range(1; stop=n_of_colors, length=length(laser_increase_range)))])

plotly()
plot(; xlabel="Δ", ylabel="Transmission")
for (l_idx, laser_increase) in enumerate(laser_increase_range)
    plot!(Δ_range, Ts[l_idx, :]; lw=4, label="w₀ = R*$(laser_increase)", c=my_colors[l_idx])
end
plot!(; title="N=$(N), ρ=$(ρ)", ylims=(0, 1.2), legend=:bottomright)
hline!([1.0]; label="", c=:black)
