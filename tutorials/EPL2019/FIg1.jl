using CoupledDipoles, Revise
using Random;
Random.seed!(1134);
using GLMakie;
Makie.inline!(true);
using PrettyNumbers

function get_my_yticks(numero)
    io = IOBuffer()
    pretty_number(io, numero, show_significand = false)
    numero_bonito = io |> take! |> String
    return numero_bonito
end

#=
    Simulation is Simple
=#
N, ρ = 2500, 0.3
atoms = Atom(CoupledDipoles.Cube(), cube_inputs(N, ρ)...)

k₀ = 1;
λ = 2π / k₀;
w₀, s, Δ = 3λ, 1e-6, 1.0
laser = Laser(Gaussian3D(w₀), s, Δ)

simulation = LinearOptics(Scalar(), atoms, laser)
ωₙ, Γₙ = get_spectrum(simulation)
modes = classify_modes(simulation)



#=
    A lot of details for the Plot
=#
f = Figure(
    backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1000, 700),
    fontsize = 25,
)
ga = f[1, 1] = GridLayout()
gb = f[1, 2] = GridLayout()
colsize!(f.layout, 1, Relative(7 / 10)) # main figure is larger than subplots

y_to_show = [1e-15, 1e-10, 1e-5, 10]
myyticks = (y_to_show, get_my_yticks.(y_to_show))

axmain = Axis(
    ga[2:4, 1],
    yticks = myyticks,
    xlabel = L"\omega_n/\Gamma",
    ylabel = L"\Gamma_n/\Gamma",
    yscale = log10,
)
axsuper = Axis(
    gb[1, 1],
    yticks = myyticks,
    xlabel = L"|r - r_{cm}|",
    ylabel = L"|\psi|^2",
    yscale = log10,
)
axsub = Axis(
    gb[2, 1],
    yticks = myyticks,
    xlabel = L"|r - r_{cm}|",
    ylabel = L"|\psi|^2",
    yscale = log10,
)
axloc = Axis(
    gb[3, 1],
    yticks = myyticks,
    xlabel = L"|r - r_{cm}|",
    ylabel = L"|\psi|^2",
    yscale = log10,
)

# remove more spacing, and set the yticks
for ax in [axsuper, axsub, axloc]
    ax.xlabelpadding = -15
    ax.ylabelpadding = 4
    ax.yticks = myyticks
end

# overlay each mode type
scatter!(axmain, ωₙ[modes.super], Γₙ[modes.super], color = (:blue, 0.5))
scatter!(axmain, ωₙ[modes.sub], Γₙ[modes.sub], color = (:red, 0.5))
scatter!(axmain, ωₙ[modes.loc], Γₙ[modes.loc], color = (:green, 0.5))
# "-1.5" is not realistic, is just a example to test the vline function
vlines!(axmain, [-1.5, 1.5], color = :black, linestyle = :dash)

# hightlight 3 modes to see spatial profile
scatter!(
    axmain,
    [ωₙ[modes.super[1]]],
    [Γₙ[modes.super[1]]],
    color = :blue,
    markersize = 30,
    strokewidth = 5,
    strokecolor = :white,
)
scatter!(
    axmain,
    [ωₙ[modes.sub[1]]],
    [Γₙ[modes.sub[1]]],
    color = :red,
    markersize = 30,
    strokewidth = 5,
    strokecolor = :white,
)
scatter!(
    axmain,
    [ωₙ[N]],
    [Γₙ[N]],
    color = :green,
    markersize = 30,
    strokewidth = 5,
    strokecolor = :white,
)

# show spatial profile in the subplot
DCM, ψ² = get_spatial_profile_single_mode(simulation, modes.super[1])
l1 = scatter!(axsuper, DCM, ψ², color = :blue)

DCM, ψ² = get_spatial_profile_single_mode(simulation, modes.sub[1])
l2 = scatter!(axsub, DCM, ψ², color = :red)

DCM, ψ² = get_spatial_profile_single_mode(simulation, N)
l3 = scatter!(axloc, DCM, ψ², color = :green)

# add the legend 
labels = ["super", "sub", "loc"]
leg = Legend(ga[1, :], [l1, l2, l3], labels)
leg.tellheight = true
leg.orientation = :horizontal
rowsize!(ga, 1, Relative(1 / 50))

# remove spacing between figures
colgap!(f.layout, 3)
rowgap!(f.layout, 0)
rowgap!(ga, 25)
rowgap!(gb, 1)

# finish with subplot labels
for (label, layout) in zip(["(a)", "(b)", "(c)"], [gb[1, 1], gb[2, 1], gb[3, 1]])
    Label(
        layout[1, 1, TopLeft()],
        label,
        textsize = 20,
        padding = (0, 0, 0, 0),
        halign = :right,
    )
end
f
