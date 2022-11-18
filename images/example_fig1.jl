using CoupledDipoles, Revise
using Random;
Random.seed!(1135);
using Plots, PlotThemes, LaTeXStrings
theme(:vibrant)

#= Spherical Cloud and Plane Wave laser at z-direction =#
N, ρ = 2_000, 1.0
atoms = Atom(Sphere(), sphere_inputs(N, ρ)...)

s, Δ = 1e-5, 0.0
laser = Laser(PlaneWave3D(), s, Δ)
simulation = LinearOptics(Scalar(), atoms, laser)

ωₙ, Γₙ = get_spectrum(simulation)
modes = classify_modes(simulation)

### ------ Choose beter looking modes--------------
## I choose these modes because they got clear interpretation 
loc_idx = modes.loc[end - 5]
sub_idx = modes.sub[4]
super_idx = modes.super[4]

## Small figures
fig_loc = plot_mode(simulation, loc_idx)
plot!(; title="(a) Loc", guidefont=25, title_position=:left)
fig_sub = plot_mode(simulation, sub_idx)
plot!(; title="(b) Sub", guidefont=25, title_position=:left)
fig_super = plot_mode(simulation, super_idx)
plot!(; title="(c) SR", guidefont=25, title_position=:left)

## Main Figure
fig_big = scatter(ωₙ[modes.sub], Γₙ[modes.sub]; markershape=:square, label="Sub")
scatter!(ωₙ[modes.loc], Γₙ[modes.loc]; markershape=:circle, label="Loc")
scatter!(ωₙ[modes.super], Γₙ[modes.super]; markershape=:utriangle, label="SR")

## Add black markers to indicate the modes shown
scatter!(
    [ωₙ[loc_idx]], [Γₙ[loc_idx]]; label="", color=:black, markershape=:square, markersize=7
)
scatter!(
    [ωₙ[sub_idx]], [Γₙ[sub_idx]]; label="", color=:black, markershape=:circle, markersize=7
)
scatter!(
    [ωₙ[super_idx]],
    [Γₙ[super_idx]];
    label="",
    color=:black,
    markershape=:utriangle,
    markersize=7,
)
plot!(;
    left_margin=5Plots.mm,
    right_margin=5Plots.mm,
    top_margin=5Plots.mm,
    bottom_margin=5Plots.mm,
    yscale=:log10,
    framestyle=:box,
    tickfont=20,
    guidefont=25,
    legendfont=17,
    xlabel=L"\omega_n/\Gamma",
    ylabel=L"\Gamma_n/\Gamma",
    legend=:bottomleft,
)

### --------Edges of Lorentz-Lorenz shift ---------
## Original formula of Lorentz-Lorenz shift cames from Eq 4 in :
## Microscopic and Macroscopic Signatures of 3D Anderson Localization of Light (2019)
## (using α = 0.5)
## α = 0.5
## λ_temp = 2π/k₀
## δ_cP = (λ_temp^3*ρ/(8*π^2)) + 0.5*sqrt(3*α*ρ*λ_temp^3/(4*π^2) -1)
## δ_cM = (λ_temp^3*ρ/(8*π^2)) - 0.5*sqrt(3*α*ρ*λ_temp^3/(4*π^2) -1)

k₀ = CoupledDipoles.k₀
## EPL Paper's formula
δ_cP = (π * ρ / (k₀^3)) + 0.5 * sqrt(3π * ρ / (k₀^3) - 1)
δ_cM = (π * ρ / (k₀^3)) - 0.5 * sqrt(3π * ρ / (k₀^3) - 1)

vline!([-δ_cP, -δ_cP]; label="", linestyle=:dashdot, color=:red, linewidth=5)
vline!([-δ_cM, -δ_cM]; label="", linestyle=:dashdot, color=:red, linewidth=5)

### ------------------ Figure 1--------------------
l = @layout([a{0.75w} [b; c; d]])
plot(fig_big, fig_super, fig_sub, fig_loc; layout=l, size=(1500, 900))
plot!(;
    left_margin=10Plots.mm,
    right_margin=5Plots.mm,
    top_margin=5Plots.mm,
    bottom_margin=5Plots.mm,
)

for fmt in ["png", "svg"]
    savefig("Fig1.$(fmt)")
end