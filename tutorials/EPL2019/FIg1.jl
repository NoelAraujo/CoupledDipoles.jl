using CoupledDipoles, Revise
using GLMakie;
Makie.inline!(true);


include("params.jl")
atoms = Atom(CoupledDipoles.Cube(), cube_inputs(N, ρ)...)
laser = Laser(Gaussian3D(w₀), s, Δ)
simulation = LinearOptics(Scalar(), atoms, laser)

@time ωₙ, Γₙ = get_spectrum(simulation)
modes = classify_modes(simulation)
ξₙ, R¹ₙ = get_localization_length(simulation)



#=
    Plot
=#
include("my_theme.jl")
function whatToSee(modetype; kwargs...)
    R¹ₙ = kwargs[:R¹ₙ]
    _, idxMax = findmax(R¹ₙ[modetype])
    return modetype[idxMax]
end


fig1_EPL = with_theme(Fig1_theme) do
    f = Figure(
        backgroundcolor = :white,#RGBf(0.98, 0.98, 0.98)
        resolution = (1000, 700),
        fontsize = 25,
    )
    ga = f[1, 1] = GridLayout()
    gb = f[1, 2] = GridLayout()
    colsize!(f.layout, 1, Relative(7 / 10)) # main figure is larger than subplots

    y_to_show = [1e-20, 1e-15, 1e-10, 1e-5, 10]
    myyticks = (y_to_show, get_my_yticks.(y_to_show))

    axmain = Axis(
        ga[2:4, 1],
        yscale = log10,
        yticks = myyticks,
        xlabel = L"\omega_n/\Gamma",
        ylabel = L"\Gamma_n/\Gamma",
    )
    axsubplots = [
        Axis(
            gb[nrow, 1],
            yscale = log10,
            yticks = myyticks,
            xlabel = L"|r - r_{cm}|",
            ylabel = L"|\psi|^2",
        ) for nrow = 1:3
    ]

    ll = []
    for (gbrow, modetype, color) ∈
        zip(axsubplots, [modes.super, modes.sub, modes.loc], [myColor1, myColor2, myColor3])
        oneSc = scatter!(axmain, ωₙ[modetype], Γₙ[modetype])
        push!(ll, oneSc)
        idxShow = whatToSee(modetype; R¹ₙ = R¹ₙ)

        DCM, ψ² = get_spatial_profile_single_mode(simulation, idxShow)
        scatter!(gbrow, DCM, abs.(ψ²), color = color)
    end
    for (label, layout) in zip(["(a)", "(b)", "(c)"], [gb[1, 1], gb[2, 1], gb[3, 1]])
        Label(
            layout[1, 1, TopRight()],
            label,
            textsize = 20,
            padding = (0, 0, 0, 0),
            halign = :right,
        )
    end
    δ_cP = (π * ρ / (k₀^3)) + 0.5 * sqrt(3π * ρ / (k₀^3) - 1)
    δ_cM = (π * ρ / (k₀^3)) - 0.5 * sqrt(3π * ρ / (k₀^3) - 1)
    vlines!(axmain, [-δ_cP, -δ_cM], color = :black, linestyle = :dash)

    for modetype ∈ [modes.super, modes.sub, modes.loc]
        idxShow = whatToSee(modetype, R¹ₙ = R¹ₙ)

        scatter!(
            axmain,
            [ωₙ[idxShow]],
            [Γₙ[idxShow]],
            markersize = 30,
            strokewidth = 5,
            strokecolor = :white,
        )
    end

    labels = ["super", "sub", "loc"]
    leg = Legend(ga[1, :], ll, labels)
    leg.tellheight = true
    leg.orientation = :horizontal
    rowsize!(ga, 1, Relative(1 / 50))
    f
end