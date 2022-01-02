using CoupledDipoles, Revise
using GLMakie;
Makie.inline!(true);

include("params.jl")
atoms = Atom(CoupledDipoles.Cube(), cube_inputs(N, ρ)...)
laser = Laser(Gaussian3D(w₀), s, Δ)
simulation = LinearOptics(Scalar(), atoms, laser)

modes = classify_modes(simulation)
ξₙ, R¹ₙ = get_localization_length(simulation)
IPRₙ = get_IPRs(simulation)



#=
    Plot
=#
include("my_theme.jl")
function whatToSee(modetype; kwargs... )
    R¹ₙ = kwargs[:R¹ₙ]
    _, idxMax = findmax(R¹ₙ[modetype])
    return modetype[idxMax]
end

fig2_EPL = with_theme(Fig2_theme) do
    f = Figure(
        backgroundcolor = :white,
        resolution = (800, 600),
        fontsize = 25,
    )
    ga = f[1, 1] = GridLayout()
    y_to_show = [1e-3, 1e-2, 1e-1, 1]
    myyticks = (y_to_show, get_my_yticks.(y_to_show))
    axmain = Axis(
        ga[2, 1],
        yticks = myyticks,
        xlabel = L"R^1",
        ylabel = "IPR",
        yscale = log10,
    )
    l3 = scatter!(axmain, R¹ₙ[modes.loc],   IPRₙ[modes.loc])
    l2 = scatter!(axmain, R¹ₙ[modes.sub],   IPRₙ[modes.sub] )
    l1 = scatter!(axmain, R¹ₙ[modes.super], IPRₙ[modes.super])

    vlines!(axmain, [0.65], color=:black, linestyle=:dash)
    
    labels = ["super", "sub", "loc"]
    leg = Legend(ga[1, 1], [l1, l2, l3], labels)
    leg.tellheight = true
    leg.orientation = :horizontal
    rowsize!(ga, 1, Relative(1 / 50))
    f 
end

inset_ax1 = add_box_inset(fig2_EPL; 
                left=542, right=780, 
                bottom=140, top=280)

for modetype ∈ [modes.loc, modes.sub, modes.super]    
    actual_index = whatToSee(modetype; R¹ₙ=R¹ₙ ) 
    DCM, ψ² = get_spatial_profile_single_mode(simulation, actual_index)
    scatter!(inset_ax1, DCM, log10.(ψ²))

    scatter!(
        fig2_EPL.layout.parent.content[1],
        [R¹ₙ[actual_index]],
        [IPRₙ[actual_index]],
        markersize = 30,
        strokewidth = 5,
        strokecolor = :white,
    )
end

fig2_EPL