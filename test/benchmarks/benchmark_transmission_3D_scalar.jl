using CoupledDipoles
using Plots, ProgressMeter
using Revise, ThreadPools

# cloud settings
N = 1500
ρ = 0.3

cloud = Atom(CoupledDipoles.Sphere(), sphere_inputs(N, ρ)...)
R = size(cloud)

# laser settings
radial_increase = 1.0
w₀ = R * radial_increase
s = 1e-5


# create transmission depending on detunning
Δ_range = range(-50, 50, length = 50)
T = zeros(length(Δ_range))

p = Progress(length(Δ_range); showspeed = true)
ThreadPools.@qthreads for idx ∈ 1:length(Δ_range)
    Δ = Δ_range[idx]

    _laser = Laser(Gaussian3D(w₀), s, Δ)
    _problem = LinearOptics(Scalar(), cloud, _laser)
    _βₙ = get_steady_state(_problem)

    T[idx] = get_transmission(_problem, _βₙ)
    ProgressMeter.next!(p)
end


plot(
    Δ_range,
    T,
    label = "",
    ylims = (0, 1.2),
    size = (800, 400),
    lw = 3,
    legend = :bottomright,
)
hline!([1], linestyle = :dash, c = :black, label = "")
xlabel!("Δ")
ylabel!("Transmission")
title!("Sphere : N=$(N), ρ=$(ρ),  w₀=$(radial_increase)*R")

#=
    The code below is to visualize the intensity over the space in the Far Field limit.
    
    Mesh was copied from: https://github.com/lazarusA/MakieNotebooks/blob/main/scripts/sphericalHarmonicsPlot.jl
=#

using WGLMakie
using GeometryBasics
using WGLMakie: get_dim, surface_normals
function getMesh(x, y, z)
    positions = vec(
        map(CartesianIndices(z)) do i
            GeometryBasics.Point{3,Float32}(
                get_dim(x, i, 1, size(z)),
                get_dim(y, i, 2, size(z)),
                z[i],
            )
        end,
    )
    faces = decompose(GLTriangleFace, Rect2D(0.0f0, 0.0f0, 1.0f0, 1.0f0), size(z))
    normals = surface_normals(x, y, z)
    vertices = GeometryBasics.meta(positions; normals = normals)
    meshObj = GeometryBasics.Mesh(vertices, faces)
    meshObj
end

θ = LinRange(0, π, 200)
ϕ = LinRange(0, 2π, 200)

farFieldRadius = how_far_is_FarField(cloud)
x = farFieldRadius .* [sin(θ) * sin(ϕ) for θ in θ, ϕ in ϕ]
y = farFieldRadius .* [sin(θ) * cos(ϕ) for θ in θ, ϕ in ϕ]
z = farFieldRadius .* [cos(θ) for θ in θ, ϕ in ϕ]
sensors = Array([reshape(x, 200 * 200) reshape(y, 200 * 200) reshape(z, 200 * 200)]')

_Δ = 0.0
_laser = Laser(Gaussian3D(w₀), s, _Δ)
_problem = LinearOptics(Scalar(), cloud, _laser)
_βₙ = get_steady_state(_problem)

intensities = get_intensities_over_sensors(_problem, _βₙ, sensors)
intensities[findall(intensities .< 1e-20)] .= 0.0


fig, ax, pltobj = mesh(
    getMesh(x, y, z),
    color = log10.(intensities),
    colormap = :inferno,
    shading = true,
    ambient = Vec3f0(0.85, 0.85, 0.85),
    figure = (resolution = (680, 830), fontsize = 14, backgroundcolor = :white),
    axis = (
        title = "N=$(N), ρ=$(ρ),  w₀=$(radial_increase)*R",
        type = Axis3,
        aspect = :data,
    ),
)
cbar =
    Colorbar(fig, pltobj, olormap = :inferno, label = "log10(Intensity)", flipaxis = false)
fig[1, 2] = cbar
fig
