using CoupledDipoles
using Plots, ProgressMeter
using Revise

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

# cloud settings
N = 75
ρ = 0.5

cloud = Atom(CoupledDipoles.Cylinder(), cylinder_inputs(N, ρ; h=5π)...)
R = cloud.sizes.R

# laser settings
w₀ = 1.5 * 2π
s = 1e-5

# create transmission depending on detunning
Δ_range = range(-50, 50; length=30)
T = zeros(length(Δ_range))
let
p = Progress(length(Δ_range); showspeed=true)
Threads.@threads for idx in 1:length(Δ_range)
    Δ = Δ_range[idx]

    _laser = Laser(Gaussian3D(w₀), s, Δ)
    # _problem = LinearOptics(Scalar(), cloud, _laser)
    _problem = LinearOptics(Vectorial(), cloud, _laser)
    _βₙ = steady_state(_problem)

    # _problem = NonLinearOptics(MeanField(), cloud, _laser)
    # _βₙ = steady_state(_problem)

    T[idx] = transmission(_problem, _βₙ)[1]
    ProgressMeter.next!(p)
end

plot(
    Δ_range,
    T;
    label="",
    ylims = (0, 2),
    size=(800, 400),
    lw=3,
    legend=:bottomright,
)
hline!([1]; linestyle=:dash, c=:black, label="")
xlabel!("Δ")
ylabel!("Transmission")
display(title!("Cylinder : N=$(N), ρ=$(round(ρ,digits=3)), R=$(round(R,digits=2)), w₀=$(round(w₀,digits=2)), λ=$(round(2π,digits=2))"))
end
#=
    The code below is to visualize the intensity over the space in the Far Field limit.

    Mesh was copied from: https://github.com/lazarusA/MakieNotebooks/blob/main/scripts/sphericalHarmonicsPlot.jl
=#

# using WGLMakie
# using WGLMakie: get_dim, surface_normals
using CairoMakie
CairoMakie.activate!()

using GeometryBasics
using WGLMakie: get_dim, surface_normals
function getMesh(x, y, z)
    positions = vec(
        map(CartesianIndices(z)) do i
            GeometryBasics.Point{3,Float32}(get_dim(x, i, 1, size(z)), get_dim(y, i, 2, size(z)), z[i])
        end,
    )
    faces = decompose(GLTriangleFace, Rect2D(0.0f0, 0.0f0, 1.0f0, 1.0f0), size(z))
    normals = surface_normals(x, y, z)
    vertices = GeometryBasics.meta(positions; normals=normals)
    meshObj = GeometryBasics.Mesh(vertices, faces)
    return meshObj
end

θ = LinRange(0, π, 20)
ϕ = LinRange(0, 2π, 20)

farFieldRadius = CoupledDipoles.how_far_is_FarField(cloud)
x = farFieldRadius .* [sin(θ) * sin(ϕ) for θ in θ, ϕ in ϕ]
y = farFieldRadius .* [sin(θ) * cos(ϕ) for θ in θ, ϕ in ϕ]
z = farFieldRadius .* [cos(θ) for θ in θ, ϕ in ϕ]
sensors = Array([reshape(x, 20 * 20) reshape(y, 20 * 20) reshape(z, 20 * 20)]')

_Δ = 0.0
# _laser = Laser(Gaussian3D(w₀), s, _Δ; polarization=[-1,im,0]/√2)
_laser = Laser(PlaneWave3D(), s, _Δ)
# _problem = LinearOptics(Scalar(), cloud, _laser)
_problem = LinearOptics(Vectorial(), cloud, _laser)
# _βₙ = Matrix(transpose([steady_state(_problem)]))
_βₙ =steady_state(_problem)

intensities = scattered_intensity(_problem, _βₙ, sensors; regime=:far_field)
intensities[findall(intensities .< 1e-20)] .= 0.0

fig, ax, pltobj = mesh(
    getMesh(x, y, z);
    color=log10.(intensities),
    colormap=:inferno,
    shading=false,
    ambient=Vec3f0(0.85, 0.85, 0.85),
    # colorrange=(-10, -4),
    # figure=(resolution=(680, 830), fontsize=14, backgroundcolor=:white),
    # axis=(title="N=$(N), ρ=$(ρ),  w₀=R", type=Axis3, aspect=:data),
    axis=(title="Vectorial", type=Axis3, aspect=:data),
)
cbar = Colorbar(fig, pltobj; label="log10(Intensity)", flipaxis=false)
fig[1, 2] = cbar
fig
