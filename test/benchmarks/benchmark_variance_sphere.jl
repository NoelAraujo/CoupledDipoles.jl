using CoupledDipoles
using Random

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

const k₀ = 1;
const λ₀ = 2π / k₀;

# atoms setup
atomType = CoupledDipoles.Cylinder()
N = 7500
ρλ³ = 45
ρ_over_k₀³ = ρλ³ / (2π)^3

radiusFactor = 5
kR = radiusFactor*λ₀
kh = N/(ρ_over_k₀³*π*kR^2)
atomInfo = N, kR, kh

# laser setup
waistFactor = 2
w₀, s, Δ = waistFactor*λ₀, 1e-6, 1.0

# creating problem
Random.seed!(100)
cloud = Atom(atomType, atomInfo...)
laser = Laser(Gaussian3D(w₀), s, Δ) 
problem = LinearOptics(Scalar(), cloud, laser)




# creating sensors
farFieldRadius = how_far_is_FarField(cloud)
x = farFieldRadius .* [sin(θ) * sin(ϕ) for θ in θ, ϕ in ϕ]
y = farFieldRadius .* [sin(θ) * cos(ϕ) for θ in θ, ϕ in ϕ]
z = farFieldRadius .* [cos(θ) for θ in θ, ϕ in ϕ]
sensors = Array([reshape(x, 200 * 200) reshape(y, 200 * 200) reshape(z, 200 * 200)]')




# getting light intensities
βₙ = steady_state(problem)
scattering_func = scattering_fuction(:farField, :ThreeD)
intensities = scattering_intensity(problem, βₙ,sensors, scattering_func)
intensities[findall(intensities .< 1e-20)] .= 0.0




# Plotting
fig, ax, pltobj = mesh(
    getMesh(x, y, z),
    color = log10.(intensities),
    shading = true,
    ambient = Vec3f0(0.85, 0.85, 0.85),
    figure = (resolution = (680, 830), fontsize = 15, backgroundcolor = :white),
    axis = (
        title = "w₀=$(waistFactor)λ, kR=$(radiusFactor)λ",
        type = Axis3,
        aspect = :data,
    ),
)
cbar =
    Colorbar(fig, pltobj, label = "log10(Intensity)", flipaxis = false)
fig[1, 2] = cbar
fig
