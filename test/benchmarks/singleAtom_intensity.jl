using CoupledDipoles
using Random
using CairoMakie
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

n_points = 75
θ = LinRange(0, π, n_points)
ϕ = LinRange(0, 2π, n_points)

farFieldRadius = 50
x = farFieldRadius .* [sin(θ) * sin(ϕ) for θ in θ, ϕ in ϕ]
y = farFieldRadius .* [sin(θ) * cos(ϕ) for θ in θ, ϕ in ϕ]
z = farFieldRadius .* [cos(θ) for θ in θ, ϕ in ϕ]
sensors = Array([reshape(x, n_points * n_points) reshape(y, n_points * n_points) reshape(z, n_points * n_points)]')


# cloud settings
N = 1
kR = 10
Random.seed!(2044)
cloud = Atom(CoupledDipoles.Cylinder(), N, kR, kR)

# laser settings
w₀ = 4π
s = 1e-5
Δ = 0.0

laser = Laser(Gaussian3D(w₀), s, Δ; direction=[-1,0,0], polarization=[0,0,1])
# laser = Laser(PlaneWave3D(), s, Δ; direction=[0,0,1], polarization=[1,0,0]) #

sca_problem = LinearOptics(Scalar(), cloud, laser)
vec_problem = LinearOptics(Vectorial(), cloud, laser)

sca_βₙ = steady_state(sca_problem)
vec_βₙ = steady_state(vec_problem)

## MUST BE ON NEAR FIELD TO SEE ANY DIFFERENCE
sca_intensities = laser_and_scattered_intensity(sca_problem, sca_βₙ, sensors; regime=:near_field)
vec_intensities = laser_and_scattered_intensity(vec_problem, vec_βₙ, sensors; regime=:near_field)

sca_intensities = log10.(sca_intensities)
vec_intensities = log10.(vec_intensities)

CairoMakie.activate!()
let
    f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1100, 700))

    sca_ax3d = Axis3(f[1, 1], title = "Scalar", aspect=:data)
    vec_ax3d = Axis3(f[1, 2], title = "Vectorial", aspect=:data)

    sca_m = mesh!(
        sca_ax3d,
        getMesh(x, y, z);
        color=sca_intensities,
        colormap=:plasma,
        shading=false,
        ambient=Vec3f0(0.85, 0.85, 0.85),
        colorrange=extrema(sca_intensities),
    )
    vec_m = mesh!(
        vec_ax3d,
        getMesh(x, y, z);
        color=vec_intensities,
        colormap=:plasma,
        shading=false,
        ambient=Vec3f0(0.85, 0.85, 0.85),
        colorrange=extrema(sca_intensities),
    )
    sca_cbar = Colorbar(f, sca_m; label="log10( Intensity )", flipaxis=false,  vertical = false, width = Relative(4/5),ticks=WilkinsonTicks(3))
    # vec_cbar = Colorbar(f, vec_m; label="Intensity", flipaxis=false,  vertical = false, ticklabelsize = 19, width = Relative(4/5), ticks=WilkinsonTicks(3))
    f[2, :] = sca_cbar
    # f[2, 2] = vec_cbar
    f
end




using ColorSchemes
using WGLMakie
WGLMakie.activate!()


x = LinRange(-100, 100, 100)
y = LinRange(-100, 100, 100)
z = LinRange(-100, 100, 100)

sca_vol = [log10(scattered_intensity(sca_problem, sca_βₙ, Matrix([X Y Z]');regime=:near_field)[1]) for X ∈ x, Y ∈ y, Z ∈ z]
vec_vol = [log10(scattered_intensity(vec_problem, vec_βₙ, Matrix([X Y Z]');regime=:near_field)[1]) for X ∈ x, Y ∈ y, Z ∈ z]


myColorRange = (-11, -9)
myColorRange = extrema(sca_vol)

let
    fig = Figure(resolution = (1100, 700))

    sca_ax = Axis3(fig[1, 1], title = "Scalar", aspect=:data)
    vec_ax = Axis3(fig[1, 2], title = "Vectorial", aspect=:data)

    sca_plt = volumeslices!(sca_ax, x, y, z, sca_vol,
        colormap=cgrad( ColorSchemes.linear_kryw_0_100_c71_n256, rev=false),
        # colormap = cgrad( ColorSchemes.ocean, rev=false),
        colorrange=myColorRange
        )
    sca_plt[:update_yz][](100)
    sca_plt[:update_xz][](50)
    sca_plt[:update_xy][](1)

    vec_plt = volumeslices!(vec_ax, x, y, z, vec_vol,
        colormap=cgrad( ColorSchemes.linear_kryw_0_100_c71_n256, rev=false),
        # colormap=cgrad( ColorSchemes.ocean, rev=false),
        colorrange=myColorRange
        )
    vec_plt[:update_yz][](100)
    vec_plt[:update_xz][](50)
    vec_plt[:update_xy][](1)
    # zoom!(ax.scene, cameracontrols(ax.scene), 0.65)

    sca_cbar = Colorbar(fig, sca_plt; label="log10( Intensity )", flipaxis=false,  vertical = false, width = Relative(4/5),ticks=WilkinsonTicks(3))
    fig[2, :] = sca_cbar
    # fig[2, 2] = vec_cbar
    fig
end

let
    fig = Figure(resolution = (1100, 700))

    sca_ax = Axis3(fig[1, 1], title = "Scalar", aspect=:data)
    vec_ax = Axis3(fig[1, 2], title = "Vectorial", aspect=:data)

    sca_plt = volumeslices!(sca_ax, x, y, z, sca_vol,
        colormap=cgrad( ColorSchemes.linear_kryw_0_100_c71_n256, rev=false),
        colorrange=myColorRange
        )
    sca_plt[:update_yz][](100)
    sca_plt[:update_xz][](50)
    sca_plt[:update_xy][](1)

    vec_plt = volumeslices!(vec_ax, x, y, z, vec_vol,
        colormap=cgrad( ColorSchemes.linear_kryw_0_100_c71_n256, rev=false),
        colorrange=myColorRange
        )
    vec_plt[:update_yz][](100)
    vec_plt[:update_xz][](50)
    vec_plt[:update_xy][](1)
    # zoom!(ax.scene, cameracontrols(ax.scene), 0.65)

    sca_cbar = Colorbar(fig, sca_plt; label="log10( Intensity )", flipaxis=false,  vertical = false, width = Relative(4/5),ticks=WilkinsonTicks(3))
    fig[2, :] = sca_cbar
    # fig[2, 2] = vec_cbar
    fig
end





let
    fig = Figure(resolution = (1100, 700))
    vec_ax = Axis3(fig[1, 1], title = "Vectorial", aspect=:data)

    vec_plt = volumeslices!(vec_ax, x, y, z, vec_vol,
        colormap=cgrad( ColorSchemes.ocean, rev=false),
        # colorrange=(-12, -9)
        )
    vec_plt[:update_yz][](100)
    vec_plt[:update_xz][](50)
    vec_plt[:update_xy][](1)
    # zoom!(ax.scene, cameracontrols(ax.scene), 0.65)

    sca_cbar = Colorbar(fig, vec_plt; label="log10( Intensity )", flipaxis=false,  vertical = false, width = Relative(4/5),ticks=WilkinsonTicks(3))
    fig[2, :] = sca_cbar
    # fig[2, 2] = vec_cbar
    fig
end

# CairoMakie.activate!()
let
    hist(sca_vol[:], label="Scalar"; axis = (xlabel = "normalized log10(intensity)", ylabel = "",
        xlabelsize = 22, ylabelsize = 22))
    hist!(vec_vol[:], label="Vectorial")
    axislegend(position = :lt)
    current_figure()
end