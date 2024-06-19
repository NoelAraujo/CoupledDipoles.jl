using CoupledDipoles
using CairoMakie, ProgressMeter

using LinearAlgebra, Random, CoupledDipoles, CairoMakie

Δ_range = range(-5,5; length=100)

R = 17.71830997686217301634315
L = 8.177681527782540982229875
N = 200
b0 = 4N/R^2 # for vectorial is different than scalar

# you need an initial value for Δ, just to create a problem type
s, Δ = 1e-3, 0.0 
w₀ =  R/3
polarization = [1, 0, 0]
direction = [0, 0, 1]

T = map(1:20) do j
    Random.seed!(j)
    println(j)
    # atoms = Atom(Cylinder(), N, R, L)
    atoms = Atom(CoupledDipoles.Sphere(gaussian=true), N, R)
    laser = Laser(Gaussian3D(w₀), s, Δ; polarization=[1,0,0])
    problem = LinearOptics(Vectorial(), atoms, laser)

    transmission(problem, Δ_range)
end;
mean_transmission = sum(T)/length(T)

let 
    # b0 = 9N/R^2 # uniform
    b0 = 3N/R^2 # gaussina

    f = Figure(size = (800, 600))
    ax = Axis(f[1, 1], xlabel = "Δ", ylabel = "Transmission")
    foreach(eachindex(T)) do j
        lines!(ax, Δ_range, T[j])
    end
    lines!(ax, Δ_range, exp.(-b0 ./ (1 .+ 4 .* Δ_range.^2)), linewidth=4, color=:black, label="Beer's Law")
    lines!(ax, Δ_range, mean_transmission, linewidth=4, color=:red, linestyle=:dash, label="Average")
    axislegend(ax, position = :rb)
    ylims!(ax, 0, 1.2)
    f
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



let 
    using LinearAlgebra
    atoms = Atom(Cube(), Matrix([1.0 1.0 1.0]'), 1.0)
    sensor = Matrix([-1000 -1000 -500]')
    β = [3, 4im, 5.0]
    E_n = CoupledDipoles._vectorial_scattering_near_field(atoms, β, sensor)
    E_f = CoupledDipoles._vectorial_scattering_far_field(atoms, β, sensor)
    E_n, E_f
end


@testset "Vectorial Scattering - Single Atom" begin
    using LinearAlgebra
    atoms = Atom(Cube(), Matrix([1.0 1.0 1.0]'), 1.0)
    sensor = Matrix([-1000 -1000 -500]')
    β = [3, 4im, 5.0]
    E_μ = CoupledDipoles._vectorial_scattering_far_field(atoms, β, sensor)

    R = norm(sensor)
    n = sensor./R
    nx, ny, nz = n[1], n[2], n[3]
    C = cis(-nx - ny - nz)
    E_x_expected = +(im/2)*(3/2)*(cis(R)/R)*C*((1 - nx^2)*3 + -nx*ny*4im - nx*nz*5)
    E_y_expected = +(im/2)*(3/2)*(cis(R)/R)*C*( -ny*nx*3 + (1 -ny^2)*4im - ny*nz*5)
    E_z_expected = +(im/2)*(3/2)*(cis(R)/R)*C*( -nz*nx*3 - nz*ny*4im + (1 - nz^2)*5)
    @test all(E_μ .≈ [E_x_expected, E_y_expected, E_z_expected])

    sensor = Matrix([100 -1000 0]')
    β = [1, -2im, 15.0]
    E_μ = CoupledDipoles._vectorial_scattering_far_field(atoms, β, sensor)

    R = norm(sensor)
    n = sensor./R
    nx, ny, nz = n[1], n[2], n[3]
    C = cis(-nx - ny)
    E_x_expected = +(im/2)*(3/2)*(cis(R)/R)*C*((1 - nx^2) + nx*ny*2im - nx*nz*15)
    E_y_expected = +(im/2)*(3/2)*(cis(R)/R)*C*( -ny*nx - (1 -ny^2)*2im - ny*nz*15)
    E_z_expected = +(im/2)*(3/2)*(cis(R)/R)*C*( -nz*nx + nz*ny*2im + (1 - nz^2)*15)

    @test all(E_μ .≈ [E_x_expected, E_y_expected, E_z_expected])
end