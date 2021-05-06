"""
    get_spatial_profile_data_to_save(problem, mode_index)

return (data = ( x=DCM, y=ψ² ), fit=(x=x_fit, y=y_fit) )
"""
function get_spatial_profile_data_to_save(problem, mode_index)
    prepare_spectrum(problem)
    DCM, ψ² = get_spatial_profile_single_mode(problem, mode_index)
    x_fit, y_fit = get_decay_fit(problem, mode_index)
    return (data=(x=DCM, y=ψ²), fit=(x=x_fit, y=y_fit))
end

function plot_atoms_and_mode(problem, mode_index)
    # pyplot()
    # compute eigenvalues if not available
    prepare_spectrum(problem)
    DCM, ψ² = get_spatial_profile_single_mode(problem, mode_index)
    x_fit, y_fit = get_decay_fit(problem, mode_index)

    fig_space = plot_space(problem.atoms, ψ²)
    plot!(; guidefont=16, tickfont=15, size=(800, 600))

    fig_mode = scatter(DCM, ψ²; label="Mode Index = $(mode_index)", yscale=:log10)
    plot!(
        x_fit,
        y_fit;
        xlabel=L"|\mathbf{r}-\mathbf{r}_{cm}|",
        ylabel=L"|\psi_n|^2",
        size=(800, 600),
        guidefont=16,
        tickfont=15,
        legendfont=14,
        label="",
        lw=5,
        linestyle=:dash,
        c=:black,
        grid=false,
        framestyle=:box,
    )

    return fig_space, fig_mode
end

function plot_space(atoms::T where {T<:ThreeD}, ψ²)
    return scatter3d(
        get_X_axes(atoms),
        get_Y_axes(atoms),
        get_Z_axes(atoms);
        label="",
        xlabel="x",
        ylabel="x",
        zlabel="z",
        zcolor=ψ²,
        colorbar_title="|ψ|²",
        camera=(45, 30),
    )
end

prepare_spectrum(problem) = get_spectrum(problem)

function select_matrix_axes(static_array, component)
    N = length(static_array)
    dimensions = length(static_array[1])
    if component ≤ dimensions
        return [static_array[n][component] for n in 1:N]
    else
        return [zero(eltype(static_array[n])) for n in 1:N]
    end
end

get_X_axes(static_array::Vector{StaticArrays.SArray}) = select_matrix_axes(static_array, 1)
get_Y_axes(static_array::Vector{StaticArrays.SArray}) = select_matrix_axes(static_array, 2)
get_Z_axes(static_array::Vector{StaticArrays.SArray}) = select_matrix_axes(static_array, 3)

get_X_axes(atoms::T where {T<:ThreeD}) = select_matrix_axes(atoms.r, 1)
get_Y_axes(atoms::T where {T<:ThreeD}) = select_matrix_axes(atoms.r, 2)
get_Z_axes(atoms::T where {T<:ThreeD}) = select_matrix_axes(atoms.r, 3)

get_X_axes(atoms::T where {T<:TwoD}) = select_matrix_axes(atoms.r, 1)
get_Y_axes(atoms::T where {T<:TwoD}) = select_matrix_axes(atoms.r, 2)
get_Z_axes(atoms::T where {T<:TwoD}) = zeros(length(atoms.r))

get_atoms_matrix(atoms) = convert_StaticArray_to_matrix(atoms.r)

function get_decay_fit(problem, mode_index)
    DCM, ψ² = get_spatial_profile_single_mode(problem, mode_index)

    x_use, y_use = select_points(DCM, ψ²)
    A_L₁, ξ_L₁, R1, y_fit_L₁ = get_coeffs_L1(x_use, y_use)

    if is_R1_localized(R1)
        return x_use, y_fit_L₁
    else
        A_L₂, ξ_L₂, R2, y_fit_L₂ = get_coeffs_L2(x_use, y_use)
        return x_use, y_fit_L₂
    end
end
