"""
    get_spatial_profile_single_mode_and_fit(problem, mode_index)

return (data = ( x=DCM, y=ψ² ), fit=(x=x_fit, y=y_fit) )
"""
function get_spatial_profile_single_mode_and_fit(problem, mode_index)
    DCM, ψ² = get_spatial_profile_single_mode(problem, mode_index)
    x_fit, y_fit = get_decay_fit(problem, mode_index)
    return (DCM=DCM, ψ²=ψ², fit=(x=x_fit, y=y_fit))
end

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
