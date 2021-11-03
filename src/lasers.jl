apply_laser_over_atoms(laser, atoms)     = _apply_laser(laser, atoms)
apply_laser_over_sensors(laser, sensors) = _apply_laser(laser, sensors)

apply_laser_over_oneAtom(laser, oneAtom) = _apply_laser_over_oneCoordinate(laser, oneAtom)
apply_laser_over_oneSensor(laser, oneSensor) = _apply_laser_over_oneCoordinate(laser, oneSensor)

function _apply_laser(laser, spatial_coordinates)
    N = _get_number_elements(spatial_coordinates)
    _core_LaserFunction = _get_core_LaserFunction(laser)
    _spatial_matrix = _get_spatial_matrix(spatial_coordinates)
    E₀ = estimate_E₀(laser)

    E = zeros(ComplexF64, N)
    matrix_view = view(_spatial_matrix, :, :)

    Threads.@threads for n in 1:N
        oneCoordinate = matrix_view[:,n]
        @inbounds E[n] = _core_LaserFunction(oneCoordinate, E₀, laser)
    end
    
    return E
end
function _apply_laser_over_oneCoordinate(laser, oneCoordinate)
    _core_LaserFunction = _get_core_LaserFunction(laser)
    E₀ = estimate_E₀(laser)    
    E = _core_LaserFunction(oneCoordinate, E₀, laser)

    return E
end


_get_number_elements(r::AbstractMatrix) = size(r, 2)
_get_number_elements(atoms::Atom{T}) where T = atoms.N

_get_spatial_matrix(r::AbstractMatrix) = r
_get_spatial_matrix(atoms::Atom{T}) where T = atoms.r

_get_core_LaserFunction(laser::Laser{PlaneWave3D}) = _core_PlaneWave3D
_get_core_LaserFunction(laser::Laser{Gaussian3D})  = _core_Gaussian3D


function _core_PlaneWave3D(oneAtom, E₀, laser)
    direction = view(laser.pump.direction, :)
    Eᵢ = E₀*cis(+k₀*(  
                      oneAtom[1]*direction[1] 
                    + oneAtom[2]*direction[2] 
                    + oneAtom[3]*direction[3]
                    )
              )
    return Eᵢ
end

function _core_Gaussian3D(oneAtom, E₀, laser)
    w₀ = laser.pump.w₀
    x, y, z = oneAtom[1], oneAtom[2], oneAtom[3]

    ## This formula is stable for z==0
    # Ref: Eq 3.11 from "CHAPTER 3. PROPAGATION AND FOCUSING OF OPTICAL FIELDS"
    denominator_factor = 1 + 2im * z / (k₀ * w₀^2)
    Eᵢ = E₀ * cis(+k₀ * z)
    Eᵢ = Eᵢ * exp(-(x^2 + y^2) / (denominator_factor * w₀^2))
    Eᵢ = Eᵢ / denominator_factor
    
    return Eᵢ
end