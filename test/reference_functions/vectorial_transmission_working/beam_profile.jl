function laser_beam_profile(x, y, z, E0, w0, z0, k)
    # Calculate the radial distance from the beam axis
    r = sqrt.(x.^2 .+ y.^2)

    # Calculate the beam profile
    beam_profile = E0./ sqrt.(1 .+ (z./ z0).^2).*
                   exp.(-r.^2/ (w0^2* (1 .+ (z./ z0).^2))).*
                   exp.(im * k.* (r.^2)./ (2 * z.* (1 .+ (z0./ z).^2))).*
                   exp.(im * (k * z - atan(z./ z0)))

    return beam_profile
end