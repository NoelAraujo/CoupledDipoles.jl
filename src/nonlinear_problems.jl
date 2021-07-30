function NonLinearOptics(physic::MeanField, atoms, laser)
    @debug "start: NonLinearOptics - $( typeof(physic) )"
    
    excitations = Dict()
    data = Dict()
    
    @debug "end  : NonLinearOptics - $( typeof(physic) )"
    return NonLinearOptics(physic, atoms, laser, excitations, data)
end