function NonLinearOptics(physic::MeanField, atoms, laser)
	@debug "start: NonLinearOptics - $( typeof(physic) )"
	laser_copy = deepcopy(laser) # i don't want change the original laser
	if !all(laser_copy.polarization .== 0)
		laser_copy.polarization = [0, 0, 0] # this is an indication to use scalar functions
	end
	excitations = Dict()
	data = Dict()

	@debug "end  : NonLinearOptics - $( typeof(physic) )"
	return NonLinearOptics(physic, atoms, laser_copy, excitations, data)
end


function NonLinearOptics(physic::PairCorrelation, atoms, laser)
	@debug "start: NonLinearOptics - $( typeof(physic) )"
	laser_copy = deepcopy(laser) # i don't want change the original laser
	if !all(laser_copy.polarization .== 0)
		laser_copy.polarization = [0, 0, 0] # this is an indication to use scalar functions
	end
	excitations = Dict()
	data = Dict()

	@debug "end  : NonLinearOptics - $( typeof(physic) )"
	return NonLinearOptics(physic, atoms, laser_copy, excitations, data)
end
