
struct initial_conditions{T<:Vector{Float64}}
	
	a::T	# Rate-state-friction parameter 'a'
	b::T	# Rate-state-friction parameter 'b'
	σo::T	# Initial normal stress along depth
	τo::T	# Initial shear stress along depth 
	Lc::T	# Rate-state-friction parameter 'Lc'
end

function set_initial_conditions(fault_y, out_path)
	"""
		Set the initial conditions here 
	"""

	a = 0.015*ones(size(fault_y))
	b = 0.019*ones(size(fault_y))
	σo = 50.0e6*ones(size(fault_y))
	τo = 22.0e6*ones(size(fault_y))
	Lc = 8.0e-3*ones(size(fault_y))

	# Depth points for interpolation
	fault_depth = maximum(fault_y) - minimum(fault_y)
	p1 = 2.0e3
	p2 = 12.0e3
	p3 = 17.0e3
	p4 = 20.0e3
	
	#if p4 > fault_depth
	#	@error("Check the initial conditions wrt fault depth")
	#end

	# Transform fault_y to earth coordinates: Top surface = 0 km
	fy_trans = abs.(fault_y .- maximum(fault_y))

	# Get the fault index for each of these depths
	depth_id1 = findall(fy_trans .<= p1)
	depth_id2 = findall(p1 .< fy_trans .<= p2)
	depth_id3 = findall(p2 .< fy_trans .<= p3)
	depth_id4 = findall(p3 .< fy_trans .<= p4)
	depth_id5 = findall(fy_trans .> p4)

	# Set initial shear stress
	τo[depth_id1] .= 0.01e6 .+ 1e6*(30.0 - 0.01)/(p1 - 0.0) .* (fy_trans[depth_id1] .- 0.0)
	τo[depth_id2] .= 30.0e6
	τo[depth_id3] .= 30.0e6 .+ 1e6*(22.5 - 30.0)/(p3 - p2) .* (fy_trans[depth_id3] .- p2)
	τo[depth_id4] .= 22.5e6

	# Set initial normal stress 
	σo[depth_id1] .= 10e6 .+ 1e6*(50.0 - 10.0)/(p1 - 0.0) .* (fy_trans[depth_id1] .- 0.0)
	σo[depth_id2] .= 50.0e6
	σo[depth_id3] .= 50.0e6
	σo[depth_id4] .= 50.0e6


	# Set rate-state-friction parameters 
	# b is constant, set a values given a-b values

	# Uncomment for shallow velocity strengthening
	a[depth_id1] .= b[depth_id1] .- 0.003 .+ (-0.0041 + 0.003)/(p1 - 0.0) .* (fy_trans[depth_id1] .- 0.0)

	#a[depth_id1] .= b[depth_id1] .- 0.0041
	a[depth_id2] .= b[depth_id2] .- 0.0041
	a[depth_id3] .= b[depth_id3] .- 0.0041 .+ (0.015 + 0.0041)/(p3 - p2) .* (fy_trans[depth_id3] .- p2)
	a[depth_id4] .= b[depth_id4] .+ 0.015 .+ (0.0024 - 0.0015)/(p4 - p3) .* (fy_trans[depth_id4] .- p3)
	a[depth_id5] .= b[depth_id5] .+ 0.0047


	init_conds = hcat(fy_trans/1e3,  σo/1e6, τo/1e6, a, b, Lc)
	# Save initial conditions to file 
	open(string(out_path,"/initial_conditions.out"), "w") do io
		write(io, "fault_depth(km) normal_stress(MPa) shear_stress(MPa) friction_a friction_b friction_Lc\n")
		writedlm(io, init_conds)   
	end

	return initial_conditions(a, b, σo, τo, Lc)
end
