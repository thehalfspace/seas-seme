
function get_output_locations(fault_y, fault_id)
	"""
		Set the output locations here.

		out_loc = depth from the surface in kms

		**NOTE**: Only works with vertical fault.
	"""

	top = maximum(fault_y)
	bottom = minimum(fault_y)

	# Depth from the top surface in km 
	out_loc = [0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0].*1e3

	# Transform out_loc to my fault coordinates 
	out_fault = abs.(out_loc .- top) 
	
	# global indices for out_locations 
	out_glob_id = zeros(size(out_fault))

	# on-fault indices for out_locations: useful for on-fault variables 
	out_fault_id = zeros(size(out_fault))

	for i in 1:length(out_fault)
		id = findmin(abs.(fault_y .- out_fault[i]))[2]

		out_fault_id[i] = Int(id)
		out_glob_id[i] = Int(fault_id[id])
	end
	out_fault_id, out_glob_id, out_loc./1e3
end