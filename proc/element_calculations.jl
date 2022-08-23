# Quasi-static element computations
function qs_element_calculator(meshbox, out_dir, glob, iglob, uFault, K)

	for el in meshbox.n_elements

		# Fault Boundary stuff
		if any(meshbox.boundary_names[:,el] .== :Right)
			fault_idx = iglob[:,end,el]
			glob.u[fault_idx]
		end

	end

end