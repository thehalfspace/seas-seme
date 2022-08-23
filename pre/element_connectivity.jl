#######################################################
#######################################################
#
#	SOME INFO:
#   ----------
#
# edge = line 
# meshbox.n_boundaries 	# outside boundaries edges/faces 
# meshbox.n_corners 		# points
# meshbox.n_elements		# elements  
# meshbox.n_surfaces		# total faces 
# meshbox.n_interfaces 	# inner faces
#
# meshbox.corners = (x,y) values
#
# meshbox.element_node_ids = node (points/vertices) index
#							 for each element counterclockwise 
#
# meshbox.neighbour_information = connectivity
#
# meshbox.boundary_names = names of the boundaries 
#						(same size as element_node_ids)
##########################################################
##########################################################

function connectivity_matrix(meshbox)
	"""
		This function creates the dof mapping from local to
		global indices

		In each element el: dof_id[:,:,el] = global node numberings
		of that element
	"""
	polydeg = meshbox.polydeg
	nnodes = polydeg + 1
	nel = meshbox.n_elements
	#nsurf = meshbox.n_surfaces

	dof_id = zeros(Int, nnodes, nnodes, nel)
	#ig_surf = zeros(Int, nnodes, nsurf)

	# 1st element: all boundaries are 0 (5x5 matrix)
	ig1 = reshape(collect(1:nnodes*nnodes), nnodes, nnodes)
	
	last_dof_id = 0 

	for el in 1:nel
		if el == 1
			dof_id[:,:,el] .= ig1
		end

		# Find the surfaces for the current element to check for neighbours
		surf_id1 = findall(meshbox.neighbour_information[3,:] .== el)
		surf_id2 = findall(meshbox.neighbour_information[4,:] .== el)
		surf_id = sort(vcat(surf_id1, surf_id2))
		
		#-----------------
		# Match boundaries
		#----------------- 
		for surf in surf_id
			el1 = meshbox.neighbour_information[3,surf]
			el2 = meshbox.neighbour_information[4,surf]
			surf1 = meshbox.neighbour_information[5,surf]
			surf2 = meshbox.neighbour_information[6,surf]	
			
			if el2 == 0 || surf2 == 0
				continue
			end
			temp_var1 = 0; temp_var2 = 0
			if any(dof_id[:,:,el1] .!= 0)
				temp_var1 = el1
				el1 = el2
				el2 = temp_var1

				temp_var2 = surf1
				surf1 = surf2
				surf2 = temp_var2 
			end

				if surf1 == 1 && surf2 == 1
					dof_id[1,:,el1] = dof_id[1,:,el2]
				elseif surf1 == 1 && surf2 == 2
					dof_id[1,:,el1] = dof_id[:,nnodes,el2] 
				elseif surf1 == 1 && surf2 == 3
					dof_id[1,:,el1] = dof_id[nnodes,:,el2] 
				elseif surf1 == 1 && surf2 == 4
					dof_id[1,:,el1] = dof_id[:,1,el2] 
				elseif surf1 == 1 && surf2 == -1
					dof_id[1,:,el1] = dof_id[1,end:-1:1,el2] 
				elseif surf1 == 1 && surf2 == -2
					dof_id[1,:,el1] = dof_id[end:-1:1,nnodes,el2]
				elseif surf1 == 1 && surf2 == -3
					dof_id[1,:,el1] = dof_id[nnodes,end:-1:1,el2]
				elseif surf1 == 1 && surf2 == -4
					dof_id[1,:,el1] = dof_id[end:-1:1,1,el2]

				elseif surf1 == 2 && surf2 == 1
					dof_id[:,nnodes,el1] = dof_id[1,:,el2]  
				elseif surf1 == 2 && surf2 == 2
					dof_id[:,nnodes,el1] = dof_id[:,nnodes,el2]   
				elseif surf1 == 2 && surf2 == 3
					dof_id[:,nnodes,el1] = dof_id[nnodes,:,el2]  
				elseif surf1 == 2 && surf2 == 4
					dof_id[:,nnodes,el1] = dof_id[:,1,el2]    
				elseif surf1 == 2 && surf2 == -1
					dof_id[:,nnodes,el1] = dof_id[1,end:-1:1,el2]
				elseif surf1 == 2 && surf2 == -2
					dof_id[:,nnodes,el1] = dof_id[end:-1:1,nnodes,el2]
				elseif surf1 == 2 && surf2 == -3
					dof_id[:,nnodes,el1] = dof_id[nnodes,end:-1:1,el2]
				elseif surf1 == 2 && surf2 == -4
					dof_id[:,nnodes,el1] = dof_id[end:-1:1,1,el2]

				elseif surf1 == 3 && surf2 == 1
					dof_id[nnodes,:,el1] = dof_id[1,:,el2]   
				elseif surf1 == 3 && surf2 == 2
					dof_id[nnodes,:,el1] = dof_id[:,nnodes,el2]  
				elseif surf1 == 3 && surf2 == 3
					dof_id[nnodes,:,el1] = dof_id[nnodes,:,el2]    
				elseif surf1 == 3 && surf2 == 4
					dof_id[nnodes,:,el1] = dof_id[:,1,el2]
				elseif surf1 == 3 && surf2 == -1
					dof_id[nnodes,:,el1] = dof_id[1,end:-1:1,el2]
				elseif surf1 == 3 && surf2 == -2
					dof_id[nnodes,:,el1] = dof_id[end:-1:1,nnodes,el2]
				elseif surf1 == 3 && surf2 == -3
					dof_id[nnodes,:,el1] = dof_id[nnodes,end:-1:1,el2]
				elseif surf1 == 3 && surf2 == -4
					dof_id[nnodes,:,el1] = dof_id[end:-1:1,1,el2]
				
				elseif surf1 == 4 && surf2 == 1
					dof_id[:,1,el1] = dof_id[1,:,el2]   
				elseif surf1 == 4 && surf2 == 2
					dof_id[:,1,el1] = dof_id[:,nnodes,el2]  
				elseif surf1 == 4 && surf2 == 3
					dof_id[:,1,el1] = dof_id[nnodes,:,el2]    
				elseif surf1 == 4 && surf2 == 4
					dof_id[:,1,el1] = dof_id[:,1,el2]
				elseif surf1 == 4 && surf2 == -1
					dof_id[:,1,el1] = dof_id[1,end:-1:1,el2]
				elseif surf1 == 4 && surf2 == -2
					dof_id[:,1,el1] = dof_id[end:-1:1,nnodes,el2]
				elseif surf1 == 4 && surf2 == -3
					dof_id[:,1,el1] = dof_id[nnodes,end:-1:1,el2]
				elseif surf1 == 4 && surf2 == -4
					dof_id[:,1,el1] = dof_id[end:-1:1,1,el2]
				end # if boundary condition loop
			#end # if loop 
		end # for loop 

		#--------------------------
		# Set the rest of the nodes 
		#--------------------------
		# Find the zeros in each element
		indx = findall(dof_id[:,:,el] .== 0)

		for i in indx
			dof_id[i,el] = last_dof_id + 1
			last_dof_id = last_dof_id + 1
		end
		last_dof_id = dof_id[nnodes,nnodes,el]
	end

	dof_id 
end

