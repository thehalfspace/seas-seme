#####################################################################
#####################################################################
#
#	These scripts return the boundary information for the mesh:
#	Returns:
#		boundary_node_id: this gives the global numbering for the 
#						  boundary nodes.
#		boundary_x: global x position for each boundary node.
#		boundary_y: global y position for each boundary node.
#		boundary_mat: gives the contribution from each boundary
#					  node to the global system.
#
#	Equations from: Ampuero SEM Notes page 23 to calculate jac1D at boundaries	
#####################################################################
#####################################################################
# Example: for one quad element in the mesh:
#          boundary faces are numbered 1,2,3,4 counterclockwise from the bottom
#          called surfaces or interfaces
#
#	ELEMENT BOUNDARIES:
#		Corner Points: P1, P2, P3, P4 
#		Boundary surfaces: 1, 2, 3, 4
#     P4-------------------P3
#      |       (3)         |
#      |                   |
#      |(4)             (2)|
#      |                   |
#      |                   |
#      |       (1)         |
#     P1-------------------P2
#
#	NODE NUMBERINGS: dof_id[:,:,element_id]
#		Example: for first element, the global nodes are below
#				 stored in dof_id[:,:,1]
#     P4-------------------P3
#      | 5  10`15  20  25  |
#      | 4  9  14  19  24  |
#      | 3  8  13  18  23  |
#      | 2  7  12  17  22  |
#      | 1  6  11  16  21  |
#     P1-------------------P2
#
#
# Check pre/mesh_structured.jl and HOHQMesh documentation for more
# details about mesh numbering convention
#####################################################################
#####################################################################


# For structured cartesian mesh
function get_boundary_nodes_structured(mesh, node_coords, jac_matrix, weights, impedance, dof_id, boundary_name)
	""" 
	Get global positions of specified boundaries
	Only works for boundaries in the same coordinate system as (ξ,η):
		i.e. horizontal and vertical boundaries
	"""
	polydeg = mesh.polydeg 
	nnodes = polydeg + 1
	#fault_el = []
	if boundary_name == :absorbing
		# Left and Bottom boundaries are absorbing 
		boundary_el_id = unique(vcat(findall(mesh.boundary_names .== :Left), 
									findall(mesh.boundary_names .== :Bottom)))
	
	elseif boundary_name == :free_surface
		boundary_el_id = findall(mesh.boundary_names .== :Top)
	
	elseif boundary_name == :fault
		id_temp = findall(mesh.boundary_names .== :Right)
		boundary_el_id = id_temp[Int(length(id_temp)/2)+1:end]
		#fault_el = copy(boundary_el_id)
	
	elseif boundary_name == :creep
		id_temp = findall(mesh.boundary_names .== :Right)
		boundary_el_id = id_temp[1:Int(length(id_temp)/2)]
	end

	boundary_node_id = Int[]
	boundary_x = Float64[]
	boundary_y = Float64[]
	boundary_mat = Float64[]

	boundary_mat_local = zeros(nnodes)

	for id in boundary_el_id
		surface = id[1]
		el = id[2]

		ncx = node_coords[1,:,:,el]
		ncy = node_coords[2,:,:,el] 

		jac = jac_matrix[:,:,:,:,el]
		
		# Ampuero SEM Notes page 23 to calculate jac1D at boundaries
		# **TODO: Check if this works correctly for nonuniform mesh 
		if surface == 1		# Bottom: (dx/dξ)^2 + (dy/dξ)^2  	# bottom index [1,:]
			jac1D = sqrt.(jac[1,1,1,:].^2 .+ jac[2,1,1,:].^2)
		
		elseif surface == 2	# Right: (dx/dη)^2 + (dy/dη)^2		# right index [:,end]
			jac1D = sqrt.(jac[1,2,:,end].^2 .+ jac[2,2,:,end].^2)
		
		elseif surface == 3	# Top: (dx/dξ)^2 + (dy/dξ)^2		# top index [end,:]
			jac1D = sqrt.(jac[1,1,end,:].^2 .+ jac[2,1,end,:].^2)
		
		elseif surface == 4	# Left: (dx/dη)^2 + (dy/dη)^2		# left index [:,1]
			jac1D = sqrt.(jac[1,2,:,1].^2 .+ jac[2,2,:,1].^2)
		end

		boundary_mat_local .= weights .*jac1D .*impedance 

		append!(boundary_node_id, get_element_surface(dof_id[:,:,el], surface))
		append!(boundary_x, get_element_surface(ncx, surface))
		append!(boundary_y, get_element_surface(ncy, surface))
		append!(boundary_mat, boundary_mat_local) #get_element_surface(boundary_mat_local, surface))
	end

	# Add edge contributions
	# TODO: Possibly rewrite this to make it more efficient 
	boundary_mat_unique = []
	boundary_x_unique = []
	boundary_y_unique = []
	boundary_node_id_unique = []
	boundary_node_id_unique = unique(boundary_node_id)
	i = 0
	for n in boundary_node_id_unique
		i = i+1
		id = findall(boundary_node_id .== n)
		push!(boundary_mat_unique, sum(boundary_mat[id]))
		push!(boundary_x_unique, boundary_x[id[1]])
		push!(boundary_y_unique, boundary_y[id[1]])
		#push!(boundary_node_id_unique, boundary_node_id[id[1]])
	end


	return boundary_node_id_unique, boundary_x_unique, boundary_y_unique, boundary_mat_unique, boundary_el_id
	#return unique(boundary_node_id_unique), unique(boundary_x_unique), unique(boundary_y_unique), unique(boundary_mat_unique)
end

# For unstructured mesh with straight boundaries 
function get_boundary_nodes_unstructured(mesh, node_coords, jac_matrix, weights, impedance, dof_id, boundary_name)
	""" 
	Get global positions of specified boundaries
	WARNING: Right now the boundaries must be straight and same directions as (ξ,η) 
	"""
	polydeg = mesh.polydeg 
	nnodes = polydeg + 1
	boundary_el_id = findall(mesh.boundary_names .== boundary_name)

	boundary_node_id = Int[]
	boundary_x = Float64[]
	boundary_y = Float64[]
	boundary_mat = Float64[]

	boundary_mat_local = zeros(nnodes)

	for id in boundary_el_id
		surface = id[1]
		el = id[2]

		ncx = node_coords[1,:,:,el]
		ncy = node_coords[2,:,:,el] 

		jac = jac_matrix[:,:,:,:,el]
		
		# Ampuero SEM Notes page 23 to calculate jac1D at boundaries
		# **TODO: Check if this works correctly for nonuniform mesh 
		if surface == 1		# Bottom: (dx/dξ)^2 + (dy/dξ)^2  	# bottom index [1,:]
			jac1D = sqrt.(jac[1,1,1,:].^2 .+ jac[2,1,1,:].^2)
		
		elseif surface == 2	# Right: (dx/dη)^2 + (dy/dη)^2		# right index [:,end]
			jac1D = sqrt.(jac[1,2,:,end].^2 .+ jac[2,2,:,end].^2)
		
		elseif surface == 3	# Top: (dx/dξ)^2 + (dy/dξ)^2		# top index [end,:]
			jac1D = sqrt.(jac[1,1,end,:].^2 .+ jac[2,1,end,:].^2)
		
		elseif surface == 4	# Left: (dx/dη)^2 + (dy/dη)^2		# left index [:,1]
			jac1D = sqrt.(jac[1,2,:,1].^2 .+ jac[2,2,:,1].^2)
		end

		boundary_mat_local .= weights .*jac1D .*impedance 

		append!(boundary_node_id, get_element_surface(dof_id[:,:,el], surface))
		append!(boundary_x, get_element_surface(ncx, surface))
		append!(boundary_y, get_element_surface(ncy, surface))
		append!(boundary_mat, boundary_mat_local) #get_element_surface(boundary_mat_local, surface))
	end

	# Add edge contributions
	# TODO: Possibly rewrite this to make it more efficient 
	boundary_mat_unique = []
	boundary_x_unique = []
	boundary_y_unique = []
	boundary_node_id_unique = []
	boundary_node_id_unique = unique(boundary_node_id)
	i = 0
	for n in boundary_node_id_unique
		i = i+1
		id = findall(boundary_node_id .== n)
		push!(boundary_mat_unique, sum(boundary_mat[id]))
		push!(boundary_x_unique, boundary_x[id[1]])
		push!(boundary_y_unique, boundary_y[id[1]])
		#push!(boundary_node_id_unique, boundary_node_id[id[1]])
	end

	return boundary_node_id_unique, boundary_x_unique, boundary_y_unique, boundary_mat_unique, boundary_el_id
end

# Get the global node numberings depending on the surface
function get_element_surface(dof_id::Array{Int,2}, surface)

	ig_boundary = Int[]

	if surface == 1		# bottom
		ig_boundary = dof_id[1,:]
	elseif surface == 2	# right
		ig_boundary = dof_id[:,end]
	elseif surface == 3	# top 
		ig_boundary = dof_id[end,:]
	elseif surface == 4 # left 
		ig_boundary = dof_id[:,1]
	elseif surface == -1 # bottom reverse
		ig_boundary = dof_id[1,end:-1:1]
	elseif surface == -2 # right reverse 
		ig_boundary = dof_id[end:-1:1,end]
	elseif surface == -3 # top reverse 
		ig_boundary = dof_id[end,end:-1:1]
	elseif surface == -4 # left reverse 
		ig_boundary = dof_id[end:-1:1,1]
	else
		@error("Check surface id values")
	end

	ig_boundary 
end

# Get the global x positions 
function get_element_surface(dof_id::Array{Float64,2}, surface)

	ig_boundary = Float64[]
	
	dof_id = dof_id'

	if surface == 1		# bottom
		ig_boundary = dof_id[1,:]
	elseif surface == 2	# right
		ig_boundary = dof_id[:,end]
	elseif surface == 3	# top 
		ig_boundary = dof_id[end,:]
	elseif surface == 4 # left 
		ig_boundary = dof_id[:,1]
	elseif surface == -1 # bottom reverse
		ig_boundary = dof_id[1,end:-1:1]
	elseif surface == -2 # right reverse 
		ig_boundary = dof_id[end:-1:1,end]
	elseif surface == -3 # top reverse 
		ig_boundary = dof_id[end,end:-1:1]
	elseif surface == -4 # left reverse 
		ig_boundary = dof_id[end:-1:1,1]
	else
		@error("Check surface id values")
	end

	ig_boundary 
end

