########################################################
########################################################
#
# Structured Quadrilateral Mesh 
# You can visualize the *.t file in paraview
#
#
# add documentation here ----Meshing----
# Mesh details description: 
# https://trixi-framework.github.io/HOHQMesh/the-control-input/
# https://trixi-framework.github.io/Trixi.jl/stable/meshes/unstructured_quad_mesh/
#
########################################################
########################################################

########### MESH EXAMPLE ###############################

#			Free Surface 
#	++++++++++++++++++++++++++++++++++(40e3,40e3)
#	|								 :
#	|								 :
#	|								 :  Fault Boundary 
#	|								 :
#	|								 :
#	|								(40e3, 20e3)
#	|								 }
#	|								 }   Creep Boundary
#	|								 }
#	|								 }
# (0,0)-------------------------------(40e3,0)
# 	   Absorbing/Periodic boundaries: left and bottom 

# 	Throughout the scripts, creep and creep-id refers 
#	to the constant velocity plate loading below the fault (Creep boundary). 
# 	This is not aseismic creep due to friction. Find a better name for this later.

##########################################################

using HOHQMesh

################################################
# Interactive script in julia to generate mesh 
#
# Ref: https://trixi-framework.github.io/HOHQMesh.jl/stable/interactive_overview/
#
#################################################

p = newProject("structured", mesh_path)

setPolynomialOrder!(p, 4)
setPlotFileFormat!(p, "sem")

# Mesh properties 
# For this 2d case, the z-direction is always zero
lower_left = [0.0, 0.0, 0.0]
spacing = [1000.0, 1000.0, 0.0]
num_intervals = [40, 40, 0]

addBackgroundGrid!(p, lower_left, spacing, num_intervals)
generate_mesh(p)


#= Testing unstructured 
p2 = newProject("unstructured", mesh_path)
setPolynomialOrder!(p2, 4)
setPlotFileFormat!(p2, "sem")

absorbing_bottom = new("absorbing", [0.0,0.0,0.0], [40.0e3,0.0,0.0])
creep_right = new("creep", [40.0e3,0.0,0.0], [40.0e3,20.0e3,0.0])
fault_right = new("fault", [40.0e3,20.0e3,0.0], [40.0e3,40.0e3,0.0])
free_top = new("free", [40.0e3,40.0e3,0.0], [0.0,40.0e3,0.0])
absorbing_left = new("absorbing", [0.0,40.0e3,0.0], [0.0,0.0,0.0])


# Assemble the five boundaries in a closed chain oriented couter-clockwise. The chain
# name `unstructured` is only used internally by HOHQMesh.
addCurveToOuterBoundary!(p2, absorbing_bottom)
addCurveToOuterBoundary!(p2, creep_right)
addCurveToOuterBoundary!(p2, fault_right)
addCurveToOuterBoundary!(p2, free_top)
addCurveToOuterBoundary!(p2, absorbing_left)

addBackgroundGrid!(p2, [1.0e3,1.0e3,0.0])
generate_mesh(p2)
=#

# To Plot the mesh (Include GLMakie)
function mesh_plot(proj)
	plotProject!(proj, MODEL+GRID)
   	@info "Press enter to generate the mesh and update the plot."
   	readline()
	current_figure()
end