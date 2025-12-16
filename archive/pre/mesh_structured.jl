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

# Remove mesh_path if it already exists 
if isdir(out_path)
	if isdir(mesh_path)
		rm(mesh_path, recursive=true)
	end
	mkpath(mesh_path)
end

# Create control file
function create_mesh_control()
	
	open(joinpath(mesh_path, "structured.control"), "w") do io
		println(io, raw"""
		\begin{CONTROL_INPUT}
			\begin{RUN_PARAMETERS}
				mesh file name   = structured.me
				plot file name   = structured.te
				stats file name  = none
				mesh file format = ISM-v2
				polynomial order = 4
				plot file format = skeleton
			\end{RUN_PARAMETERS}

			\begin{BACKGROUND_GRID}
				x0 = [0.0, 0.0, 0.0]
    			dx = [1000.0, 1000.0, 0.0]
    			N  = [40,40,0]
			\end{BACKGROUND_GRID}

			\begin{SPRING_SMOOTHER}
    			smoothing            = ON
    			smoothing type       = LinearAndCrossBarSpring
				number of iterations = 50
			\end{SPRING_SMOOTHER}
		\end{CONTROL_INPUT}
		\end{FILE}
		""")
	end
end

using HOHQMesh
create_mesh_control()
control_file = joinpath(mesh_path, "structured.control")
generate_mesh(control_file, output_directory=mesh_path);