
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
#	+++++++++++++++++++++++++++++++++++++++(40e3,40e3)
#	|					:			      | 
#	|					:			      |
#	|					: Fault Boundary  | 
#	|					:			      |
#	|					:			      |
#	|				(20e3, 20e3)          |
#	|					}		          |
#	|					} Creep Boundary  |
#	|					}		          |
#	|					}		          |
# (0,0)-------------------------------------(40e3,0)
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
	
	open(joinpath(mesh_path, "split_node.control"), "w") do io
		println(io, raw"""
        \begin{MODEL}
            \begin{OUTER_BOUNDARY}
                \begin{END_POINTS_LINE}
                    name = right
                    xEnd = [40000.0,40000.0,0.0]
                    xStart = [40000.0,0.0,0.0]
                \end{END_POINTS_LINE}
                \begin{END_POINTS_LINE}
                    name = top
                    xEnd = [0.0,40000.0,0.0]
                    xStart = [40000.0,40000.0,0.0]
                \end{END_POINTS_LINE}
                \begin{END_POINTS_LINE}
                    name = left
                    xEnd = [0.0,40000.0,0.0]
                    xStart = [0.0,0.0,0.0]
                \end{END_POINTS_LINE}
                \begin{END_POINTS_LINE}
                    name = bottom
                    xEnd = [0.0,-5.0,0.0]
                    xStart = [-10.0,-5.0,0.0]
                \end{END_POINTS_LINE}
                \begin{END_POINTS_LINE}
                    name = spike
                    xEnd = [0.0,-1.0,0.0]
                    xStart = [0.0,-5.0,0.0]
                \end{END_POINTS_LINE}
                \begin{END_POINTS_LINE}
                    name = spike
                    xEnd = [1.0e-13,-5.0,0.0]
                    xStart = [0.0,-1.0,0.0]
                \end{END_POINTS_LINE}
                \begin{END_POINTS_LINE}
                    name = bottom
                    xEnd = [10.0,-5.0,0.0]
                    xStart = [1.0e-13,-5.0,0.0]
                \end{END_POINTS_LINE}
            \end{OUTER_BOUNDARY}
        \end{MODEL}

		\begin{CONTROL_INPUT}
			\begin{RUN_PARAMETERS}
				mesh file name   = split_node.mesh
				plot file name   = split_node.tec
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
control_file = joinpath(mesh_path, "split_node.control")
generate_mesh(control_file, output_directory=mesh_path);