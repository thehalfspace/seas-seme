########################################################
########################################################
#
# Unstructured Quadrilateral Mesh 
# You can visualize the *.t file in paraview
#
#
# add documentation here ----Meshing----
# Mesh details description: 
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
# 	   Absorbing/Periodic boundaries

##########################################################


struct mesh_params
	Lx::Float64
	Ly::Float64
	
	dx::Float64
	dy::Float64

	fault_p1::Tuple
	fault_p2::Tuple
end

##= Think about how to make the workflow better

# Make the outpath same here and in time_loop

out_path = joinpath(dirname(pwd()), "spear-parallel/data/test_01/")
mesh_path = joinpath(out_path, "mesh/")
if isdir(out_path)
	if isdir(mesh_path)
		rm(mesh_path, recursive=true)
	end
	mkpath(mesh_path)
else
	mkpath(out_path)
	mkpath(mesh_path)
end
# =#

# Try the entire thing with a structured mesh first and cross check the results

# Create control file
function create_mesh_control()
	
	open(joinpath(mesh_path, "unstructured.control"), "w") do io
		println(io, raw"""
		\begin{CONTROL_INPUT}
			\begin{RUN_PARAMETERS}
				mesh file name   = unstructured.mesh
				plot file name   = unstructured.tec
				stats file name  = none
				mesh file format = ISM-v2
				polynomial order = 4
				plot file format = skeleton
			\end{RUN_PARAMETERS}

			\begin{BACKGROUND_GRID}
				background grid size = [2.0E3,2.0E3,0.0]
			\end{BACKGROUND_GRID}

			\begin{SPRING_SMOOTHER}
    			smoothing            = ON
    			smoothing type       = LinearAndCrossBarSpring
				number of iterations = 50
			\end{SPRING_SMOOTHER}
			\begin{REFINEMENT_REGIONS}
			\begin{REFINEMENT_LINE}
    			type = smooth
    			x0   = [40.0E3,20.0E3,0.0]
    			x1   = [40.0E3,40.0E3,0.0]
    			h    = 0.2E3
    			w    = 2.0E3
			\end{REFINEMENT_LINE}
			\end{REFINEMENT_REGIONS}
		\end{CONTROL_INPUT}
		
		\begin{MODEL}
			\begin{OUTER_BOUNDARY}
			\begin{END_POINTS_LINE}
				name = absorbing 
				xStart = [0.0,0.0,0.0]
				xEnd   = [40.0E3,0.0,0.0]
			\end{END_POINTS_LINE}
			\begin{END_POINTS_LINE}
				name = creep
				xStart = [40.0E3,0.0,0.0]
				xEnd   = [40.0E3,20.0E3,0.0]
			\end{END_POINTS_LINE}
			\begin{END_POINTS_LINE}
				name = fault
				xStart = [40.0E3,20.0E3,0.0]
				xEnd   = [40.0E3,40.0E3,0.0]
			\end{END_POINTS_LINE}
			\begin{END_POINTS_LINE}
				name = free_surface
				xStart = [40.0E3,40.0E3,0.0]
				xEnd   = [0.0,40.0E3,0.0]
			\end{END_POINTS_LINE}
			\begin{END_POINTS_LINE}
				name = absorbing 
				xStart = [0.0,40.0E3,0.0]
				xEnd   = [0.0,0.0,0.0]
			\end{END_POINTS_LINE}
			\end{OUTER_BOUNDARY}
		\end{MODEL}
		\end{FILE}
		""")
	end
end

using HOHQMesh
create_mesh_control()
control_file = joinpath(mesh_path, "unstructured.control")
generate_mesh(control_file, output_directory=mesh_path);

#=
\begin{BACKGROUND_GRID}
				x0 = [0.0, 0.0, 0.0]
				dx = [1.0, 1.0, 0.0]
				N  = [4,4,1]
			\end{BACKGROUND_GRID}

			\begin{BACKGROUND_GRID}
    			background grid size = [1.0,1.0,0.0]
			\end{BACKGROUND_GRID}
			
=#