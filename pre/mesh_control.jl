
# Create control file (this is writing a fortran script: scroll down for interactive julia script)
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
# Uncomment these three lines to use control file from the above fortran script
create_mesh_control()
control_file = joinpath(mesh_path, "structured.control")
generate_mesh(control_file, output_directory=mesh_path);






##############################################
###### UNSTRUCTURED ##########################
##############################################

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
				background grid size = [0.5E3,0.5E3,0.0]
			\end{BACKGROUND_GRID}

			\begin{SPRING_SMOOTHER}
    			smoothing            = ON
    			smoothing type       = LinearAndCrossBarSpring
				number of iterations = 50
			\end{SPRING_SMOOTHER}
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
			
			\begin{REFINEMENT_REGIONS}
			\begin{REFINEMENT_LINE}
    			type = smooth
    			x0   = [40.0E3,20.0E3,0.0]
    			x1   = [40.0E3,40.0E3,0.0]
    			h    = 0.2E3
    			w    = 2.0E3
			\end{REFINEMENT_LINE}
			\end{REFINEMENT_REGIONS}
=#

