using Trixi, LinearAlgebra, DelimitedFiles, SparseArrays, Printf, StatsBase


#----------------------------#
# Set simulation name here 
#----------------------------#
simulation_name = "test_01/"
folder_name = "seas-seme/"
#----------------------------#

# Create paths to save stuff  
pre_path = joinpath(dirname(pwd()), folder_name, "pre/")
proc_path = joinpath(dirname(pwd()), folder_name, "proc/")
out_path = joinpath(dirname(pwd()), folder_name, "data/", simulation_name)
mesh_path = joinpath(out_path, "mesh/")


# =
# Comment this block if mesh is already generated
##############################################################
##############################################################
# Clean older files and paths for simulation output
if isdir(out_path)
	rm(out_path, recursive=true)
end 

mkpath(out_path)
mkpath(mesh_path)

#-------------------------------------------------#
# 					Create Mesh					  #
#-------------------------------------------------#
# Structured mesh is cartesian grid: ignore 
# include(joinpath(pre_path, "mesh_structured.jl"))
include(joinpath(pre_path, "mesh_split_node.jl"))
#-------------------------------------------------#

# Note: structured cartesian mesh will give warning: MoDEL block missing
#       you can safely ignore that
##############################################################
##############################################################
# =#

# Include other functions
include(joinpath(pre_path, "element_connectivity.jl"))
include(joinpath(pre_path, "set_paramaters.jl"))
include(joinpath(pre_path, "set_initial_conditions.jl"))
include(joinpath(pre_path, "local_mass_stiffness.jl"))
include(joinpath(pre_path, "boundary_information.jl"))
include(joinpath(pre_path, "output_locations.jl"))
include(joinpath(pre_path, "mappings_geometry_straight_2d.jl"))
include(joinpath(proc_path, "element_calculations.jl"))
include(joinpath(proc_path, "compute_timestep.jl"))
include(joinpath(proc_path, "fault_solvers.jl"))
include(joinpath(proc_path, "time_loop_structured.jl"))


#-------------------------------------------------#
# 					Time Loop					  #
#-------------------------------------------------#
include(joinpath(proc_path, "time_loop_structured.jl"))
#-------------------------------------------------#
@elapsed @time main_loop();

@info("\n!!!!!!!Yay Simulation Complete!!!!!!!")