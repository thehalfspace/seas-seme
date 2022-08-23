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


