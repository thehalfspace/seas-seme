using PyPlot
#using StatsBase
#using LaTeXStrings
#using PyCall
#mpl = pyimport("matplotlib")

# Default plot params
function plot_params()
  plt.rc("xtick", labelsize=16)
  plt.rc("ytick", labelsize=16)
  plt.rc("xtick", direction="in")
  plt.rc("ytick", direction="in")
  plt.rc("font", size=15)
  plt.rc("figure", autolayout="True")
  plt.rc("axes", titlesize=16)
  plt.rc("axes", labelsize=17)
  plt.rc("xtick.major", width=1.5)
  plt.rc("xtick.major", size=5)
  plt.rc("ytick.major", width=1.5)
  plt.rc("ytick.major", size=5)
  plt.rc("lines", linewidth=2)
  plt.rc("axes", linewidth=1.5)
  plt.rc("legend", fontsize=13)
  plt.rc("mathtext", fontset="stix")
  plt.rc("font", family="STIXGeneral")

  # Default width for Nature is 7.2 inches, 
  # height can be anything
  #plt.rc("figure", figsize=(7.2, 4.5))
end



#mesh_file = joinpath(mesh_path, "structured.mes")
#m = UnstructuredMesh2D(mesh_file)


#=
# Derivative matrix for a reference element: 
# H[i,j] = l'_i(ξ_j)	(where l_i(ξ) is the shape function)
H = [-5.0 6.7565 -2.66667 1.41016 -0.5;
	 -1.24099 0.0 1.74574 -0.763763 0.25901;
	 0.375 -1.33658 0.0 1.33658 -0.375;
	 -0.25901 0.763763 -1.74574 0.0 1.24099;
	 0.5 -1.41016 2.66667 -6.7565 5.0]

# This is the material property:
# W[i,j] = shear_modulu*w[i]*w[j] (where w is the weights for quadrature)
W = [3.20381e8 1.7443e9 2.27827e9 1.7443e9 3.20381e8;
	 1.7443e9 9.49673e9 1.24039e10 9.49673e9 1.7443e9;
	 2.27827e9 1.24039e10 1.6201e10 1.24039e10 2.27827e9;
	 1.7443e9 9.49673e9 1.24039e10 9.49673e9 1.7443e9;
	 3.20381e8 1.7443e9 2.27827e9 1.7443e9 3.20381e8]

function temp2(Wξξ, H)
	nnodes = 5
	Ke_temp::Array{Float64,4} = zeros(nnodes, nnodes, nnodes, nnodes)	
	Ke_temp2::Array{Float64,4} = zeros(nnodes, nnodes, nnodes, nnodes)	
	δ = Matrix(I, size(H))
	for i in 1:nnodes, j in 1:nnodes
		for k in 1:nnodes, l in 1:nnodes
			Ke_ξξ = 0.; Ke_2 = 0.
			for p in 1:nnodes 
				for q in 1:nnodes
					Ke_ξξ +=  Wξξ[p,q]*δ[j,q]*δ[l,q]*H[i,p] * H[k,p]
				end
				Ke_2 += δ[j,l] * Wξξ[p,j] * H[i,p] * H[k,p]		# This gives me correct results 
			end
			Ke_temp[i,j,k,l] = Ke_ξξ
			Ke_temp2[i,j,k,l] = Ke_2
		end
	end
	reshape(Ke_temp, nnodes*nnodes, nnodes*nnodes), reshape(Ke_temp2, nnodes*nnodes, nnodes*nnodes)
end
=#
