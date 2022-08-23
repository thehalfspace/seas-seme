###################################################################### 
# Scripts to calculate the mass and stiffness matrix for an element
# Since mass matrix is diagonal, we assemble it into a global matrix
#
#	These equations are from SEM Notes by Ampuero:
# 	
###################################################################### 

###########################
##### TODO: ###############
###########################
# 1. Check the stiffness matrix for unstructured meshes: the assembly might not be correct 
# 2. Try without assembling the mass and stiffness matrix 


function elemental_mass_matrix(rho, jac, basis)
	#x = basis.nodes
	w = basis.weights
	Me = (w*w') .* rho .* jac
	return Me
end

function global_mass_matrix!(Mglob, mesh, Me, dof_id)
	for el in 1:mesh.n_elements
		local_index = dof_id[:,:,el]
		Mglob[local_index] .+= Me[:,:,el]
	end
	Mglob 
end

function elemental_stiffness_matrix(mu, jac_mat, jac_det, basis)
	x = basis.nodes
	w = basis.weights
	H = basis.derivative_matrix'

	w2 = w*w'

	nnodes, = size(x)

	Wξξ = zeros(nnodes, nnodes)
	Wξη = zeros(nnodes, nnodes)
	Wηξ = zeros(nnodes, nnodes)
	Wηη = zeros(nnodes, nnodes)

	for i in 1:nnodes
		for j in 1:nnodes
			δx_δξ = jac_mat[1,1, i, j]
			δx_δη = jac_mat[1,2, i, j]
			δy_δξ = jac_mat[2,1, i, j]
			δy_δη = jac_mat[2,2, i, j]
			
			#inv_jac_mat = inv(jac_mat[:,:,i,j])
			#δξ_δx = inv_jac_mat[1,1]
			#δη_δx = inv_jac_mat[1,2]
			#δξ_δy = inv_jac_mat[2,1]
			#δη_δy = inv_jac_mat[2,2]

			Wξξ[i,j] = w2[i,j] * mu * jac_det[i,j] / (δx_δξ*δx_δξ + δy_δξ*δy_δξ)
			Wξη[i,j] = w2[i,j] * mu * jac_det[i,j] / (δx_δξ*δx_δη + δy_δξ*δy_δη)
			Wηξ[i,j] = w2[i,j] * mu * jac_det[i,j] / (δx_δη*δx_δξ + δy_δη*δy_δξ)
			Wηη[i,j] = w2[i,j] * mu * jac_det[i,j] / (δx_δη*δx_δη + δy_δη*δy_δη)

			# Check for divide by zero 
			if isnan(Wξξ[i,j]) || isinf(Wξξ[i,j])
				Wξξ[i,j] = 0.
			end
			if isnan(Wξη[i,j]) || isinf(Wξη[i,j])
				Wξη[i,j] = 0.
			end
			if isnan(Wηξ[i,j]) || isinf(Wηξ[i,j])
				Wηξ[i,j] = 0.
			end
			if isnan(Wηη[i,j]) || isinf(Wηη[i,j])
				Wηη[i,j] = 0.
			end

		end
	end

	Ke_temp::Array{Float64,4} = zeros(nnodes, nnodes, nnodes, nnodes)	
	Ke::Matrix{Float64} = zeros(size(H))

	δ = Matrix(I, size(H))
	for i in 1:nnodes, j in 1:nnodes
		for k in 1:nnodes, l in 1:nnodes
			Ke_ξξ = 0.; Ke_ξη = 0.; Ke_ηξ = 0.; Ke_ηη = 0.
			for p in 1:nnodes
				for q in 1:nnodes
					# Double check these equations
					Ke_ξξ += Wξξ[p,q]*δ[j,q]*δ[l,q]*H[i,p]*H[k,p]
					Ke_ξη += Wξη[p,q]*δ[j,q]*δ[k,p]*H[i,p]*H[l,q]
					Ke_ηξ += Wηξ[p,q]*δ[i,p]*δ[l,q]*H[j,q]*H[k,p]
					Ke_ηη += Wηη[p,q]*δ[i,p]*δ[k,p]*H[j,q]*H[l,q]
				end
				
				# With Kronecker reduction (You can reduce a loop for cartesian mesh)
				#Ke_ξξ += δ[j,l] * Wξξ[p,j] * H[i,p] * H[k,p]
				#Ke_ηη += δ[i,k] * Wηη[k,p] * H[j,p] * H[l,p]
			end
			Ke_temp[i,j,k,l] = Ke_ξξ +  Ke_ηη + Ke_ξη + Ke_ηξ
		end
	end
	
	Ke = reshape(Ke_temp, nnodes*nnodes, nnodes*nnodes)
	return Ke
end