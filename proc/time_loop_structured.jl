
# Rough work 
function stiffness_assembly(K_el, dof_id)
	#  Ksparse = spzeros(ndof,ndof) 
	I = Vector{Int}(undef, length(K_el))
	J = Vector{Int}(undef, length(K_el))
	V = Vector{Float64}(undef, length(K_el))
	ct = 1

	Nel = size(dof_id)[3]	# number of elements 
	ndof = maximum(dof_id)
	for eo in 1:Nel
		v = view(dof_id,:,:,eo)
		for j in 1:length(v)
		for i in 1:length(v)
			I[ct] = v[i]
			J[ct] = v[j]
			V[ct] = K_el[i, j, eo]
			ct += 1
		end
		end
	end

	return sparse(I,J,V,ndof,ndof,+)
end

# Create global variable structures
mutable struct global_vars{T<:Vector{Float64}}
	u::T
	v::T
	a::T
end

function main_loop()
	polydeg = 4
	nnodes = polydeg + 1
	basis = LobattoLegendreBasis(polydeg)

	####### MESH ##########
	mesh_file = joinpath(mesh_path, "structured.mesh")
	mesh = UnstructuredMesh2D(mesh_file)

	# Global to local index numbering
	dof_id = connectivity_matrix(mesh)

	# Number of global nodes 
	ndof = maximum(dof_id)

	# Temporary material properties for NOW
	# TODO: Add tags in mesh for different material properties 
	ρ = 2670 
	vs = 3464
	μ = ρ*vs^2

	# Get global node coordinates and the jacobian matrix 
	node_coordinates = zeros(2, nnodes, nnodes, mesh.n_elements)
	jac_matrix = zeros(2, 2, nnodes, nnodes, mesh.n_elements)
	jac_det = zeros(nnodes, nnodes, mesh.n_elements)
	jac_1d = zeros(nnodes, nnodes, mesh.n_elements)		# Check page 22 SEM notes 
	normal_dirs = zeros(2,nnodes, 4, mesh.n_elements)	# 2 directions, 4 sides in each element

	# Mass and stiffness elemental matrix: is constant for all elements for now
	M_el = zeros(nnodes, nnodes, mesh.n_elements)
	K_el = zeros(nnodes*nnodes, nnodes*nnodes, mesh.n_elements)

	for el in 1:mesh.n_elements
		corners = mesh.corners[:,mesh.element_node_ids[:,el]]'
		calc_node_coordinates!(node_coordinates, el, basis.nodes,corners)
		calc_metric_terms!(jac_matrix, el, basis.nodes, corners)
		calc_normal_directions!(normal_dirs, el, basis.nodes, corners)

		for i in 1:nnodes, j in 1:nnodes
			jac_det[i,j,el] = det(jac_matrix[:,:,i,j,el])
			jac_1d[i,j,el] = jac_matrix[1,1,i,j,el]^2 + jac_matrix[2,1,i,j,el]^2 # Check Ampuero SEM notes 
		end

		M_el[:,:,el] = elemental_mass_matrix(ρ, jac_det[:,:,el], basis)
		K_el[:,:,el] = elemental_stiffness_matrix(μ, jac_matrix[:,:,:,:,el], jac_det[:,:,el], basis) 
	end

	# Global Mass Matrix: Since it is diagonal, we can store globally 
	# TODO: don't store M_el, calculate dynamically in the global_mass_matrix function 
	Mglob = zeros(maximum(dof_id))
	global_mass_matrix!(Mglob, mesh, M_el, dof_id)

	# Get parameters and initial conditions
	#yr2sec, Total_Time, dtmax, Vpl, fo, Vo, Vthres, CFL
	par = set_default_parameters()


	#= BOUNDARY INFO FOR UNSTRUCTURED MESH (** REWRITE THIS **)
	# Boundary nodes and info  
	# Impedance = ρ*vs (on absorbing and free boundaries), = 1 on fault and creep boundary 
	fault_id, fault_x, fault_y, fault_matrix = 
		get_boundary_nodes(mesh, node_coordinates, normal_dirs, basis.weights, 1, dof_id, :fault)
	creep_id, creep_x, creep_y, creep_matrix = 
		get_boundary_nodes(mesh, node_coordinates, normal_dirs, basis.weights, 1, dof_id, :creep)
	absorbing_id, absorbing_x, absorbing_y, absorb_matrix = 
		get_boundary_nodes(mesh, node_coordinates, normal_dirs, basis.weights, ρ*vs, dof_id, :absorbing)
	#free_id, free_x, free_y, free_matrix = 
	#	get_boundary_nodes(mesh, node_coordinates, jac_det, basis.weights, ρ*vs, dof_id, :free_surface)
	=#

	## BOUNDARY INFO FOR STRUCTURED MESH: I am assuming that the right side boundary is half-fault and half-plate loading.
	# Impedance = ρ*vs (on absorbing and free boundaries), = 1 on fault and creep boundary 
	# Left and Bottom boundaries are absorbing 
	# Right boundary is fault and creep
	fault_id, fault_x, fault_y, fault_matrix = 
		get_boundary_nodes_structured(mesh, node_coordinates, jac_matrix, basis.weights, 1, dof_id, :fault)
	
	creep_id, creep_x, creep_y, creep_matrix = 
		get_boundary_nodes_structured(mesh, node_coordinates, jac_matrix, basis.weights, 1, dof_id, :creep)
	
	absorbing_id, absorbing_x, absorbing_y, absorb_matrix = 
		get_boundary_nodes_structured(mesh, node_coordinates, jac_matrix, basis.weights, ρ*vs, dof_id, :absorbing)

	# Sparse global stiffness assembly 
	Ksparse = spzeros(ndof, ndof)
	Ksparse = stiffness_assembly(K_el, dof_id)

	# Creep and fault id share one node here: 
	# add the contributions to fault matrix and remove the node from creep id.
	## Find a better workflow to tackle this later.
	fault_matrix[1] = fault_matrix[1] + creep_matrix[1]
	creep_id = creep_id[1:end-1]
	creep_matrix = creep_matrix[1:end-1]
	creep_y = creep_y[1:end-1]

	# Interface id: right side boundary: fault+creep
	# 0-20 km = creep, 20-40 km = fault 
	interface_id = vcat(creep_id, fault_id)

	# These are the global dofs not on the fault (fault + creep):
	# used to split stiffness into on-fault and off-fault parts.
	fltni = unique(dof_id[:])
	deleteat!(fltni, interface_id)
	
	# Some constants for timestep evolution calculation 
	hcell = mean(diff(fault_y))
	μmax = μ 	# Update this when adding material tags 
	
	# minimum timestep is c*dx/vs, c is constant: use between 0.1 - 0.5
	# CFL criteria is usually 0.6, but we want the minimum timestep to be a little lower
	# Check SCEC SEAS benchmark papers for discussions
	dt_min = 0.2*hcell/vs 

	# Update Mass at boundaries for neumann conditions: implicit treatment of terms
	# Absorbing boundary = inhomogeneous neumann boundary
	# Free surface = homogeneous neumann boundary: it is implicitly satisfied so nothing to do
	# Fault boundary = Robin boundary condition: solving ode on fault, so nothing to do here
	Mglob[absorbing_id] .= Mglob[absorbing_id] .+ (0.5*dt_min).*absorb_matrix 
	
	# Fault boundary inertia terms (stick-traction): Check Kaneko (2008 and 2011) for details
	fault_z::Vector{Float64} = Mglob[fault_id]./(fault_matrix.*dt_min)
	fault_vfree::Vector{Float64} = zeros(size(fault_id))

	# Set initial conditions (Fault is vertical => only along y directions)
	ics = set_initial_conditions(fault_y, out_path)

	
	# Initialize global variables
	# TODO: Distribute this using mesh coloring algorithm
	glob = global_vars(zeros(ndof), zeros(ndof), zeros(ndof))

	# Temporary variables to store previous timestep information
	u_pre = zeros(ndof)
	v_pre = zeros(ndof)

	# f = Kf.uf + Km.um (Introduce a new term from Kaneko 2011: Eq 4)
	f = zeros(ndof)	

	# On-fault variables: for frictional solver
	# Shear Stress 
	τn = zeros(size(fault_id))
	
	# State Variable 
	θo = zeros(size(fault_id))
	θn = zeros(size(fault_id))
	
	# Transformed state variable: ψ = log(θ.Vo/Lc) 
	# (as used by Kaneko's matlab code)
	ψo = zeros(size(fault_id))
	ψn = zeros(size(fault_id))

	# On-fault Slip-rate = δ_dot in Kaneko (2011)
	Vfo = zeros(size(fault_id))
	Vfn = zeros(size(fault_id)) 
	Vfnp1 = zeros(size(fault_id))

	# Initialize velocity (and fault-slip-rate) fields
	v_init = 0.5e-3
	glob.v .= v_init .- 0.5*par.Vpl  
	Vfn .= 2*glob.v[fault_id]

	#----------------SAVE OUTPUT TO FILE-----------------------------------------------
	# Output variables to save: Change output locations in pre/output_locations.jl
	out_fault_id::Vector{Int}, out_glob_id::Vector{Int}, out_loc = get_output_locations(fault_y, fault_id)

	# Open output files for writing: Time series at different locations on the fault 
	time_series = []
	for (i,o) in enumerate(out_loc)
		push!(time_series, open(joinpath(out_path, "time_series_depth_$(o)km.out"), "w"))
		write(time_series[i], "t slip sliprate shear_stress state \n")
	end

	# Save Scalar/1D parameters (e.g., peak slip-rate, shear modulus change, etc) through time
	time_series_1d = open(joinpath(out_path, "time_series_1d.out"), "w")
	write(time_series_1d, "t max_sliprate \n")
	#-----------------------------------------------------------------------------------

	# Initial state variable: Not using this right now
	θo = begin 
		(ics.Lc ./ par.Vo) .* exp.(
			(ics.a ./ ics.b) .* log.((2*par.Vo./par.Vpl) .* 
			sinh.(ics.τo ./ (ics.a .* ics.σo))) .- 
			(par.fo ./ ics.b)
		)
	end
	
	# Change variable ψ = log(θ.Vo/Lc)
	ψo = ics.τo ./(ics.σo.*ics.b) .- par.fo./ics.b .- (ics.a./ics.b).*log.(2*glob.v[fault_id]./par.Vo)
	ψn .= ψo

	glob.v[creep_id] .= 0.

	# Iterators
	it = 0
	t = 0.
	Vfmax = 0.
	dt = dt_min

	# Debug variables
	tsave::Vector{Float64} = []
	Vfmaxsave::Vector{Float64} = []

	# Display some information before the time loop 
	@printf("\nMinimum dt: %1.03f s\n", dt)
	@printf("Average node spacing along fault: %1.02f m\n", hcell)
	@printf("Simulation Time: %1.02f years\n", par.Total_Time/par.yr2sec)

	# Solver is initially quasistatic. 
	# But because of the overstressed patch and high initial slip-rates,
	# we have an artifical earthquake at time = 0.
	solver = :quasistatic

	##################################################################################
	# The steps in the time-loop are from Kaneko et al (2011) and Kaneko et al (2008)
	##################################################################################

	while t < par.Total_Time
		it = it+1
		t = t + dt 

		if solver == :quasistatic 

			# Store previous global displacement and velocity
			u_pre .= glob.u
			v_pre .= glob.v 

			# Velocity (slip-rate) on fault 
			Vfn .= 2*v_pre[fault_id] .+ par.Vpl
			Vfo .= Vfn 

			for p1 = 1:2
			
				# Step 1: Compute on-fault displacement 
				f .= 0.

				f[fault_id] .= u_pre[fault_id] .+ glob.v[fault_id]*dt
				f[creep_id] .= u_pre[creep_id] .+ glob.v[creep_id]*dt

				# Step 2: Quasi-static linear solver
				# Try the distributed linear solve algorithm here.
				# Check julia linearsolve package.
				#=
				for el in 1:mesh.n_elements
					local_index = dof_id[:,:,el]

					# Ask on julia discourse how to solve this 
					# Create an MWE 
					u_el = K_el[:,:,el]*f[local_index][:] 
					#if  iszero(u_el)
					#	u_solve .= u_el 
					#else
					u_solve .=  qr(K_el[:,:,el], Val(true)) \ u_el
					#end
					glob.u[local_index] .= glob.u[local_index] .- reshape(u_solve, nnodes, nnodes)
				end
				=#

				# Step 2: Quasi-static linear solver
				# Direct sparse solver 
				rhs = (Ksparse*f)[fltni] 
		
				glob.u[fltni] .= -(Ksparse[fltni, fltni]\rhs)

				# Make u = f on fault 
				glob.u[fault_id] .= u_pre[fault_id] .+ glob.v[fault_id]*dt
				glob.u[creep_id] .= u_pre[creep_id] .+ glob.v[creep_id]*dt

				# Step 3: Compute on-fault stresses 
				# (f* is glob.a in my case: to save space and reuse variables): a is typically zero in quasi-static algorithm
				
				# Try to do this in a distributed loop 
				#=
				for el in 1:mesh.n_elements
					local_index = dof_id[:,:,el]	
					a_el = K_el[:,:,el]*glob.u[local_index][:]
					glob.a[local_index] .= glob.a[local_index] .+ reshape(a_el, nnodes, nnodes) 
				end =#

				# Step 3: Compute on-fault traction 
				# (f* is glob.a in my case: to save space and reuse variables): a is typically zero in quasi-static algorithm
				glob.a .= 0.
				glob.a .= Ksparse*glob.u
				
				# Enforce a = 0 on creep boundary
				glob.a[creep_id] .= 0	

				# Step 3: Calculate on-fault traction  
				τn .= - glob.a[fault_id]./fault_matrix

				# Step 4: Determine first prediction of state variable (ψn) using aging or slip law
				# Step 5: compute slip-rate from friction law
				# TODO: Add slip law as well
				
				# Steps 4 and 5: at each point on fault
				for i in 1:length(fault_id)
					# Update the value of ψn in place (last variable factor is 1)
					ψn[i] =	state_time_evolution!(ψn[i], Vfn[i], dt, ics.Lc[i], par.Vo)
					
					# Calculate on-fault slip rate
					Vfn[i] = fault_slip_rate!(ψn[i], τn[i], ics.τo[i], ics.σo[i], 
										ics.a[i], ics.b[i], par.Vo, par.fo, 1)
				end
				
				#Vfnp1 .= Vfn	# used for nr-search during the dynamic part 

				# These equations are because we only have one side of the fault.
				Vfn .= 0.5*(Vfo .+ Vfn)
				glob.v[fault_id] .= 0.5*(Vfn .- par.Vpl)
				glob.v[creep_id] .= 0. 

			end #for p1 = 1:2

			# Steps 6,7,8,9 = repeat steps 2,3,4,5

			glob.v .= (glob.u .- u_pre)./dt 
			glob.v[fault_id] .= 0.5*(Vfn .- par.Vpl)
			glob.v[creep_id] .= 0. 

			glob.a .= 0.
			glob.u[creep_id] .= 0.
			glob.v[creep_id] .= 0.

		elseif solver == :dynamic

			# These steps are from Kaneko et al (2008). Also in Appendix B of Kaneko et al., 2011 
			# The step numbers below correspond to Appendix B.

			# Store previous displacement and velocity values
			u_pre .= glob.u
			v_pre .= glob.v 
			
			# Step 1: Update displacement and velocity 
			glob.u .= glob.u .+ dt*glob.v .+ (0.5*dt_min^2).*glob.a 
			
			# Step 1: Partial update of velocity 
			glob.v .= glob.v + (0.5*dt_min).*glob.a 
			glob.a .= 0

			# Step 2: Internal forces -K*d[t+1] stored in global array 'a'
			#=
			for el in 1:mesh.n_elements
				local_index = dof_id[:,:,el]	
				a_el = K_el[:,:,el]*glob.u[local_index][:]
				glob.a[local_index] .= glob.a[local_index] .- reshape(a_el, nnodes, nnodes) 
			end =#
			glob.a .= -Ksparse*glob.u 

			# Enforce a = 0 on creep boundary
			glob.a[creep_id] .= 0

			# Enforce absorbing boundaries
			glob.a[absorbing_id] .= glob.a[absorbing_id] .- (absorb_matrix .* glob.v[absorbing_id])

			#----------------------------------------------
			# Rate-state fault boundary
			#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

			# Step 3: Compute the stick traction (FaultVFree is stick traction)
			fault_vfree .= 2*glob.v[fault_id] .+ dt_min.*glob.a[fault_id] ./ Mglob[fault_id]

			# compute state variable using Vf from the previous time step
			Vfn .= 2*v_pre[fault_id] .+ par.Vpl 

			# Step 4: solve 2 nonlinear system of equations at each point on the fault
			#			while updating the state variable through a timestep
			#
			#			Eqn. 1: stick traction from step 4.
			#			Eqn. 2: Fault slip velocity from rate-state-friction
			# Dynamic fault boundary 
			for i in 1:length(fault_id)
				# Update state variable for one timestep
				ψn[i] = state_time_evolution!(ψn[i], Vfn[i], dt, ics.Lc[i], par.Vo)

				Vfnp1[i] = Vfn[i]	# for NR-search 2nd iteration initial guess

				# NR search: 1st iteration 
				Vfn[i], τn[i] =	nr_search!(τn[i], par.fo, par.Vo, ics.a[i], ics.b[i], 
									ics.σo[i], ics.τo[i], ψn[i], fault_z[i], fault_vfree[i])

				if isnan(Vfo[i]) || isnan(τn[i])
					println("Fault Location = ", i)
					println("Vf = ", Vfo[i])
					println("τn = ", τn[i])
					println("it = ", it)
					@error("NR SEARCH FAILED!")
					return
				end

				# 2nd iteration: only one iteration gives less accurate results (see Kaneko et al., 2008)
				ψn[i] = state_time_evolution!(ψn[i], 0.5*(Vfn[i] + Vfnp1[i]), dt, ics.Lc[i], par.Vo)
				Vfn[i], τn[i] =	nr_search!(τn[i], par.fo, par.Vo, ics.a[i], ics.b[i], 
									ics.σo[i], ics.τo[i], ψn[i], fault_z[i], fault_vfree[i])
					
			end
			
			# The rate-state equation above uses total stress = initial stress + on-fault stress 
			τn .= τn .- ics.τo

			# Step 5: update the internal forces with fault-boundary terms  
			glob.a[fault_id] .= glob.a[fault_id] .- fault_matrix.*τn
			
			#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
			# End of Fault Boundary
			#------------------------- --------------------

			# Step 6: Solve for acceleration
			glob.a .= glob.a./Mglob 

			# Ste 7: Complete the velocity update
			glob.v .= glob.v .+ (0.5*dt_min).*glob.a

			glob.v[creep_id] .= 0.
			glob.a[creep_id] .= 0.

		else
			@error("Check solver status")
		end

		# Current Sliprate 
		Vf_current = 2*glob.v[fault_id] .+ par.Vpl

		Vfmax = maximum(Vf_current)

		#-----------------------------------------
		# Output variables before and after events
		#------------------------------------------
		#Time-series variables to save 
		for (i,o) in enumerate(out_fault_id)
			slip = 2*glob.u[out_glob_id[i]] + par.Vpl*t
			sliprate = 2*glob.v[out_glob_id[i]] + par.Vpl
			shear_stress = (τn[o] + ics.τo[o])/1e6
			state = ψn[o]
			println(time_series[i], "$t $slip $sliprate $shear_stress $state")
		end
		
		println(time_series_1d, "$t $Vfmax \n")
		#------------------------------------------

		# Determine quasi-static or dynamic regime based on max-slip velocity
        if solver == :quasistatic && Vfmax < 5e-3 || solver == :dynamic && Vfmax < 2e-3
            solver = :quasistatic 
		else #solver == :quasistatic && Vfmax >= 5e-3 || solver == :dynamic && Vfmax >= 2e-3
            solver = :dynamic 
        end

		# Display some stuff to keep track
		if mod(it, 500) == 0
			println(" ")
			println("Solver = ", solver)
			println("Iteration = ", it)
			println("Vfmax = ", Vfmax)
			println("time = ", t/(365*24*60*60))
			println("dt = ", dt)
		end
		# Compute the next timestep 
		#return hcell, μmax, dt_min, par.dt_max, dt, fault_y, Vf_current, ics, solver
		dt = dt_evol(hcell, μmax, dt_min, par.dt_max, dt, fault_y, Vf_current, ics, solver)

	end # time loop 

	#println(tsave)
	# Close time-series output files 
	for i in 1:length(time_series)
		close(time_series[i])
	end
	close(time_series_1d)

	return tsave, Vfmaxsave
end # main function 