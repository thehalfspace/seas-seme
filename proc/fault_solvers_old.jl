
# These don't work for some reason. Have Peng try the nr-search with sem2dpack using state variable θ

# On-fault slip-rate: Equation 14 from Kaneko et al. (2011)
function fault_slip_rate(τn, θnp1, par, ics)
	"""
	Equation 14 Kaneko et al. (2011)
	TODO: Add slip law here as well
	"""
	#2 .* par.Vo .* sinh.(τn ./ (ics.a .* ics.σo)) .* 
	#exp.((-par.fo .- ics.b .* log.(abs.(par.Vo .* θnp1 ./ ics.Lc))) ./ ics.a)
	psi = log.(par.Vo .* θnp1 ./ ics.Lc)
	term1 = -(par.fo .+ ics.b .*psi)./ics.a 
	term2 = τn ./(ics.σo .* ics.a)
	
	return par.Vo .* (exp.(term1 .+ term2) .- exp.(term1 .- term2))
end

# State variable timestep evolution 
function state_time_evolution!(θnp1, θn, Vn, dt, par, ics, solver)
	"""
		Equation 20 Kaneko et al. (2011)
		TODO: Add slip law here as well
	"""
	if solver == :quasistatic
		#temp = ics.Lc ./ abs.(Vf)
		#return (θn .- temp).* exp.(-dt ./ temp) .+ temp 

		# This is as is from Kaneko 2011: Equation 20
		#return θn .* exp.(-(abs.(Vn) ./ ics.Lc) .* dt) .+ 
		#		(ics.Lc ./ abs.(Vn)) .* (1 .- exp.(-(abs.(Vn) ./ ics.Lc) .*dt))


		# = This is from Kaneko's matlab code: 
		# I think they use taylor expansion for exponential term in the above equation
		# ψ = log(θV/Lc)
		#--------------
		for i in 1:length(Vn)
			VdtL = abs(Vn[i]*dt / ics.Lc[i])
			if VdtL < 1.0e-6
				θnp1[i] = θn[i]*exp(-VdtL) + (ics.Lc[i]/abs(Vn[i])) *
							(par.Vo*dt/ics.Lc[i] - 0.5*par.Vo*abs(Vn[i])*dt^2/(ics.Lc[i]^2))
			else
				θnp1[i] = θn[i]*exp(-VdtL) + (par.Vo*ics.Lc[i]/Vn[i]^2)*(1 - exp(-VdtL))
			end
		end
		# =#

	elseif solver == :dynamic
		#temp = ics.Lc ./ abs.(Vf)
		#return (θn .- temp).* exp.(-dt ./ temp) .+ temp 
		
		#return θn .* exp.(-(abs.(Vf) ./ ics.Lc) .* dt) .+ 
		#		(ics.Lc ./ abs.(Vf)) .* (1 .- exp.(-(abs.(Vf) ./ ics.Lc) .*dt))

		# This is from Kaneko's matlab code: 
		# I think they use taylor expansion for exponential term in the above equation
		# ψ = log(θV/Lc)
		#--------------
		for i in 1:length(Vn)
			VdtL = abs(Vn[i]*dt / ics.Lc[i])
			if VdtL < 1.0e-6
				θnp1[i] = θn[i]*exp(-VdtL) + (ics.Lc[i]/abs(Vn[i])) *
							(par.Vo*dt/ics.Lc[i] - 0.5*par.Vo*abs(Vn[i])*dt^2/(ics.Lc[i]^2))
			else
				θnp1[i] = θn[i]*exp(-VdtL) + (par.Vo*ics.Lc[i]/Vn[i]^2)*(1 - exp(-VdtL))
			end
		end
	else
		@error("Check solver tag")
	end
	θnp1 
end

# Non-linear friction solver on the fault : functions to pass into NLsolve.jl 
# 	Equations from Kaneko et al. (2008): Eqns 11 and 15 (rearranged)
#	x[1] = Vfn = 2Vo sinh(τn/aσ) exp((-fo - b.log(Voθ/Lc)))
#	x[2] = τn = τo + Z*Vfree - Z*Vfn 
function friction_equations!(F, x, fault_z, fault_vfree, θnp1, Vo, a, b, σo, τo, fo, Lc)
	F[1] = 2*Vo * sinh(x[2]/(a * σo)) * exp((-fo - b*log(abs(Vo*θnp1/Lc))) /a) - x[1]
	F[2] = fault_z*fault_vfree - fault_z*x[1] + τo - x[2]
end

# The jacobian derivatives for the above equations to pass into nlsolve
# Cross check with sem2dpack/SRC/bc_dynflt_rsf.f90  
function friction_jacobians(J, x, fault_z)
	return 0
end


############# NLSOLVE for fault solver stuff #############################
			# Step 4: Determine first prediction of state variable (θnp1) using aging or slip law
			#=
			state_time_evolution!(θnp1,	θn, Vfn, dt, par, ics, solver) 




			# Appendix B: Step 4: Find fault traction and slip-velocity statisfying friction law
			# Nonlinear solver: NR Search algorithm
			for i in 1:length(fault_id)
				initial_guess = [0.;0.]
				tol = 0.001*ics.a[i]*ics.σo[i]
				#initial_guess = [Vfn[i];τn[i]]

			 	temp = nlsolve( (F,x) -> friction_equations!(F, x, fault_z[i], fault_vfree[i], 
			 								θnp1[i], par.Vo, ics.a[i], ics.b[i], ics.σo[i], 
			 								ics.τo[i], par.fo, ics.Lc[i]), 
			 								initial_guess, ftol = tol)

			 	Vfnp1[i] = temp.zero[1]
			 	τnp1[i]  = temp.zero[2]
			end

			θn .= θnp1

			state_time_evolution!(θnp1,	θn, 0.5*(Vfnp1 .+ Vfn), dt, par, ics, solver) 

			# Nonlinear search: 2nd iteration 
			for i in 1:length(fault_id)
				initial_guess = [0.;0.]
				tol = 0.001*ics.a[i]*ics.σo[i]
				#initial_guess = [Vfnp1[i];τnp1[i]]
				temp = nlsolve( (F,x) -> friction_equations!(F, x, fault_z[i], fault_vfree[i], 
											θnp1[i], par.Vo, ics.a[i], ics.b[i], ics.σo[i], 
											ics.τo[i], par.fo, ics.Lc[i]), 
											initial_guess, ftol = tol)
				
				Vfnp1[i] = temp.zero[1]
				τnp1[i]  = temp.zero[2]
			end	

			τn .= τnp1 .- ics.τo
			θn .= θnp1 
			Vfn .= Vfnp1
			=#