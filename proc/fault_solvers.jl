# Solvers on the fault
# All these are directly from Kaneko et al's matlab code. Try to solve this without the transformation 
# from θ to ψ.

# Reduce the Vf variables in original code and see what works, then copy here 
function fault_slip_rate!(ψn, τn, τo, σo, a, b, Vo, fo, fact)
	"""
		Quasi-static fault slip-rate calculator
		Returns Vn at fault location
	"""
	#state_time_evolution_2!(ψn, Vn, dt, Lc, Vo) # do this in the main loop 

	term1 = fact*(τn + τo)/(σo*a)
	term2 = -(fo + b*ψn)/a
	return Vo *(exp(term2 + term1) - exp(term2 - term1))
end

# Transformed state variable time evolution 
function state_time_evolution!(ψn, Vn, dt, Lc, Vo)
	"""
		Assuming IDState = 2 in Kaneko's code
		Scalar function: point-by-point calculation on the fault 
		
		This time-integration is a different scheme not mentioned in Kaneko's paper. 
	"""
	Vdt_L = abs(Vn)*dt/Lc # Vn.dt/Lc
	if Vdt_L < Vo
		return log( exp(ψn-Vdt_L) + Vo*dt/Lc - 0.5*Vo*abs(Vn)*dt*dt/(Lc^2) )
	else
		return log( exp(ψn-Vdt_L) + (Vo/abs(Vn))*(1 - exp(-Vdt_L)) )
	end
end

# Newton-Raphson search for dynamic fault boundary
function nr_search!(τn, fo, Vo, a, b, σo, τo, ψn, fault_z, fault_vfree)
	Vw = 1.0e10
	fact = 1 + (Vo/Vw)*exp(-ψn)

	#fa::BigFloat = 0.
    #help1::BigFloat = 0.
    #help2::BigFloat = 0.
    #delta::BigFloat = 0.
	#Vnprime::BigFloat = 0.

	tol = 0.001*a*σo	# tolerance (from Kaneko's Matlab code)
	k = 0
	delta = Inf 

	while abs(delta) > tol
		fa = fact*τn/(σo*a)
        help = -(fo + b*ψn)/a

        help1 = exp(help + fa)
        help2 = exp(help - fa)

        Vn = Vo*(help1 - help2)
		#Vn = fault_slip_rate_qs!(ψn, τn, 0., σo, Vo, fact)

		Vnprime = fact*(Vo/(a*σo))*(help1 + help2)
		delta = (fault_z*fault_vfree - fault_z*Vn + τo - τn)/(1 + fault_z*Vnprime)
		τn = τn + delta 
		k = k + 1
		if abs(delta) > 1e10 || k == 1000
            println("k = ", k)
            @error("NR search fails to converge")

            return Float64(Vn), Float64(τn)
        end
	end 
	#Vn = Float64(Vn)
	#Vn = fault_slip_rate!(ψn, τn, τo, σo, a, b, Vo, fo, fact)

	# Calculate the sliprate again
	fa = fact*τn/(σo*a)
    help =  -(fo + b*ψn)/a
    help1 = exp(help + fa)
    help2 = exp(help - fa)
    Vn = Vo*(help1 - help2)

	return Vn, τn
end
