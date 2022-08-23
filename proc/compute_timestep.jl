function dt_evol(hcell, μmax, dt_min, dt_max, dt, fault_y, Vfn, ics, solver)
	"""
		Algorithm from Lapusta et al. (2000)
		Check Kaneko et al., 2011 Appendix C
	"""
	if solver == :quasistatic 
		ξthf = 1
		ξmax = 0.5
		ξLf = 0.
		ξth = 0.

		dt_incf = 1.2

		dtnx = dt_max # initial value of dt

		for i in 1:length(fault_y)
			
			# Compute time restricting parameters
			expr1 = -(ics.a[i] - ics.b[i])/ics.a[i]
			expr2 = 0.25π*(μmax/hcell) * (ics.Lc[i]/(ics.a[i]*ics.σo[i]))
			ro = expr2 - expr1

			if (0.25*ro*ro - expr2) >= 0
				ξth = 1/ro
			else
				ξth = 1 - expr1/expr2
			end

			if ξthf*ξth > ξmax 
				ξLf = ξmax*ics.Lc[i]
			else
				ξLf = ξthf*ξth*ics.Lc[i]
			end

			#println(ξLf)
			if abs(Vfn[i])*dt_max > ξLf 
				dt_cell = ξLf/abs(Vfn[i])
				if dt_cell < dtnx
					dtnx = dt_cell 
				end
			end
		end
		
		if dt_min > dtnx
            dtnx = dt_min
        end

        if dtnx > dt_incf*dt
            dtnx = dt_incf*dt
        end

        dt = dtnx

	elseif solver == :dynamic 
		dt = dt_min
	else 
		@error("Check solver tag")
	end
	dt
end