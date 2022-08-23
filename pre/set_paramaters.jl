#######################################################################
#######################################################################

#	PARAMETER FILE: SET THE PHYSICAL PARAMETERS FOR THE SIMULATION

#######################################################################
#######################################################################

struct parameters
	"""
		-----------
		PARAMETERS
		-----------
	"""

	yr2sec::Int
	Total_Time::Int
	
	dt_max::Float64		# Maximum Timestep 
	Vpl::Float64		# Plate Loading Rate
	fo::Float64			# Reference Friction Coefficient
	Vo::Float64			# Reference Slip Rate
	Vthres::Float64		# Velocity Threshold to Capture Earthquakes
	CFL::Float64		# CFL Stability Criteria
end

function set_default_parameters()
	"""
		Unpack and initialize the parameters here
	"""

	yr2sec = 365*24*60*60
	Total_Time = 50*yr2sec
	dt_max = 100 * 24 * 60 * 60 # days

	Vpl = 35e-3/yr2sec	# Plate Loading Rate
	fo = 0.6			# Reference Friction Coefficient
	Vo = 1.0e-6			# Reference Slip Rate

	Vthres = 1.0e-3		# Velocity Threshold to Capture Earthquakes
	CFL = 0.6			# CFL Stability Criteria

	return parameters(yr2sec, Total_Time, dt_max, Vpl, fo, Vo, Vthres, CFL)
end