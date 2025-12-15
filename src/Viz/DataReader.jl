"""
    DataReader

Utilities for reading HDF5 simulation output data for visualization.
"""

using HDF5


"""
    read_timeseries_data(h5_file::String) -> (times, Vfmax, fault_depths, slip_rates)

Read time series data from HDF5 output file.

# Arguments
- `h5_file::String`: Path to HDF5 output file

# Returns
- `times::Vector`: Simulation times [s]
- `Vfmax::Vector`: Maximum fault slip rate at each timestep [m/s]
- `fault_depths::Vector`: Fault node depths [m]
- `slip_rates::Matrix`: Slip rate at each depth and time [m/s]

# Example
```julia
times, Vfmax, depths, slip_rates = read_timeseries_data("outputs/simulation.h5")
```
