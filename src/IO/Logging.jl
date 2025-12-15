"""
    Logging

Structured logging utilities for earthquake cycle simulations.

Provides logging levels and formatted output for simulation events.
"""

using Logging
using Printf
using Dates


"""
    setup_logging(config::SimulationConfig)

Configure logging for simulation.

# Arguments
- `config::SimulationConfig`: Simulation configuration

# Sets up:
- Log level based on configuration (Info, Debug, Warn)
- Log file in output directory (optional)
- Console logging with timestamps

# Example
```julia
setup_logging(config)
@info "Simulation started"
```
"""
function setup_logging(config::SimulationConfig)
    # Default to Info level
    log_level = Logging.Info

    # Create log file in output directory
    log_file = joinpath(config.simulation.output_dir,
                       "$(config.simulation.name).log")

    # Set up logging
    # For now, use default console logging
    # TODO: Add file logging if needed

    @info "Logging configured" level=log_level

    return nothing
end


"""
    log_event(event_type::Symbol, state::SimulationState, details::Dict)

Log a simulation event with structured data.

# Arguments
- `event_type::Symbol`: Event type (`:earthquake`, `:mode_switch`, `:checkpoint`, etc.)
- `state::SimulationState`: Current simulation state
- `details::Dict`: Additional event-specific details

# Example
```julia
log_event(:earthquake, state, Dict(
    "magnitude" => 6.5,
    "nucleation_depth" => 12.5
))
```
"""
function log_event(event_type::Symbol, state::SimulationState, details::Dict=Dict())
    timestamp = now()
    iter = state.iteration
    time = state.time

    if event_type == :earthquake
        @info "Earthquake event" iter=iter time=time details...

    elseif event_type == :mode_switch
        @info "Solver mode switch" iter=iter time=time details...

    elseif event_type == :checkpoint
        @info "Checkpoint saved" iter=iter time=time details...

    else
        @info "Simulation event" type=event_type iter=iter time=time details...
    end

    return nothing
end
