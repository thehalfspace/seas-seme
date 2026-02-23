"""
    Logging

Structured logging utilities for earthquake cycle simulations.

Provides logging levels and formatted output for simulation events.
"""

using Logging
using Printf
using Dates
using LoggingExtras


"""
    setup_logging(config::SimulationConfig; console_level=Logging.Info, file_level=Logging.Info)

Configure logging for simulation to both console and file.

# Arguments
- `config::SimulationConfig`: Simulation configuration
- `console_level`: Minimum log level for console output (default: Info)
- `file_level`: Minimum log level for file output (default: Info)

# Sets up:
- Dual logging to console and file
- Log file: `data/{simulation_name}/output.log`
- Timestamped log entries
- Formatted output

# Example
```julia
setup_logging(config)
@info "Simulation started"
```

# Returns
- `log_file::String`: Path to the log file
"""
function setup_logging(config::SimulationConfig;
                      console_level=Logging.Info,
                      file_level=Logging.Info)
    # Create log file in simulation directory
    log_file = joinpath(config.simulation.output_dir, "output.log")

    # Ensure directory exists
    mkpath(dirname(log_file))

    # Open log file for writing
    log_io = open(log_file, "w")

    # Create timestamp logger for file
    timestamp_logger_file = TransformerLogger(SimpleLogger(log_io, file_level)) do log
        merge(log, (; message = "$(Dates.format(now(), "yyyy-mm-dd HH:MM:SS")) | $(log.message)"))
    end

    # Create timestamp logger for console
    timestamp_logger_console = TransformerLogger(ConsoleLogger(stderr, console_level)) do log
        merge(log, (; message = "$(Dates.format(now(), "HH:MM:SS")) | $(log.message)"))
    end

    # Create a TeeLogger to write to both console and file
    logger = TeeLogger(
        MinLevelLogger(timestamp_logger_console, console_level),
        MinLevelLogger(timestamp_logger_file, file_level)
    )

    # Set as global logger
    global_logger(logger)

    @info "Logging configured" log_file=log_file console_level=console_level file_level=file_level

    return log_file, log_io
end


"""
    close_logging(log_io::IO)

Close the log file stream.

# Arguments
- `log_io::IO`: Log file IO stream

# Example
```julia
_, log_io = setup_logging(config)
# ... simulation runs ...
close_logging(log_io)
```
"""
function close_logging(log_io::IO)
    flush(log_io)
    close(log_io)
    @info "Log file closed"
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
function log_event(event_type::Symbol, state::AbstractSimulationState, details::Dict=Dict())
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
