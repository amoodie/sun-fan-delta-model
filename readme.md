# Implementation of adapted Sun fan delta model

## Running the model.

Model was developed and run with Matlab R2021a.

In summary, the model is run using wrappers. The wrapper defines
parameters for the model run(s), and then executes simulations in parallel.
The main simulations are executed via `wrappers/wrapper_simulation_set.m`.


## Analysis of outputs
The model outputs a series of .mat files that contain the full output of
the state fields of the model (e.g., `eta`, `Qw`). We process these files
into a netCDF file with coordinates for use with [DeltaMetrics](https://github.com/DeltaRCM/DeltaMetrics).
