Julia code for simulating TCSPC measurements on proteins.

To run:

- copy `parameters.def.json` to a new file (it's just a set of defaults) and modify it as you need
- some parameters in there are given as lists of possible values; pick one value for each
- then either enter the Julia REPL and `include("src/Simulations.jl")` etc or just run `julia single_core.jl` in a terminal

This script will run your simulation, output histograms of all the decay processes happening in the system, plot them, and finally perform a reconvolution fit (by default it'll try to do monoexponential and biexponential fits).

The output will be placed in a directory labelled with some relevant parameters along with a copy of the parameter file.
