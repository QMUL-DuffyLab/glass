using Distributed

n_procs = 2
addprocs(n_procs)

@everywhere begin
    using Pkg; Pkg.activate(".")
    Pkg.instantiate(); Pkg.precompile()
end

# NB: gonna need more setup stuff here:
# - protein setup since some rates will change
# - pulse and simulation parameters (via JSON?)
@everywhere begin
    include("src/Simulations.jl"); using .Simulations
end
