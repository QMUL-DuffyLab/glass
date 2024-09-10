using Pkg; Pkg.activate("."); Pkg.instantiate(); Pkg.precompile()

include("src/Simulations.jl"); using .Simulations
include("src/ReconvolutionFits.jl"); using .ReconvolutionFits

hist_files = Simulations.run("parameters.json", 0)

for file in hist_files
    irf_file = joinpath(dirname(file), "pulse.txt")
    τᵢ= [1.1, 0.1]
    ReconvolutionFits.fit(file, τᵢ, irf_file)
end
