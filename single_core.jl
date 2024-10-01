using Pkg; Pkg.activate("."); Pkg.instantiate(); Pkg.precompile()

using ArgParse

function parse_cmds()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--parameter_file"
            help = "Filename for parameters (JSON)"
            arg_type = String
            required = false
            default = "parameters.json"
    end
    return parse_args(s)
end

args = parse_cmds()

include("src/Simulations.jl"); using .Simulations
include("src/ReconvolutionFits.jl"); using .ReconvolutionFits

outpath = joinpath("out", "test")
hist_files = Simulations.run(args["parameter_file"], 0, outpath)
# copy the paramter file to the output dir
cp(args["parameter_file"],
   joinpath(dirname(hist_files[1]),
            args["parameter_file"]), force=true)

for file in hist_files
    irf_file = joinpath(dirname(file), "pulse.txt")
    τᵢ= [1.1e-9]
    try
      ReconvolutionFits.fit(file, τᵢ, irf_file)
    catch e
      println("Monoexponential fit didn't work")
    end
    τᵢ= [1.1e-9, 0.1e-9]
    try
      ReconvolutionFits.fit(file, τᵢ, irf_file)
    catch e
      println("Monoexponential fit didn't work")
    end

end
