using Pkg; Pkg.activate("."); Pkg.resolve(); Pkg.instantiate(); Pkg.precompile()

using ArgParse

function parse_cmds()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--filename"
            help = "Data file"
            arg_type = String
            required = true
        "--irf_file"
            help = "IRF file"
            arg_type = String
            required = false
        "--type"
            help = "Data type"
            arg_type = String
            required = true
    end
    return parse_args(s)
end

args = parse_cmds()

include("src/ReconvolutionFits.jl"); using .ReconvolutionFits

τᵢ= [1.1e-9, 0.5e-9, 0.05e-9]

for i = 1:length(τᵢ)
    taus = τᵢ[1:i]
    try
        ReconvolutionFits.fit(args["filename"], taus, args["type"], args["irf_file"])
    catch e
        println("$(i)-exponential fit didn't work")
    end
end
