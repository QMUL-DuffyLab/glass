using Distributed

# still testing stuff - hardcoding this for now
addprocs(4)

@everywhere begin
  using Pkg; Pkg.activate(@__DIR__)
  Pkg.instantiate(); Pkg.precompile()
end
# this has to be in a separate everywhere block, for some reason
# if you put this line in the above begin/end block it doesn't work.
# i have no idea why this should be. julia syntax! :)
@everywhere using ArgParse

@everywhere begin
  function parse_cmds()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--parameter_file"
            help = "Filename for parameters (JSON)"
            arg_type = String
            required = false
            default = "parameters.json"
        "--outpath"
            help = "Filename for parameters (JSON)"
            arg_type = String
            required = false
            default = "out"
    end
    return parse_args(s)
  end
end

# same as with the ArgParse include above :)))))
@everywhere include("src/Simulations.jl")
@everywhere using .Simulations

@everywhere begin
    args = parse_cmds()
    outpath = joinpath(pwd(), args["outpath"])
    mkpath(outpath)
    hist_files = Simulations.run(args["parameter_file"], myid(),
    # 5 here because the main process will run plus n_procs workers
                                 5, outpath)
end
