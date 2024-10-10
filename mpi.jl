using MPI
using ArgParse

function parse_cmds()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--n_procs"
            help = "Number of cores to use"
            arg_type = Int
            required = false
            default = 4
    end
    return parse_args(s)
end

args = parse_cmds()
addprocs(args["n_procs"])

"""
mpi broadcast the name of the parameter file?
then get each proc to print out the simulation params
etc and check it's in the right form?
then check the sum over histogram in the code and
figure out how to do that in the reduce call
"""
