using Pkg; Pkg.activate("."); Pkg.resolve(); Pkg.instantiate(); Pkg.precompile()

using MPI, ArgParse, DelimitedFiles

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

MPI.Init()
const comm = MPI.COMM_WORLD
const root = 0

include("src/Simulations.jl"); using .Simulations
include("src/ReconvolutionFits.jl"); using .ReconvolutionFits

outpath = joinpath("out", "test_MPI")


proc_files, total_counts = Simulations.run(args["parameter_file"],
                             MPI.Comm_rank(comm),
                             MPI.Comm_size(comm),
                             outpath)
MPI.Barrier(comm)
hist_files = MPI.Reduce(proc_files, vcat, comm; root)

if MPI.Comm_rank(comm) == root
  """
  Loop over hist files, pick them up by run, read back in,
  sum, then fit the total files. Hopefully it should be clear
  from the names of the individual files etc what's happening.
  """
  # copy the paramter file to the output dir
  cp(args["parameter_file"],
     joinpath(dirname(hist_files[1]),
              args["parameter_file"]), force=true)
  sim, lattice, outdir = Simulations.setup(args["parameter_file"],
                                 MPI.Comm_size(comm), outpath)
  run_files = [[] for i = 1:sim.repeats]
  total_files = []
  for file in hist_files
    run = parse(Int, split(basename(file), "_")[2])
    push!(run_files[run], file)
    total_file = join(dirname(file), "hist_$(run)_total.dat")
    push!(total_files, total_file)
  end
  println("run files = $(run_files)")
  println("total files = $(total_files)")

  for (tf, rfs) in zip(total_files, run_files)
    (b, c, l, e) = Simulations.generate_histogram(sim, lattice)
    ec = zeros(Int, size(b))
    for rf in rfs
      raw = readdlm(rf)
      emissive = raw[2, 1:end]
      bincounts = raw[3:end, 1:end]
      ec += sum(bincounts[:, emissive .> 0], dims=2)
    end
    Simulations.write_hist(b, c, l, ec, tf)
  end
  println("total files: $(total_files)")

  for file in total_files
      irf_file = joinpath(dirname(file), "pulse.txt")

      τᵢ= [1.1e-9, 0.5e-9, 0.05e-9]
      for i = 1:length(τᵢ)
          taus = τᵢ[1:i]
          try
              ReconvolutionFits.fit(file, taus, "simulated", irf_file)
          catch e
              println("$(i)-exponential fit didn't work")
          end
      end
  end
end
