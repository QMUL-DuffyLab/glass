using Distributed
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

@everywhere begin
    using Pkg; Pkg.activate(".")
    Pkg.instantiate(); Pkg.precompile()
end

# NB: gonna need more setup stuff here:
# - protein setup since some rates will change
# - pulse and simulation parameters (via JSON?)
@everywhere begin
    include("src/Simulations.jl"); using .Simulations
    proteins = [make_lh2(args["t_t_rate"])]
    rho = [1.0]
    lattice = lattice_generation(get_lattice("hex"), nmax,
            proteins, rho)

    pulse_params = PulseParams(200e-12, 50e-12, 800e-9, args["fluence"])
    sim = SimulationParams(15e-9, 1e-12, args["rep_rate"] * 1e6,
                           1e-9, 25e-12, args["pulse_params"],
                           args["n_counts"], args["n_repeats"])
    protein_names = [p.name for p in proteins]
    outdir = joinpath("out",
        join([protein_names..., sim.rep_rate/1e6, "MHz_fluence_",
              pulse_params.f, "_T_T_rate_", t_t_rate, "_"]))
    mkpath(outdir)
    lattice_plot_file = joinpath(outdir, "lattice.svg")
    pulse_file = joinpath(outdir, "pulse.txt")
    hist_path = joinpath(outdir, "hist")

    pulse = construct_pulse(pulse_params, sim.dt1)
    x = (0:ceil((2.0 * pulse_params.μ)/sim.dt1)) * sim.dt1

    open(pulse_file, "w") do io
        writedlm(io, hcat(x, pulse))
    end

    plot_lattice(lattice, lattice_plot_file)

    (bins, counts, labels, ec) = generate_histogram(sim, lattice)

    run = 1
    while run <= sim.repeats
        hist_file = "$(hist_path)_$(run).txt"
        # need to reduce over this
        counts += one_run(sim, lattice, seed_start + run, hist_file)
        run += 1
    end

end
