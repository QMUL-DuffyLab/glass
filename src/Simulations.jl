module Simulations
include("Lattices.jl")

using JSON, Random, DelimitedFiles, PyPlot

struct PulseParams
    μ::Float64
    fwhm::Float64
    # unsure if we'll need this or not
    λ::Float64
    f::Float64
end

struct SimulationParams
    tmax::Float64
    dt1::Float64
    rep_rate::Float64
    dt2::Float64
    binwidth::Float64
    pulse_params::PulseParams
    maxcount::Int
    repeats::Int
end

trapz(A) = sum(A) + (A[begin] + A[end]) / 2.0
trapz(A, dx) = dx * trapz(A)

function construct_pulse(p, dt)
    tmax = 2.0 * p.μ
    σ = p.fwhm / (2.0 * sqrt(2.0 * log(2.0)))
    pulse = [(1.0 / (σ * sqrt(2.0 * pi))) * 
             exp(-((i * dt) - p.μ)^2 / (sqrt(2.0) * σ)^2)
             for i=0:ceil(tmax / dt)]
    pulse .*= (p.f / trapz(pulse, dt))
end

"""
use some pretty simple criteria to determine whether a given type of move
is currently possible on site i.
the annihilation criteria's less trivial, not fixed it yet
"""
function possible_moves(l::Lattice, i::Int, n, ft::Float64)
    move_types = []
    # pulse has to be on and there has to be a state with a nonzero cross section
    nₛ= l.proteins[l.sites[i].ip].nₛ
    xsecs = l.proteins[l.sites[i].ip].xsec
    if ft > zero(ft) && any(xsecs .> 0.0)
        append!(move_types, ["gen"])
    end
    # product of xsec and current population must be > 0 in order for a second
    # photon to be absorbed and cause stimulated emission
    # println("xsecs = ", xsecs)
    # println("i = ", " protein = ",
    #         l.proteins[l.sites[i].ip].name, " n = ", n[i, 1:nₛ])
    if ft > zero(ft) && any(xsecs .* n[i, 1:nₛ] .> 0.0)
        append!(move_types, ["se"])
    end
    if any(n[i, 1:nₛ] .> 0)
        append!(move_types, ["hop", "decay"])
        if nₛ > 1
            # can be cleverer about this too - there must be a state
            # with available population and a nonzero transfer rate
            append!(move_types, ["transfer"])
        end
    end
    # this isn't always gonna work. need to think about it
    if sum(n[i, 1:nₛ]) > 1
        append!(move_types, ["ann"])
    end
    return move_types
end

"""
Based on the current state of the pulse and current populations,
figure out which moves are possible and whether 
excitation energy will be lost as a result of the move. Return the rate of
the process and information about it, so we can then run it through
Metropolis and bin any decays if the move's accepted.

We don't need to figure out whether it's emissive here, so long as we 
figure out at the start of the simulation which decays are emissive, 
then we can just index back into that.

- `n`: matrix of occupations of each state on each site.
- `i`: current site index.
- `l`: the lattice.
- `ft`: the current intensity of the pulse (only considering excitation
  at one wavelength, at least currently.)

NB: this can be sped up by only picking states with nonzero populations.
Figure out how to do that.
"""
function propose_move(l::Lattice, i::Int, n, ft::Float64)
    
    rate = 0.0
    ip = l.sites[i].ip
    p = l.proteins[ip]
    if ip > 1
        loss_start_index = sum(p.nₛ * (p.nₛ + 2)
                               for p in l.proteins[1:(ip - 1)])
    else
        loss_start_index = 0
    end
    loss_index = 0
    isi = i
    fsi = i
    ist = 0
    fst = 0

    # get possible move types. this won't be empty because we check
    # whether there are any in the monte carlo step before this function's called
    move_types = possible_moves(l, i, n, ft)

    # pick a type of move
    mt = move_types[rand(1:length(move_types))]
    # println(mt)
    if mt == "gen"
        s = rand(findall(>(0), p.xsec))
        rate = (ft * p.xsec[s] *
                (p.n_tot[p.ps[s]] - n[i, s]) / p.n_tot[p.ps[s]])
        # rate = ft * p.xsec[s] * (p.n_tot[p.ps[s]] - n[i, s])
        ist = s
        fst = s
    elseif mt == "se"
        s = rand(findall(>(0), p.xsec))
        rate = ft * p.xsec[s] * n[i, s] / p.n_tot[p.ps[s]]
        loss_index = loss_start_index + s
        ist = s
        fst = s
    elseif mt == "hop"
        s = rand(findall(>(0), n[i, 1:p.nₛ]))
        nn = rand(l.nn[:, i])
        while nn == 0
            nn = rand(l.nn[:, i])
        end
        rate = n[i, s] * p.hop[s]
        # entropy factor - rate is reduced if hopping to a higher-occupied state
        if n[i, s] < n[nn, s]
            rate *= (n[i, s] * (p.n_thermal[p.ps[s]] - n[nn, s])) / 
            ((n[nn, s] + 1) * (p.n_thermal[p.ps[s]] - n[i, s] + 1))
        end
        ist = s
        fsi = nn
        fst = s
    elseif mt == "transfer"
        # 
        si = rand(findall(>(0), n[i, 1:p.nₛ]))
        sf = rand(1:p.nₛ)
        while sf == si
            # this would be a decay process - reroll
            sf = rand(1:p.nₛ)
        end
        rate = n[i, si] * p.intra[si, sf] # is this right? check
        ist = si
        fst = sf
    elseif mt == "decay"
        s = rand(findall(>(0), n[i, 1:p.nₛ]))
        rate = n[i, s] * p.intra[s, s]
        # println("DECAY: ", n[i, s], " ", p.intra[s,s],
        #         " typeof(s) = ", typeof(s),
        #         " typeof(",p.intra,") = ", typeof(p.intra),
        #         " typeof(",p.intra[s, :],") = ", typeof(p.intra[s, :]),
        #         " typeof(",p.intra[:, s],") = ", typeof(p.intra[:, s]),
        #         " typeof(",p.intra[s, s],") = ", typeof(p.intra[s, s]),
        #        )
        emissive = p.emissive[s]
        loss_index = loss_start_index + p.nₛ + s
        ist = s
        fst = s
    elseif mt == "ann"
        """
        this will probably work but be inefficient. the cleverer way
        would be to pick one randomly, then select a second species
        for which the corresponding annihilation rate > 0
        """
        # si = rand(findall(>(0), n[i, 1:p.nₛ]))
        # sf = rand(findall(>(0), n[i, 1:p.nₛ]))
        # sf = rand(1:p.nₛ)
        # ordering ensures annihilation counting is commutative
        pair = sort!(rand(findall(>(0), n[i, 1:p.nₛ]), 2))
        # println("p.dist = ", p.dist, " type = ", typeof(p.dist),
        #         " p.ann = ", p.ann, " type = ", typeof(p.ann),
        #         " pair = ", pair, " p.dist[pair] = ", p.dist[pair],
        #         " p.dist[pair expl] = ", p.dist[pair[1], pair[2]])
        if p.dist[pair...]
            # distinguishable
            rate = n[i, pair[1]] * n[i, pair[2]] * p.ann[pair...]
        else
            neff = n[i, pair[1]] + n[i, pair[2]]
            rate = p.ann[pair...] * (neff * (neff - 1)) / 2.0
        end
        # one column for every combination
        loss_index = loss_start_index + (2 + (pair[1] - 1)) * p.nₛ + pair[2]
        ist = pair[1]
        fst = pair[2]
    else
        println("shouldn't be here! mt = ", mt)
    end
    return (rate, mt, loss_index, isi, ist, fsi, fst)
end

""" 
carry out an accepted move.
- `n`: matrix of occupations of states on each site
- `move_type`: string to tell us what the effect of the move is
- `isi`: index of initial site
- `ist`: index of initial state on `sites[isi]`
- `fsi`: index of final site
- `fst`: index of final state on `sites[fsi]`
- `which_ann`: integer to tell us which of the states loses population.
  Check which site this corresponds to and decrement that one.
"""
function move!(n, move_type, isi, ist,
    fsi, fst, which_ann)
    if move_type == "hop"
        n[isi, ist] -= 1
        n[fsi, fst] += 1
    elseif move_type == "transfer"
        n[isi, ist] -= 1
        n[fsi, fst] += 1
    elseif move_type == "decay"
        n[isi, ist] -= 1
    elseif move_type == "gen"
        n[isi, ist] += 1
    elseif move_type == "ann"
        if which_ann == ist
            n[isi, ist] -= 1
        else
            n[fsi, fst] -= 1
        end
    end
    if n[isi, ist] < 0
        println("n[$(isi), $(ist)] = $(n[isi, ist])")
        println("$(move_type), [fsi, fst] = $(fsi), $(fst)")
        throw(DomainError(n[isi, ist],
              "populations must be non-negative"))
    elseif n[fsi, fst] < 0
        println("n[$(fsi), $(fst)] = $(n[fsi, fst])")
        println("$(move_type), [isi, ist] = $(isi), $(ist)")
        throw(DomainError(n[fsi, fst],
              "populations must be non-negative"))
    end
end

"""
one monte carlo step. 
attempt one of each type of possible move on each site on average;
`propose_move` picks a possible move.

once a move's been picked, run it through Metropolis; if it's accepted,
carry it out and increment the histogram iff 1.) an excitation's been lost
and 2.) we are currently binning decays.
NB: currently i pass pulse_params to check time but it'd be easier
to just send a 0 or an array of zeroes so we don't have to do that check
"""
function mc_step!(l::Lattice, n, pulse, pulse_params,
    t::Real, dt::Real, bin::Bool, counts::Matrix{Integer},
    binwidth::Real, print_debug::Bool)
    
    # current intensity of pulse
    ft = zero(eltype(pulse))
    if t <= 2.0 * pulse_params.μ
        pind = round(Int, t / dt) + 1
        if pind <= length(pulse)
            ft = pulse[pind]
        end
    end
    
    # pick each site once on average
    for i=1:length(l.sites)
        s = rand(1:length(l.sites))
        move_types = possible_moves(l, s, n, ft)
        if isempty(move_types)
            continue
        end
        # attempt each possible move type once on average
        for j=1:length(move_types)
            # check whether any moves are possible now
            move_types = possible_moves(l, s, n, ft)
            if isempty(move_types)
                break
            end
            
            (rate, move_type, loss_index, isi, ist, fsi, fst) = propose_move(l, s, n, ft)
        
            # metropolis
            prob = rate * dt * exp(-1.0 * rate * dt)
            if rand() < prob
                # carry out the move - this changes population, which
                # is why we have to check possible moves again above
                if print_debug
                    if move_type == "gen"
                        println(
                                "t = $(t), dt = $(dt), bin = $(bin), ",
                                "ft = $(ft), rate = $(rate) ",
                                "prob = $(prob), ",
                                "n[isi, ist] = $(n[isi, ist])",
                                " n[fsi, fst] = $(n[fsi, fst])",
                               )
                    end
                end
                p = l.proteins[l.sites[s].ip]
                move!(n, move_type, isi, ist, fsi, fst,
                      p.which_ann[ist, fst])
                if bin && loss_index > 0
                    # increment the correct column of the histogram
                    counts[loss_index, floor(Int, t / binwidth) + 1] += 1
                end
            end
        end
    end
end


"""
create the histogram array we'll bin counts into, the bin values themselves,
and the labels for each column of the count array.
this won't work with `propose_move` yet: we need to add a counter for which
protein we're on.
"""
function generate_histogram(s::SimulationParams, l::Lattice)
    n = ceil(Int, (s.tmax / s.binwidth) + 1)
    bins = [(i - 1) * s.binwidth for i=1:n]
    labels = []
    emissive_columns = []
    n_losses = 0
    for p in l.proteins
        for s1 in p.states
            n_losses += 1
            append!(labels, [p.name * "_se_" * s1])
        end
        for (i, s1) in enumerate(p.states)
            n_losses += 1
            append!(labels, [p.name * "_decay_" * s1])
            if p.emissive[i]
                append!(emissive_columns, [n_losses])
            end
        end
        for s1 in p.states
            for s2 in p.states
                n_losses += 1
                append!(labels, [p.name * "_ann_" * s1 * "_" * s2])
            end
        end
    end
    counts = zeros(Integer, n_losses, n)
    return (bins, counts, labels, emissive_columns)
end

function write_hist(bins, counts, labels, ec, file)
    open(file, "w") do io
        # first column is always bins
        # first row labels, second row emissive
        write(io, join(["bins", labels...], '\t') * '\n')
        emissive = [false for l in labels]
        emissive[ec...] = true
        write(io, join([false, emissive...], '\t') * '\n')
        writedlm(io, hcat(bins, transpose(counts)))
    end
end

function plot_counts(bins, counts, labels, maxcount, outfile)
    # only plot columns with nonzero counts
    cols = [sum(c).>0 for c in eachcol(counts)]
    cplot = counts[:, cols]
    lplot = labels[cols]
    fig, ax = plt.subplots(figsize=(12,8))
    for i = 1:sum(cols)
        plot(bins, cplot[:, i], label=permutedims(lplot)[i])
    end
    ax.set_ylim(1, 2.0 * maxcount)
    ax.set_yscale("log")
    plt.grid(visible=true)
    ax.legend()
    ax.set_xlabel("time (s)")
    ax.set_ylabel("counts")
    savefig(outfile)
end

function one_run(sim, lattice, seed, hist_file)

    Random.seed!(seed)
    (bins, counts, labels, ec) = generate_histogram(sim, lattice)
    n = zeros(Int, length(lattice.sites),
              maximum([p.nₛ for p in lattice.proteins]))
    pulse = construct_pulse(sim.pulse_params, sim.dt1)

    # time from the end of one data collection period to the next pulse
    pulse_interval = (1.0 / sim.rep_rate)

    curr_maxcount = 0
    rep = 1
    print_debug = false
    while curr_maxcount < sim.maxcount
        t = 0.0
        while t < sim.tmax
            mc_step!(lattice, n, pulse, sim.pulse_params, t, sim.dt1,
                     true, counts, sim.binwidth, print_debug)
            t += sim.dt1
            if t >= 2.0 * sim.pulse_params.μ && sum(n) == 0
                # no excitations and no pulse - no need to simulate
                break
            end
        end

        # now up to rep rate do the other bit
        while t < pulse_interval
            mc_step!(lattice, n, pulse, sim.pulse_params, t, sim.dt2,
                     false, counts, sim.binwidth, print_debug)
            t += sim.dt2
            if sum(n) == 0
                # no excitations - no need to simulate
                break
            end
        end
        curr_maxcount = maximum(counts[ec..., :])
        if rep // 100 == 0
            println("run number ", run,
                    " rep ", rep, " max count = ", curr_maxcount)
        end
        rep += 1
    end

    write_hist(bins, counts, labels, ec, hist_file)
    plot_counts(bins, transpose(counts), labels, maximum(counts),
                splitext(hist_file)[1] * ".png")
    counts
end

function setup(json_file)
    params = JSON.parsefile(json_file)
    pd = params["protein"]
    ld = params["lattice"]
    proteins = [make_lh2(pd)]
    lattice = lattice_generation(get_lattice(ld["type"]),
                                 ld["n_sites"], proteins, [1.0])
    pd = params["pulse"]
    pulse_params = PulseParams(pd["peak_time"],
        pd["fwhm"], pd["wavelength"],
        pd["fluence"])
    sd = params["simulation"]
    sim = SimulationParams(sd["tmax"], sd["dt1"], sd["rep_rate"] * 1e6,
                           sd["dt2"], sd["binwidth"], pulse_params,
                           sd["n_counts"], sd["n_repeats"])
    protein_names = [p.name for p in proteins]
    outdir = joinpath("out",
        join([protein_names..., "_", sim.rep_rate/1e6, "MHz_fluence_",
              pulse_params.f, "_T_T_rate_", 
              params["protein"]["T_T_transfer"], "_"]))
    mkpath(outdir)
    (sim, lattice, outdir)
end

function run(json_file, seed_start = 0)

    sim, lattice, outdir = setup(json_file)
    lattice_plot_file = joinpath(outdir, "lattice.svg")
    pulse_file = joinpath(outdir, "pulse.txt")
    hist_path = joinpath(outdir, "hist")

    pulse = construct_pulse(sim.pulse_params, sim.dt1)
    x = (0:ceil((2.0 * sim.pulse_params.μ)/sim.dt1)) * sim.dt1

    open(pulse_file, "w") do io
        writedlm(io, hcat(x, pulse))
    end

    plot_lattice(lattice, lattice_plot_file)

    (bins, total_counts, labels, ec) = generate_histogram(sim, lattice)

    run = 1
    while run <= sim.repeats
        hist_file = "$(hist_path)_$(run).txt"
        total_counts += one_run(sim, lattice, seed_start + run, hist_file)
        run += 1
    end

    total_hist_file = "$(hist_path)_total.txt"
    write_hist(bins, total_counts, labels, ec, total_hist_file)
    plot_counts(bins, transpose(total_counts), labels,
                maximum(total_counts),
                splitext(total_hist_file)[1] * ".png")
    
    (bins, total_counts, labels, ec)
end

end
