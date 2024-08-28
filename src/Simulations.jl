module Simulations
include("Lattices.jl")

using Distributed, Random, DelimitedFiles

struct PulseParams
    μ::Real
    σ::Real
    # unsure if we'll need this or not
    λ::Real
    f::Real
end

struct SimulationParams
    tmax::Real
    dt1::Real
    rep_rate::Real
    dt2::Real
    binwidth::Real
    pulse::PulseParams
    maxcount::Integer
    repeats::Integer
end

function construct_pulse(p, dt)
    tmax = 2.0 * p.μ
    pulse = [(1.0 / (p.σ * sqrt(2.0 * pi))) * 
             exp(-((i * dt) - p.μ)^2 / (sqrt(2.0) * p.σ)^2)
             for i=1:ceil(tmax / dt)]
    pulse *= p.f
end

"""
use some pretty simple criteria to determine whether a given type of move
is currently possible on site i.
the annihilation criteria's less trivial, not fixed it yet
"""
function possible_moves(l::Lattice, i::Integer, n, ft::Real)
    move_types = []
    # pulse has to be on and there has to be a state with a nonzero cross section
    nₛ= l.proteins[l.sites[i].ip].nₛ
    xsecs = l.proteins[l.sites[i].ip].xsec
    if ft > 0 && any(xsecs .> 0.0)
        append!(move_types, ["gen"])
    end
    # product of xsec and current population must be > 0 in order for a second
    # photon to be absorbed and cause stimulated emission
    # println("xsecs = ", xsecs)
    # println("i = ", " protein = ",
    #         l.proteins[l.sites[i].ip].name, " n = ", n[i, 1:nₛ])
    if ft > 0 && any(xsecs .* n[i, 1:nₛ] .> 0.0)
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
function propose_move(l::Lattice, i::Integer, n, ft::Real)
    
    rate = 0.0
    ip = l.sites[i].ip
    p = l.proteins[ip]
    if ip > 1
        loss_start_index = sum(p.nₛ * (p.nₛ + 2) for p in l.proteins[1:(ip - 1)])
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
        # rate = ft * p.xsec[s] * (p.n_tot[p.ps[s]] - n[i, s]) / p.n_tot[p.ps[s]]
        rate = ft * p.xsec[s] * (p.n_tot[p.ps[s]] - n[i, s])
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
        # this currently fails if there's only one species
        # because apparently A[1, :] is still a matrix for a 1x1
        # matrix but not for larger than that. fucking julia, man
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
        si = rand(findall(>(0), n[i, 1:p.nₛ]))
        sf = rand(1:p.nₛ)
        # ordering ensures annihilation counting is commutative
        pair = (si < sf) ? [si, sf] : [sf, si]
        # println("p.dist = ", p.dist, " type = ", typeof(p.dist),
        #         " p.ann = ", p.ann, " type = ", typeof(p.ann),
        #         " pair = ", pair, " p.dist[pair] = ", p.dist[pair],
        #         " p.dist[pair expl] = ", p.dist[pair[1], pair[2]])
        if p.dist[pair...]
            # distinguishable
            rate = n[i, pair[1]] * n[i, pair[2]] * p.ann[pair...]
        else
            if pair[1] == pair[2] && n[i, pair[1]] == 1 && n[i, pair[2]] == 1
                println("ILLEGAL ANNIHILATION")
            end
            neff = n[i, pair[1]] + n[i, pair[2]]
            rate = p.ann[pair...] * (neff * (neff - 1)) / 2.0
        end
        # one column for every combination
        loss_index = loss_start_index + (2 + (pair[1] - 1)) * p.nₛ + pair[2]
        ist = si
        fst = sf
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
end

"""
one monte carlo step. 
attempt one of each type of possible move on each site on average;
`propose_move` picks a possible move.

once a move's been picked, run it through Metropolis; if it's accepted,
carry it out and increment the histogram iff 1.) an excitation's been lost
and 2.) we are currently binning decays.
"""
function mc_step!(l::Lattice, n, pulse,
    t::Real, dt::Real, bin::Bool, counts::Matrix{Integer}, binwidth::Real)
    
    # current intensity of pulse
    pind = round(Int, t / dt) + 1
    if pind <= length(pulse)
        ft = pulse[pind]
    else 
        ft = 0
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
            # println(ft, " ", rate, " ", prob, " ", move_type, " ", loss_index, " ",
            #         isi, " ", ist, " ", fsi, " ", fst)
            if rand() < prob
                # carry out the move - this changes population, which
                # is why we have to check possible moves again above
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

function one_run()
    # length(Sys.cpu_info()) looks awful but should give a reasonable
    # guess of how many processes to start - if not you'll have to
    # check how many cores your machine's got and put that in explicitly
    # t = addprocs(length(Sys.cpu_info()) - 1)

    # @everywhere begin
    # begin
        # using Random
        # Random.seed!(myid())
    Random.seed!()
    nmax = 200
    max_count = 1000
    rep = 1
    # proteins = [get_protein("lh2"), get_protein("lhcii")]
    # rho = [0.5, 0.5]
    proteins = [get_protein("chl_det")]
    rho = [1.0]
    lattice = lattice_generation(get_lattice("hex"), nmax,
            proteins, rho)
    plot_lattice(lattice)
    pulse_params = PulseParams(200e-12, 50e-12, 485e-9, 1e14)
    sim = SimulationParams(15e-9, 1e-12, 1e6, 1e-9, 25e-12,
                           pulse_params, 2000, 1)
    pulse = construct_pulse(pulse_params, sim.dt1)
    println(pulse)
    (bins, counts, labels, ec) = generate_histogram(sim, lattice)
    n = zeros(Int, nmax, maximum([p.nₛ for p in proteins]))
    run = 1
    curr_maxcount = 0
    while run <= sim.repeats
        while curr_maxcount < sim.maxcount
            t = 0.0
            while t < sim.tmax
                mc_step!(lattice, n, pulse, t, sim.dt1, true, counts, sim.binwidth)
                t += sim.dt1
                if t >= 2.0 * pulse_params.μ && sum(n) == 0
                    break
                end
            end
            curr_maxcount = maximum(counts[ec..., :])
            println("run ", run, " rep ", rep, " cmc = ", curr_maxcount)
            rep += 1
            # now up to rep rate do the other bit
        end
        run += 1
    end

    open("out/bincounts.txt", "w") do io
        writedlm(io, hcat(bins, transpose(counts)))
    end
    
    (bins, counts, labels, ec)
end

end
