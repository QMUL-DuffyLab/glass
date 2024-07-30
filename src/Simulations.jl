module Simulations
include("Lattices.jl")

using Distributed, Random

struct PulseParams
    μ::Real
    σ::Real
    # unsure if we'll need this or not
    λ::Real
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

function construct_pulse(p::PulseParams, dt::Real)
    tmax = 2.0 * p.μ
    pulse = [(1.0 / (p.σ * sqrt(2.0 * pi))) * 
             exp(-1.0 * ((i * dt) - p.μ)^2 / 
                 (sqrt(2.0) * p.σ)^2) for i=1:ceil(tmax / dt)]
end

"""
use some pretty simple criteria to determine whether a given type of move
is currently possible on site i.
the annihilation criteria's less trivial, not fixed it yet
"""
function possible_moves(l::Lattice, i::Integer, n, ft::Real)
    move_types = []
    # pulse has to be on and there has to be a state with a nonzero cross section
    xsecs = l.proteins[l.sites[i].ip].xsec
    if ft > 0 && any(xsecs .> 0.0)
        append!(move_types, ["gen"])
    end
    # product of xsec and current population must be > 0 in order for a second
    # photon to be absorbed and cause stimulated emission
    if ft > 0 && any(xsecs .* n[i, :] .> 0.0)
        append!(move_types, ["se"])
    end
    if any(n[i, :] .> 0)
        append!(move_types, ["hop", "transfer", "decay"])
    end
    # this isn't always gonna work. need to think about it
    if sum(n[i, :] > 1)
        append!(move_types, "ann")
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
    
    move_types = []
    rate = 0.0
    ip = l.sites[i].ip
    protein = l.proteins[ip]
    if ip > 1
        loss_start_index = sum(p.nₛ * (p.nₛ + 2) for p in l.proteins[1:(ip - 1)])
    else
        loss_start_index = 0
    end
    isi = i
    fsi = i
    ist = 0
    fst = 0

    # get possible move types. this won't be empty because we check
    # whether there are any in the monte carlo step before this function's called
    move_types = possible_moves(l, i, n, ft)

    # pick a type of move
    mt = move_types[rand(1:length(move_types))]
    if mt == "gen"
        s = rand(p.xsec[p.xsec .> 0.0])
        rate = ft * p.xsec[s] * (p.n_tot[p.ps[s]] - n[i, s]) / p.n_tot[p.ps[s]]
        ist = s
        fst = s
    elseif mt == "se"
        s = rand(p.xsec[p.xsec .> 0.0])
        rate = ft * p.xsec[s] * n[i, s] / p.n_tot[p.ps[s]]
        loss_index = lost_start_index + s
        ist = s
        fst = s
    elseif mt == "hop"
        s = rand(1:p.nₛ)
        nn = rand(l.sites[i].nn)
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
        si = rand(1:p.nₛ)
        sf = rand(1:p.nₛ)
        while sf == si
            # this would be a decay process - reroll
            sf = rand(1:p.nₛ)
        end
        rate = n[i, si] * intra[si, sf] # is this right? check
        ist = si
        fst = sf
    elseif mt == "decay"
        s = rand(1:p.nₛ)
        rate = n[i, si] * intra[si, si]
        emissive = p.emissive[s]
        loss_index = loss_start_index + p.nₛ + s
        ist = s
        fst = s
    elseif mt == "ann"
        si = rand(1:p.nₛ)
        sf = rand(1:p.nₛ)
        # ordering ensures annihilation counting is commutative
        pair = (si < sf) ? [si, sf] : [sf, si]
        if p.dist[pair]
            # distinguishable
            rate = n[i, pair[1]] * n[i, pair[2]] * ann[pair]
        else
            neff = n[i, pair[1]] + n[i, pair[2]]
            rate = ann[pair] * (neff * (neff - 1)) / 2.0
        end
        # one column for every combination
        loss_index = loss_start_index + (2 + (pair[1] - 1)) * p.nₛ + pair[2]
        ist = si
        fst = sf
    end
    return (rate, move_type, loss_index, isi, ist, fsi, fst)
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
    if t < length(pulse) * dt
        ft = pulse[round(Int, t / dt) + 1]
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
            prob = rate * exp(-1.0 * rate * dt)
            if rand() < prob
                # carry out the move - this changes population, which
                # is why we have to check possible moves again above
                move!(n, move_type, isi, ist, fsi, fst, l.sites[s].p.which_ann[ist, fst])
                if bin && loss_index > 0
                    # increment the correct column of the histogram
                    counts[loss_index, floor(t / binwidth) + 1] += 1
                end
            end
        end
    end
    t += dt
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
    counts = zeros(Integer, n, n_losses)
    return (bins, counts, labels, emissive_columns)
end

function one_run()
    # length(Sys.cpu_info()) looks awful but should give a reasonable
    # guess of how many processes to start - if not you'll have to
    # check how many cores your machine's got and put that in explicitly
    # t = addprocs(length(Sys.cpu_info()) - 1)

    # @everywhere begin
    begin
        # using Random
        # Random.seed!(myid())
        Random.seed!()
        nmax = 200
        proteins = [get_protein("lh2"), get_protein("lhcii")]
        rho = [0.5, 0.5]
	lattice = lattice_generation(get_lattice("hex"), nmax,
	        proteins, rho)
	plot_lattice(lattice)
        pulse_params = PulseParams(50e-12, 200e-12, 485e-9)
        sim = SimulationParams(15e-9, 1e-12, 1e6, 1e-9, 25e-12,
                               pulse_params, 10000, 5)
        pulse = construct_pulse(pulse_params, sim.dt1)
        (bins, counts, labels, ec) = generate_histogram(sim, lattice)
        n = zeros(Int, nmax, sum([p.nₛ for p in proteins]))
        t = 0.0
        while t < sim.tmax
            mc_step!(lattice, n, pulse, t, sim.dt1, true, counts, sim.binwidth)
            if !(round(Int, t % sim.tmax / 100))
                print("t = ", t)
                print(counts)
            end
        end

    end
    
    # rmprocs(workers())
end

end
