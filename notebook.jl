### A Pluto.jl notebook ###
# v0.19.43

using Markdown
using InteractiveUtils

# ╔═╡ 948ebdca-a3ce-464e-98a3-637d1ea1a070
begin
	using Graphs, GraphPlot, StatsBase, Colors
	include("src/Lattices.jl")
end

# ╔═╡ 3f795568-9894-4292-88ab-2a85c28f2c68
using Distributed

# ╔═╡ b1b44fd9-127d-4009-9cc1-993b357a6690
"""
use some pretty simple criteria to determine whether a given type of move
is currently possible on site i.
the annihilation criteria's less trivial, not fixed it yet
"""
function possible_moves(l::Lattice, i::Integer, n::Matrix{Integer}, ft::Real)
	move_types = []
	# pulse has to be on and there has to be a state with a nonzero cross section
	if ft > 0 && any(l.sites[i].p.xsec .> 0.0)
		append!(move_types, ["gen"])
	end
	# product of xsec and current population must be > 0 in order for a second
	# photon to be absorbed and cause stimulated emission
	if ft > 0 && any(l.sites[i].p.xsec .* n[i, :] .> 0.0)
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

# ╔═╡ 038c9e9c-5252-4801-9259-b669d42e2ccd
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
function propose_move(l::Lattice, i::Integer, n::Matrix{Int}, ft::Real)
	
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
		pair = (si < sf) ? [si, sf] ? [sf, si]
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

# ╔═╡ 652428a5-4618-4986-bf87-cf7819b4359e
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

# ╔═╡ 9fbd7017-7bbb-4057-a609-1b079afb4207
"""
one monte carlo step. 
attempt one of each type of possible move on each site on average;
`propose_move` picks a possible move.

once a move's been picked, run it through Metropolis; if it's accepted,
carry it out and increment the histogram iff 1.) an excitation's been lost
and 2.) we are currently binning decays.
"""
function mc_step!(l::Lattice, n::Matrix{Integer}, pulse::Vector{Real},
	t::Real, dt::Real, bin::Bool, counts::Matrix{Integer}, binwidth::Real)
	
	# current intensity of pulse
	if t < length(pulse) * dt
		ft = pulse[int(t / dt) + 1]
	else 
		ft = 0
	end
	
	# pick each site once on average
	for i=1:length(l)
		s = rand(length(l))
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

# ╔═╡ 59b19e4d-ec01-46b2-8585-e4a82f467f13
struct Pulse
	μ::Real
	σ::Real
	# unsure if we'll need this or not
	λ::Real
end

# ╔═╡ 2640ff2c-003d-4831-a6d5-a2f0048bc96c
struct SimulationParams
	tmax::Real
	dt1::Real
	rep_rate::Real
	dt2::Real
	binwidth::Real
	pulse::Pulse
	maxcount::Integer
	repeats::Integer
end

# ╔═╡ caa3d2dc-dae8-4ef7-8013-57e55a0b8b9d
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

# ╔═╡ 724ebebb-540b-4a09-8d78-6339df8dfa5f
begin
	# various random testing bits
	xsec = [1., 0., 0.]
	n = [0, 1, 1]
	ann = [[1 1 1]; [0 2 0]; [0 0 3]]
	wh = findall(n .> 0)
	unique(collect(Iterators.product(wh, wh)))
	s = [10, 40, 37]
	rho = [0.8, 0.1, 0.1]
	sample([1:length(s);], Weights(rho), 50)

	# testing the actual code
	lattice = lattice_generation(get_lattice("hex"), 200, [get_protein("lh2"), get_protein("lhcii")], [0.5, 0.5])
	plot_lattice(lattice)
end

# ╔═╡ bdf32183-e30b-40ac-b861-8461d1ebb8c2
begin
	pulse = Pulse(50e-12, 200e-12, 485e-9)
	sim = SimulationParams(15e-9, 1e-12, 1e6, 1e-9, 25e-12, pulse, 10000, 5)
	(bins, counts, labels, ec) = generate_histogram(sim, lattice)
end

# ╔═╡ e4163c69-0c72-4b02-8b17-0bc9f1a0ce57
one_run()

# ╔═╡ c67abef7-0a18-4a42-88fe-45a1b04b88fc
using Random

# ╔═╡ 79b36702-8c8d-46a5-8f39-4486380b9a0e
function one_run()
	# length(Sys.cpu_info()) looks awful but should give a reasonable
	# guess of how many processes to start - if not you'll have to
	# check how many cores your machine's got and put that in explicitly
	t = addprocs(length(Sys.cpu_info()) - 1)

	@everywhere begin
		using Random
		Random.seed!(myid())
		pulse = Pulse(50e-12, 200e-12, 485e-9)
		sim = SimulationParams(15e-9, 1e-12, 1e6, 1e-9, 25e-12, pulse, 10000, 5)
		(bins, counts, labels, ec) = generate_histogram(sim, lattice)
	end
	
	rmprocs(workers())
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
Distributed = "8ba89e20-285c-5b6f-9357-94700520ee1b"
GraphPlot = "a2cc645c-3eea-5389-862e-a155d0052231"
Graphs = "86223c79-3864-5bf0-83f7-82e725a168b6"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
Colors = "~0.12.11"
GraphPlot = "~0.6.0"
Graphs = "~1.9.0"
StatsBase = "~0.34.3"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.3"
manifest_format = "2.0"
project_hash = "021557df86d641f870a293dec791dfb40999c8ba"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "b1c55339b7c6c350ee89f2c1604299660525b248"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.15.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Compose]]
deps = ["Base64", "Colors", "DataStructures", "Dates", "IterTools", "JSON", "LinearAlgebra", "Measures", "Printf", "Random", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "bf6570a34c850f99407b494757f5d7ad233a7257"
uuid = "a81c6b42-2e10-5240-aca2-a61377ecd94b"
version = "0.9.5"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.GraphPlot]]
deps = ["ArnoldiMethod", "Colors", "Compose", "DelimitedFiles", "Graphs", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "f76a7a0f10af6ce7f227b7a921bfe351f628ed45"
uuid = "a2cc645c-3eea-5389-862e-a155d0052231"
version = "0.6.0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "899050ace26649433ef1af25bc17a815b3db52b7"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.9.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "eeafab08ae20c62c44c8399ccb9354a04b80db50"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.7"

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

    [deps.StaticArrays.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "5cf7606d6cef84b543b483848d4ae08ad9832b21"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.3"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"
"""

# ╔═╡ Cell order:
# ╠═948ebdca-a3ce-464e-98a3-637d1ea1a070
# ╠═b1b44fd9-127d-4009-9cc1-993b357a6690
# ╠═038c9e9c-5252-4801-9259-b669d42e2ccd
# ╠═652428a5-4618-4986-bf87-cf7819b4359e
# ╠═9fbd7017-7bbb-4057-a609-1b079afb4207
# ╠═59b19e4d-ec01-46b2-8585-e4a82f467f13
# ╠═2640ff2c-003d-4831-a6d5-a2f0048bc96c
# ╠═caa3d2dc-dae8-4ef7-8013-57e55a0b8b9d
# ╠═724ebebb-540b-4a09-8d78-6339df8dfa5f
# ╠═bdf32183-e30b-40ac-b861-8461d1ebb8c2
# ╠═3f795568-9894-4292-88ab-2a85c28f2c68
# ╠═c67abef7-0a18-4a42-88fe-45a1b04b88fc
# ╠═79b36702-8c8d-46a5-8f39-4486380b9a0e
# ╠═e4163c69-0c72-4b02-8b17-0bc9f1a0ce57
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
