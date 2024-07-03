### A Pluto.jl notebook ###
# v0.19.43

using Markdown
using InteractiveUtils

# ╔═╡ 8309ef8f-23d5-49a6-9cc2-cb210302c555
using LinearAlgebra

# ╔═╡ d536667a-0b7a-45ce-9471-5fe8d1b654d1
using Graphs

# ╔═╡ 38ecc0bf-1fb3-4c5e-b7d9-5ec608277104
using GraphPlot

# ╔═╡ 18eb45f4-c158-4080-a0e9-dcc5e01c479e
# this might be overkill, honestly, but we can use it to define arbitrary lattices
struct LatticeParams
	basis::Vector{Vector{Real}}
	radius::Real
	coordination::Integer
end

# ╔═╡ 19bbc9c5-18e8-43af-a307-90e133fe5200
struct Lattice
	params::LatticeParams
	sites::Matrix{Real}
	A::Matrix{Integer}
	nn::Matrix{Integer}
end

# ╔═╡ 2aaf887d-6984-4403-a8f6-bc1d36710899
neighbours(p1, p2, spacing) = isapprox(norm(p1 - p2), spacing, atol=1e-6)

# ╔═╡ f3955bc5-f06f-4956-b5aa-1f7f0cd29f19
"""
generate the near-neighbour lattice vectors for a given set of
lattice parameters. need to check this for honeycomb?
"""
function lattice_vectors(params)
	nn = zeros(Real, 2, params.coordination)
	for i=1:params.coordination
		ϕ = 2π * i/params.coordination
		r = [[cos(ϕ) sin(ϕ)]; [-sin(ϕ) cos(ϕ)]] * [params.radius, 0]
		nn[:, i] = r
	end
	return nn
end

# ╔═╡ 29ce2bef-bd26-4d96-971f-fc44acb687eb
"""
generate a lattice of size nmax according to params.
iterate over current sites and add if the current point is new.
check for neighbours while checking overlap; we have to iterate twice anyway.
"""
function lattice_generation(params::LatticeParams, nmax::Int)
	lv = lattice_vectors(params)
	nn = zeros(Integer, params.coordination, nmax)
	spacing = norm(lv[:, 1])
	sites = zeros(Real, 2, nmax)
	sites[:, 1] = params.basis[1]
	adj = zeros(Integer, nmax, nmax)
	n_sites = 1
	while n_sites < nmax
		for site in eachcol(sites)
			for n in eachcol(lv)
				for b in params.basis
					ri = site + n + b
					add = true
					for (i, site2) in enumerate(eachcol(sites[:, 1:n_sites]))
						if isapprox(ri, site2, atol=1e-6)
							add = false
							break
						end
					end
					if add && n_sites < nmax
						n_sites += 1
						sites[:, n_sites] = ri
					end
				end
			end
		end
	end
	for (i, s1) in enumerate(eachcol(sites))
		k = 1
		for (j, s2) in enumerate(eachcol(sites))
			if neighbours([s1[1], s1[2]], [s2[1], s2[2]], spacing)
				adj[i, j] = 1
				nn[k, i] = j
				k += 1
			end
		end
	end
	Lattice(params, sites, adj, nn)
end

# ╔═╡ d4fb1c0c-5fd5-4fb9-809a-741c1952095b
function plot_lattice(l::Lattice)
	g = SimpleGraph(l.A)
	gplot(g, l.sites[1, :], l.sites[2, :])
end

# ╔═╡ 97eadd96-16a8-4b63-81a3-09cc59cecbf8
begin
	hex = LatticeParams([[0.0, 0.0]], 1.0, 6)
	square = LatticeParams([[0.0, 0.0]], 1.0, 4)
	line = LatticeParams([[0.0, 0.0]], 1.0, 2)
	# honeycomb doesn't work atm. need to fix
	honeycomb = LatticeParams([[0.0, 0.0], [0.0, -1.0]], 1.0, 6)

	lattice = lattice_generation(hex, 200)
	plot_lattice(lattice)
end

# ╔═╡ cd0a1835-cee5-43d5-b2f4-b09309363b4c
"""
now we need to define a protein.
parameters here are as follows:
- ``n_p`` = the total number of different pigments we're considering
- ``n_s`` = total number of different states on those pigments
- ``n_{tot}`` = vector of length ``n_p`` = number of each pigment per protein
- hop = vector of length ``n_s`` = intercomplex hopping rate per state
- intra = matrix (``n_s``, ``n_s``) = decay rates on the diagonal,
  transfer rates between states on the off-diagonal
- ann = matrix (``n_s``, ``n_s``) = annihilation rates between states
- xsec = vector (``n_s``)? unsure about this yet. need it to get absorption
  and stimulated emission rates
- emissive = bool vector to tell us which of the possible decay pathways are
  emissive, i.e. which we'd see in a TCSPC experiment. Unsure how long this will
  have to be as yet; requies a bit of thought

Putting all these together we can define a set of base rates for all the possible
processes that can happen on the protein, and then define which of the 
population losses are emissive; a lattice of these then gives us all the 
data we need to run our simulated TCSPC experiment.

There are some assumptions here - firstly, we assume that excitations can't hop
between complexes and transfer between states at the same time, which I think
should be reasonable. Also the treatment of absorption and stimulated emission
might need to be fleshed out, depending on what we want to do going forward.
"""
struct Protein
	nₚ::Integer
	nₛ::Integer
	n_tot::Vector{Integer}
	hop::Vector{Real}
	intra::Matrix{Real}
	ann::Matrix{Real}
	xsec::Vector{Real}
	emissive::Vector{Bool}
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
GraphPlot = "a2cc645c-3eea-5389-862e-a155d0052231"
Graphs = "86223c79-3864-5bf0-83f7-82e725a168b6"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[compat]
GraphPlot = "~0.6.0"
Graphs = "~1.9.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.4"
manifest_format = "2.0"
project_hash = "a0bd92f3b844cbfd97daa9564f0fad17013832bb"

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

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

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
# ╠═8309ef8f-23d5-49a6-9cc2-cb210302c555
# ╠═d536667a-0b7a-45ce-9471-5fe8d1b654d1
# ╠═38ecc0bf-1fb3-4c5e-b7d9-5ec608277104
# ╠═18eb45f4-c158-4080-a0e9-dcc5e01c479e
# ╠═19bbc9c5-18e8-43af-a307-90e133fe5200
# ╠═2aaf887d-6984-4403-a8f6-bc1d36710899
# ╠═f3955bc5-f06f-4956-b5aa-1f7f0cd29f19
# ╠═29ce2bef-bd26-4d96-971f-fc44acb687eb
# ╠═d4fb1c0c-5fd5-4fb9-809a-741c1952095b
# ╠═97eadd96-16a8-4b63-81a3-09cc59cecbf8
# ╠═cd0a1835-cee5-43d5-b2f4-b09309363b4c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
