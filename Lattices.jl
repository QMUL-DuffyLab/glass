module Lattices

using Proteins, LinearAlgebra, Graphs, GraphPlot

# can use this to define arbitrary lattices
struct LatticeParams
    basis::Vector{Vector{Real}}
    radius::Real
    coordination::Integer
end

struct Site
    r::Vector{Real}
    p::Protein
end

struct Lattice
    params::LatticeParams
    sites::Vector{Site}
    A::Matrix{Integer}
    nn::Matrix{Integer}
end

neighbours(p1, p2, spacing) = isapprox(norm(p1 - p2), spacing, atol=1e-6)

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

"""
generate a lattice of size nmax according to params.
iterate over current sites, add if the current point is new.
then loop again to add indices of neighbours.
"""
function lattice_generation(params::LatticeParams, nmax::Int, proteins::Vector{Protein}, rho)
    lv = lattice_vectors(params)
    nn = zeros(Integer, params.coordination, nmax)
    spacing = norm(lv[:, 1])
    ps = sample(proteins, Weights(rho), nmax)
    sites = Vector{Site}(undef, nmax)
    sites[1] = Site(params.basis[1], ps[1])
    adj = zeros(Integer, nmax, nmax)
    n_sites = 1
    while n_sites < nmax
        for site in sites
            for n in eachcol(lv)
                for b in params.basis
                    ri = site.r + n + b
                    add = true
                    for (i, site2) in enumerate(sites[1:n_sites])
                        if isapprox(ri, site2.r, atol=1e-6)
                            add = false
                            break
			end
                    end
                    if add && n_sites < nmax n_sites += 1
			sites[n_sites] = Site(ri, ps[n_sites])
                    end
		end
            end
	end
    end
    # in theory there should be a way to generate neighbour data while
    # adding sites, but in practice it's easier just to do this
    for (i, s1) in enumerate(sites)
        k = 1
        for (j, s2) in enumerate(sites)
            if neighbours(s1.r, s2.r, spacing)
                adj[i, j] = 1
                nn[k, i] = j
                k += 1
            end
        end
    end
    Lattice(params, sites, adj, nn)
end

function plot_lattice(l::Lattice)
    g = SimpleGraph(l.A)
    xs = [s.r[1] for s in l.sites]
    ys = [s.r[2] for s in l.sites]
    gplot(g, xs, ys)
end

hex = LatticeParams([[0.0, 0.0]], 1.0, 6)
square = LatticeParams([[0.0, 0.0]], 1.0, 4)
line = LatticeParams([[0.0, 0.0]], 1.0, 2)
# honeycomb doesn't work atm. need to fix
honeycomb = LatticeParams([[0.0, 0.0], [0.0, -1.0]], 1.0, 6)

end
