"""
- pigments = list of names of pigments
- states = list of names of states on those pigments
- ``n_p`` = `length(pigments)`
- ``n_s`` = `length(states)`
- ps = vector (``n_s``) labelling which pigment each state is on
- dist = matrix (``n_s, n_s``) - which states are distinguishable for annihilation
- ``n_{tot}`` = vector (``n_p``) = number of each pigment per protein
- ``n_{thermal}`` = vector (``n_p``) = number of each pigment whose excited
  state is thermally available
- hop = vector of length ``n_s`` = intercomplex hopping rate per state
- intra = matrix (``n_s``, ``n_s``) = decay rates on the diagonal,
  transfer rates between states on the off-diagonal
- ann = matrix (``n_s``, ``n_s``) = annihilation rates between states
- `which_ann` = matrix (``n_s``, ``n_s``) = which species (index) is annihilated
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
mutable struct Protein
    name
    pigments
    states
    nₚ
    nₛ
    ps
    dist
    n_tot
    n_thermal
    hop
    intra
    ann
    which_ann
    xsec
    emissive
end

function intersystem_rate(e_chl, e_car, k_init)
    # detailed balance
    # this is not correct but it's something like this, check
    k = k_init * (e_chl / e_car)
end

"""
note - this will need changing so that we can mess with the
carotenoid rates more easily.
"""
function get_protein(name)
    if name == "lh2"
        p = Protein("LH2",
            ["BChl", "Car"],
            ["BChl_S", "BChl_T", "Car_T"],
            2, 3, [1, 1, 2], # nₚ, nₛ, which pigment is each state on
            # following line is distinguishability - Chl states aren't
            [[false false true]; [false false true]; [true true false]],
            [20, 4], # number of pigments total
            [20, 4], # number of states accessible thermally
            [1.0/10e-12, 0.0, 0.0], # intercomplex hopping rates
            [[1.0/1e-9 1.0/15e-9 0.0]; # intracomplex transfer
             [0.0 1.0/1e-7 1.0/15e-9]; # diagonal -> decay rate
             [0.0 0.0 1.0/1e-7]],
            [[1.0/16e-9 1.0/16e-9 0.0]; # annihilation rates
             [1.0/16e-9 1.0/16e-9 0.0];
             [0.0 0.0 1.0/16e-9]],
            [[1 1 0]; [1 1 0]; [0 0 3]], # which state gets annihilated
            [1e-16, 0.0, 0.0], # cross-section of each state
            [true, false, false] # which decays are emissive
    )
    elseif name == "lhcii"
        p = Protein("LHCII",
            ["Chl", "Car"],
            ["Chl_S", "Car_S"],
            2, 2, [1, 2], 
            [[false true]; [true false]],
            [14, 1],
            [5, 1],
            [1.0/10e-12, 0.0],
            [[1.0/4e-9 1.0/100.0e-12];
             [0.0 1.0/10.0e-12]],
            [[1.0/16e-12 0.0];
             [0.0 0.0]],
            [[1 0]; [0 1]],
            [2.74e-16, 0.0],
            [true, false]
    )
    elseif name == "chl"
        p = Protein("Chl",
            ["Chl"],
            ["Chl_S"],
            1, 1, [1], 
            false * ones(Bool, 1, 1),
            [20],
            [20],
            [1.0/10e-12],
            1.0/4e-9 * ones(1,1),
            1.0/16e-12 * ones(1,1),
            1 * ones(1,1),
            [2.74e-16, 0.0],
            [true]
    )
    elseif name == "chl_det"
        p = Protein("Chl",
            ["Chl"],
            ["Chl_S"],
            1, 1, [1], 
            false * ones(Bool, 1, 1),
            [20],
            [20],
            [0.0],
            1.0/4e-9 * ones(1,1),
            1.0/16e-12 * ones(1,1),
            1 * ones(1,1),
            [2.74e-16, 0.0],
            [true]
    )
    end
    p
end
