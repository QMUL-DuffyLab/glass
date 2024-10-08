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
            [27, 9], # number of pigments total
            [18, 9], # number of states accessible thermally
            [1.0/10e-12, 0.0, 0.0], # intercomplex hopping rates
            [[1.0/1e-9 1.0/15e-9 0.0]; # intracomplex transfer
             [0.0 1.0/1e-7 1.0/1e-12]; # diagonal -> decay rate
             [0.0 0.0 1.0/1e-7]],
            [[1.0/16e-12 1.0/16e-12 1.0/16e-12]; # annihilation rates
             [1.0/16e-12 1.0/16e-12 0.0];
             [1.0/16e-12 0.0 1.0/16e-12]],
            [[1 1 1]; [1 1 0]; [1 0 3]], # which state gets annihilated
            [(1.0/27.0)*1e-14, 0.0, 0.0], # cross-section of each state
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

"""
convenience function to change LH2 rates based on carotenoid energy.
does a quick detailed balance calculation to get the forward and backward
rates of transfer between the chlorophyll and carotenoid triplets.
assumes (among other things) that BChl singlets can annihilate with
themselves, BChl triplets and Car triplets, but BChl triplets cannot
annihilate with Car triplets. It's possible to change various other
parameters here - hopefully it's clear what they all are.
singlet xsec here taken from σ(800) in paper, which i assume is for
one LH2 complex
"""
function make_lh2(t_t_time;
    hop_time = 10.0e-12,
    isc_time = 18.6e-9,
    ann_time = 16.0e-12,
    bchls_decay = 1.14e-9,
    bchlt_decay = 5.5e-6,
    cart_decay = 5.5e-6,
    xsec = 1e-14
    )
    hop_rates = 1.0 ./ [hop_time, 0., 0.]
    intra_rates = (1.0 ./ 
                   [[bchls_decay isc_time 0.];
                    [0.0 bchlt_decay t_t_time]; 
                    [0.0 0.0 cart_decay]])
    ann_rates = (1.0 ./ 
                   [[ann_time ann_time ann_time];
                    [ann_time ann_time 0.0]; 
                    [0.0 0.0 ann_time]])
    p = Protein("LH2",
        ["BChl", "Car"],
        ["BChl_S", "BChl_T", "Car_T"],
        2, 3, [1, 1, 2], # nₚ, nₛ, which pigment is each state on
        # following line is distinguishability - Chl states aren't
        [[false false true]; [false false true]; [true true false]],
        [27, 9], # number of pigments total
        [18, 9], # number of states accessible thermally
        hop_rates,
        intra_rates,
        ann_rates,
        [[1 1 1]; [1 1 0]; [1 0 3]], # which state gets annihilated
        [xsec, 0.0, 0.0], # cross-section of each state
        [true, false, false] # which decays are emissive
    )
end

function make_lh2(dict::Dict)
    hop_rates = 1.0 ./ [dict["hop"], 0., 0.]
    intra_rates = (1.0 ./ 
                   [[dict["BChl_S_decay"] dict["ISC"] 0.];
                    [0.0 dict["BChl_T_decay"] dict["T_T_transfer"]]; 
                    [0.0 0.0 dict["Car_T_decay"]]])
    s_s_ann = dict["BChl_S_BChl_S_ann"]
    s_t_ann = dict["BChl_S_BChl_T_ann"]
    s_car_ann = dict["BChl_S_Car_T_ann"]
    t_car_ann = dict["BChl_T_Car_T_ann"]
    t_t_ann = dict["BChl_T_BChl_T_ann"]
    car_car_ann = dict["Car_T_Car_T_ann"]
    ann_rates = (1.0 ./ 
                 [[s_s_ann s_t_ann s_car_ann];
                    [s_t_ann t_t_ann t_car_ann]; 
                    [s_car_ann t_car_ann car_car_ann]])
    p = Protein("LH2",
        ["BChl", "Car"],
        ["BChl_S", "BChl_T", "Car_T"],
        2, 3, [1, 1, 2], # nₚ, nₛ, which pigment is each state on
        # following line is distinguishability - Chl states aren't
        [[false false true]; [false false true]; [true true false]],
        [27, 9], # number of pigments total
        [18, 9], # number of states accessible thermally
        hop_rates,
        intra_rates,
        ann_rates,
        [[1 1 1]; [1 1 0]; [1 0 3]], # which state gets annihilated
        [dict["xsec"], 0.0, 0.0], # cross-section of each state
        [true, false, false] # which decays are emissive
    )
end
