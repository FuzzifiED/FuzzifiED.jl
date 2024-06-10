using FuzzifiED
using ITensors
using ITensorMPOConstruction
const σx = [  0  1 ;  1  0 ]
const σz = [  1  0 ;  0 -1 ]

function MyMPO(os, sites)
    operatorNames = [ "I", "C", "Cdag", "N" ]
    opCacheVec = [ [OpInfo(ITensors.Op(name, n), sites[n]) for name in operatorNames] for n in eachindex(sites)  ]
    return MPO_new(os, sites ; basisOpCacheVec = opCacheVec)
end

nm = 12
nf = 2
no = nm * nf

path = "nm_$(nm)/"
mkpath(path)

ps_pot = [4.75, 1.] ./ 2
tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, 2 ; ps_pot) - 
    GetDenIntTerms(nm, 2 ; ps_pot, mat_a = σx) - 
    3.16 * GetPolTerms(nm, 2 ; mat = σz)
)
qnd = [ 
    GetNeQNDiag(no), 
    GetLz2QNDiag(nm, nf), 
    GetZnfChargeQNDiag(nm, nf) 
]
hmt, sites = GetMPOSites("hmt", tms_hmt, qnd ; path, mpo_method = MyMPO)

cf0 = [ isodd(o) ? 1 : 0 for o = 1 : no ]
st0 = MPS(sites, string.(cf0))

Eg, stg = EasySweep("g", hmt, st0 ; path)
Ee, ste = EasySweep("e", hmt, st0 ; path, proj = ["g"])

cf1 = cf0
cf1[1] = 0
cf1[2] = 1
st1 = MPS(sites, string.(cf1))
Es, sts = EasySweep("s", hmt, st1 ; path)

tms_l2 = GetL2Terms(nm, 2)
l2 = GetMPO("l2", tms_l2, sites ; path, mpo_method = MyMPO)
@show inner(stg', l2, stg)