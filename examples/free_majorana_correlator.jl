# This example calculates the 2-pt correlator of free Majorana fermion.
# This example reproduces Figure 7 in arXiv : 2508.*****. 
# On my portable computer, this calculation takes 3.338 s

using FuzzifiED
using FuzzifiED.Fuzzifino
using Jacobi
FuzzifiED.ElementType = Float64
≈(x, y) = abs(x - y) < √eps(Float64)
outround(v) = round.(v .+ 1E-7 ; digits = 6) 

let ovl_χηI

nmf = 10
nof = nmf 
nmb = nmf - 1
nob = nmb
qnd = [
    GetNeSQNDiag(nof, nob),
    GetBosonLz2SQNDiag(nof, nmb, 1) + SQNDiag(GetLz2QNDiag(nmf, 1), nob)
]

cfs = Dict{Int64, SConfs}()
for lz = 0 : 1 
    cfs[lz] = SConfs(nof, nob, nof, [nof, lz], qnd)
end 

amd_ff = GetFermionSMod(nmf, 1, 1) * GetFermionSMod(nmf, 1, 1)
amd_fb = GetFermionSMod(nmf, 1, 1) * GetBosonSMod(nmb, 1, 1)
amd_bb = GetBosonSMod(nmb, 1, 1) * GetBosonSMod(nmb, 1, 1)

tms_hop = ContractMod(amd_ff', amd_bb, nmf - 2)
tms_fb  = ContractMod(amd_fb', amd_fb, nmf - 3/2)
tms_bb  = ContractMod(amd_bb', amd_bb, nmf - 2)
tms_pol = STerms(GetPolTerms(nof, 1))

tms_hmt = SimplifyTerms(
    2.0 * tms_fb + 1.0 * tms_bb
    - 0.3 * (tms_hop + tms_hop')
    + 0.0 * tms_pol
) 
tms_l2 = GetL2STerms(nmf, 1, nmb, 1) 

result = []
sts = Dict{Int64, Matrix{ComplexF64}}()
for lz = 0 : 1
    bs = SBasis(cfs[lz])
    hmt = SOperator(bs, tms_hmt)
    hmt_mat = OpMat(hmt)
    enrg, st = GetEigensystem(hmt_mat, 3)
    sts[lz] = st
end

stI = sts[0][:, 1]
stχ = sts[1][:, 1]

FuzzifiED.ObsNormRadSq = nmf
obs_η = StoreComps(GetFermionSObs(nmf, 1, 1)' * GetBosonSObs(nmb, 1, 1))
cor_χ = Dict()
for l = 0.5 : nmf - 1.5
    st_ηI = SOperator(SBasis(cfs[0]), SBasis(cfs[1]), GetComponent(obs_η, l, -0.5)) * stI
    if (l == 0.5)
        ovl_χηI = stχ' * st_ηI
    end
    cor_χ[l] = st_ηI' * st_ηI / abs(ovl_χηI) ^ 2
end

cor_χ_fn(θ :: Float64) = sum([(2 * l + 1) * sin(θ / 2) / 2 * jacobi(cos(θ), l - 1/2, 1, 0) * cor for (l, cor) in cor_χ])
display(real(permutedims(hcat([[θ / π, cor_χ_fn(θ)] for θ = 0 : π / 10 : 2 * π]...))))

end