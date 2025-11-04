# This example calculates the spectrum of free Majorana fermion.
# This example reproduces Figure 5 in arXiv : 2509.08038
# On my portable computer, this calculation takes 7.718 s

using FuzzifiED
using FuzzifiED.Fuzzifino
FuzzifiED.ElementType = Float64
≈(x, y) = abs(x - y) < √eps(Float64)

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
for lz = 0 : 1
    bs = SBasis(cfs[lz])
    hmt = SOperator(bs, tms_hmt)
    hmt_mat = OpMat(hmt)
    enrg, st = GetEigensystem(hmt_mat, 20)

    l2 = SOperator(bs, tms_l2)
    l2_mat = OpMat(l2)
    l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)] 

    for i in eachindex(enrg)
        push!(result, [enrg[i], l2_val[i]])
    end
end
sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_T = filter(st -> st[2] ≈ 6, result)[1][1]
spec = [ round.([ 3 * (st[1] - enrg_0) / (enrg_T - enrg_0) ; st] .+ √eps(Float64), digits = 6) for st in result ] 
display(permutedims(hcat(spec...)))
