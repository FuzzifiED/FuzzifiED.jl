# This example calculates the spectrum of 3d Ising model on fuzzy sphere
# for bosons at fractional filling ν = 1/2.
# This example reproduces Figure 12a,b in Phys. Rev. X 15, 031007 (2025).
# On my portable computer, this calculation takes 10.186 s.
# We acknowlege Cristian Voinea for his help in reproducing the results. 

using FuzzifiED
using FuzzifiED.Fuzzifino
const σ1 = [ 1 0 ; 0 0 ]
const σ2 = [ 0 0 ; 0 1 ]
const σ0 = [ 1 0 ; 0 1 ]
const σx = [ 0 1 ; 1 0 ]
FuzzifiED.ElementType = Float64
≈(x, y) = abs(x - y) < eps(Float32)

ne = 7
nm = 2 * ne - 1
nf = 2 
nof = 1 # Fuzzifino can only deal with mixture of bosons and fermions, so we put a single site of fermion and keep it empty.
nob = nm * nf 
qnd = [
    GetNeSQNDiag(nof, nob),
    SQNDiag(GetNeQNDiag(nof), nob),
    GetBosonLz2SQNDiag(nof, nm, nf)
]
qnf = [
    GetBosonFlavPermSQNOffd(nof, nm, nf, [2, 1]),
    GetBosonRotySQNOffd(nof, nm, nf)
]
cfs = SConfs(nof, nob, ne, [ne, 0, 0], qnd) ;

tms_hmt = SimplifyTerms(
    GetBosonDenIntSTerms(nm, 2, [1.0])
    + GetBosonDenIntSTerms(nm, 2, 2 .* [0.0, 0.53, 0.09],  σ1, σ2)
    - 0.25 * GetBosonPolSTerms(nm, 2, σx) 
)
tms_l2 = GetBosonL2STerms(nm, nf) ; 

result = []
for Z in [1, -1], R in [1, -1]
    bs = SBasis(cfs, [Z, R], qnf)
    hmt = SOperator(bs, tms_hmt)
    hmt_mat = OpMat(hmt)
    enrg, st = GetEigensystem(hmt_mat, 20)

    l2 = SOperator(bs, tms_l2)
    l2_mat = OpMat(l2)
    l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)]

    for i in eachindex(enrg)
        push!(result, round.([enrg[i], l2_val[i], Z], digits = 6))
    end
end

sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_T = filter(st -> st[2] ≈ 6 && st[3] ≈ 1, result)[1][1]
result_dim = [ [ 3 * (st[1] - enrg_0) / (enrg_T - enrg_0) ; st] for st in result ]
display(permutedims(hcat(result_dim...)))
