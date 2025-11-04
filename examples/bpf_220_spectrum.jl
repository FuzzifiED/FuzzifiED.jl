# This example calculates the spectrum of the transition between
# bosonic Pfaffian and Halperin 220 described by gauged Majorana fermion CFT
# This example reproduces Figure 4a in arXiv:2509.08036
# On my portable computer, this calculation takes 15.886 s.

using FuzzifiED
using FuzzifiED.Fuzzifino
const σ1 = [ 1 0 ; 0 0 ]
const σ2 = [ 0 0 ; 0 1 ]
const σ0 = [ 1 0 ; 0 1 ]
const σx = [ 0 1 ; 1 0 ]
FuzzifiED.ElementType = Float64
≈(x, y) = abs(x - y) < eps(Float32)

nm = 9
ne = nm + 1
nf = 2 
nof = 1
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
    GetBosonDenIntSTerms(nm, 2, [1.0], [σ1, σ2])
    + GetBosonDenIntSTerms(nm, 2, 2 .* [0.48],  σ1, σ2)
    - 0.58 * GetBosonPolSTerms(nm, 2, σx) 
)
tms_l2 = GetL2STerms(0, 0, nm, nf) ; 

result = []
for Z in [1], R in [1, -1]
    bs = SBasis(cfs, [Z, R], qnf)
    hmt = SOperator(bs, tms_hmt)
    hmt_mat = OpMat(hmt)
    enrg, st = GetEigensystem(hmt_mat, 20)

    l2 = SOperator(bs, tms_l2)
    l2_mat = OpMat(l2)
    l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)]

    for i in eachindex(enrg)
        push!(result, [enrg[i], l2_val[i], Z])
    end
end

sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_T = filter(st -> st[2] ≈ 6 && st[3] ≈ 1, result)[1][1]
spec = [ round.([ 3 * (st[1] - enrg_0) / (enrg_T - enrg_0) ; st] .+ √eps(Float64), digits = 6) for st in result if st[2] <= 12 ]
display(permutedims(hcat(spec...)))
