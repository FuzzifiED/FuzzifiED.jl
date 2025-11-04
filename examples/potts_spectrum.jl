# This example calculates the spectrum of 3d Potts model on fuzzy sphere.
# This example reproduces Table I and Figure 4 in arXiv : 2501.14320
# On my table computer, this calculation takes 1.544 s

using FuzzifiED
using LinearAlgebra
≈(x, y) = abs(x - y) < eps(Float32)
const ω = (-1 + √3im) / 2
const mat_Z = diagm([1, ω, ω'])
const mat_X = [ 0 0 1 ; 1 0 0 ; 0 1 0 ]

nm = 8
nf = 3
no = nm * nf
qnd = [
    GetNeQNDiag(no)
    GetLz2QNDiag(nm, nf)
]
qnf = [
    GetRotyQNOffd(nm, nf),
    GetFlavPermQNOffd(nm, nf, [[2,3]]),
    GetFlavPermQNOffd(nm, nf, [[1,2,3]], 3)
]
tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, nf, [1, 0.2])
    - GetDenIntTerms(nm, nf, [1, 0.2], mat_Z)
    - 0.54 * GetPolTerms(nm, nf, mat_X + mat_X')
)
tms_l2 = SimplifyTerms(GetL2Terms(nm, nf))

cfs = Confs(no, [nm, 0], qnd)

result = []
for R in [1, -1], (Z2, Z3) in [(1, 1), (-1, 1), (0, ω)]
    bs = Basis(cfs, [R, Z2, Z3], qnf)
    hmt = Operator(bs, tms_hmt)
    hmt_mat = OpMat(hmt)
    enrg, st = GetEigensystem(hmt_mat, 20)

    l2 = Operator(bs, tms_l2)
    l2_mat = OpMat(l2)
    l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)]

    for i in eachindex(enrg)
        push!(result, real.([enrg[i], l2_val[i], Z2, Z3 == 1 ? 0 : 1]))
    end
end

sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_T = filter(st -> st[2] ≈ 6 && st[3] ≈ 1 && st[4] ≈ 0, result)[1][1]
spec = [ round.([ 3 * (st[1] - enrg_0) / (enrg_T - enrg_0) ; st] .+ √eps(Float64), digits = 6) for st in result ]
display(permutedims(hcat(spec...)))