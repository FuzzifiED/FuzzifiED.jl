# This example calculates the spectrum of O(2) Wilson-Fisher CFT.
# This example reproduces Table 1, etc. in arXiv:2604.18705.
# On my table computer, this calculation takes 8.823 s

using FuzzifiED
using LinearAlgebra
FuzzifiED.ElementType = Float64
≈(x, y) = abs(x - y) < √eps(Float64)

nm = 9
nf = 3
no = nm * nf 

qnd = [
    GetNeQNDiag(no),
    GetLz2QNDiag(nm, nf),
    GetFlavQNDiag(nm, nf, Dict([1 => 1, 2 => -1])) 
]
qnf = [
    GetRotyQNOffd(nm, nf),
    GetFlavPermQNOffd(nm, nf, Dict([1 => 2, 2 => 1]))
] 

tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, nf, [4.0, 1.0])
    - GetDenIntTerms(nm, nf, [4.0, 1.0], [0 0 1 ; 0 0 0 ; 0 1 0])
    - 5.796 * GetPolTerms(nm, nf, diagm([0, 0, 1])) 
) # The D factor is different from the paper due to normal ordering
tms_l2 = GetL2Terms(nm, nf)

cfs = Dict{Int64, Confs}()
for s = 0 : 3 
    cfs[s] = Confs(no, [nm, 0, s], qnd)
end

result = []
for (Sz, X) in [(0, 1), (0,-1), (1, 0), (2, 0), (3, 0)], R in [1, -1]
    bs = Basis(cfs[Sz], [R, X], qnf)
    hmt = Operator(bs, tms_hmt)
    hmt_mat = OpMat(hmt)
    enrg, st = GetEigensystem(hmt_mat, 10)

    l2 = Operator(bs, tms_l2)
    l2_mat = OpMat(l2)
    l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)]

    for i in eachindex(enrg)
        push!(result, [enrg[i], l2_val[i], Sz, X])
    end
end

sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_σ = filter(st -> st[2] ≈ 0 && st[3] ≈ 1, result)[1][1]
spec = [ round.([ 0.519088 * (st[1] - enrg_0) / (enrg_σ - enrg_0) ; st] .+ √eps(Float64), digits = 6) for st in result ]
display(permutedims(hcat(spec...)))
