# This example calculates the spectrum of Sp(3) CFT on fuzzy sphere.
# This example reproduces Table I in Phys. Rev. Lett. 135, 026504 (2025).
# On my portable computer, this calculation takes 3.551 s

using FuzzifiED
using LinearAlgebra
FuzzifiED.ElementType = Float64
≈(x, y) = abs(x - y) < eps(Float32)

nm = 5
nf = 6
no = nf * nm
ne = 2 * nm

qnd = [
    GetNeQNDiag(no), 
    GetLz2QNDiag(nm, nf), 
    GetFlavQNDiag(nm, nf, Dict(1=>1, 4=>-1)),
    GetFlavQNDiag(nm, nf, Dict(2=>1, 5=>-1), 2),
    GetFlavQNDiag(nm, nf, Dict(3=>1, 6=>-1), 3)
]
qnf = [
    GetRotyQNOffd(nm, nf),
    GetFlavPermQNOffd(nm, nf, [[1,4]], Dict(1=>-1)),
    GetFlavPermQNOffd(nm, nf, [[2,5]], Dict(2=>-1)),
    GetFlavPermQNOffd(nm, nf, [[3,6]], Dict(3=>-1)),
    GetFlavPermQNOffd(nm, nf, [[1,2],[4,5]])
]
cfs = Confs(no, [ne, 0, 0, 0, 0], qnd)

Ω = [ i - j == nf / 2 ? 1 : 0 for i = 1 : nf, j = 1 : nf ]
tms_l2 = GetL2Terms(nm, nf)
tms_c2 = SimplifyTerms(
    GetDenIntTerms(nm, nf, [isodd(m) ? -.5 : 0 for m = 1 : nm]) + 
    GetPairIntTerms(nm, nf, [isodd(m) ? -1.0 : 0 for m = 1 : nm], Ω)
)
push!(tms_c2, Term(ne * (ne + nf) / 4, [-1, -1]))
tms_c2 = SimplifyTerms(tms_c2)
tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, nf, [1.0, 0.26435143448887427, 0.06520453215339449]) -
    GetPairIntTerms(nm, nf, [0.3798431681470266, 0.0, -0.02191234730375343], Ω)
)
result = []
for R in (1, -1), (Z1, Z2, Z3) in ((1,1,1), (1,1,-1), (-1,-1,1), (-1,-1,-1)), X in (1,-1)
    bs = Basis(cfs, [R, Z1, Z2, Z3, X], qnf)
    hmt = Operator(bs, tms_hmt)
    hmt_mat = OpMat(hmt)
    nst = ((R, Z1, Z2, Z3, X) == (1, 1, 1, 1, 1)) ? 20 : 10
    enrg, st = GetEigensystem(hmt_mat, nst)

    l2 = Operator(bs, tms_l2)
    l2_mat = OpMat(l2)
    l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)]
    
    c2 = Operator(bs, tms_c2)
    c2_mat = OpMat(c2)
    c2_val = [ st[:, i]' * c2_mat * st[:, i] for i in eachindex(enrg)]
    
    for i in eachindex(enrg)
        push!(result, round.([enrg[i], l2_val[i], c2_val[i]], digits = 6))
    end
end

sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_T = filter(st -> st[2] ≈ 6 && st[3] ≈ 0, result)[1][1]
result_dim = [ [ 3 * (st[1] - enrg_0) / (enrg_T - enrg_0) ; st] for st in result ]
display(permutedims(hcat(result_dim...)))
