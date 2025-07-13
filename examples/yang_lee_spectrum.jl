# This example calculates the spectrum of 3d non-unitary Yang-Lee CFT on fuzzy sphere at nm = 12.
# This example reproduces the finite size data of Figure 11 in arXiv : 2505.06369.
# On my portable computer, this calculation takes 10.714 s

using FuzzifiED
const σ1 = [  1  0 ;  0  0 ]
const σ2 = [  0  0 ;  0  1 ]
const σz = [  1  0 ;  0 -1 ]
const σx = [  0  1 ;  1  0 ]
FuzzifiED.ElementType = ComplexF64
≈(x, y) = abs(x - y) < eps(Float32)

nm = 12
qnd = [ 
    GetNeQNDiag(2 * nm), 
    GetLz2QNDiag(nm, 2) ]
qnf = [ 
    GetParityQNOffd(nm, 2, [2, 1], [-1, 1]), 
    GetRotyQNOffd(nm, 2) ]

tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, 2, 2 .* [4.75, 1.0], σ1, σ2)
    - 10.0 * GetPolTerms(nm, 2, σx) 
    + 4.41im * GetPolTerms(nm, 2, σz))
tms_l2 = GetL2Terms(nm, 2)

cfs = Confs(2 * nm, [nm, 0], qnd)

result = []
for P in [1, -1], R in [1, -1]
    bs = Basis(cfs, [P, R], qnf)
    hmt = Operator(bs, tms_hmt ; red_q = 0, sym_q = 2)
    hmt_mat = OpMat(hmt)
    enrg, st = GetEigensystem(hmt_mat, 20)

    l2 = Operator(bs, tms_l2)
    l2_mat = OpMat(l2)
    l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)]

    for i in eachindex(enrg)
        push!(result, round.([enrg[i], l2_val[i], P], digits = 6))
    end
end

sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_T = filter(st -> st[2] ≈ 6 && st[3] ≈ 1, result)[2][1]
result_dim = [ [ 3 * (st[1] - enrg_0) / (enrg_T - enrg_0) ; st] for st in result ]
for P in [1, -1]
    display(permutedims(hcat(
        filter(st -> st[4] ≈ P, result_dim)...
    )))
end