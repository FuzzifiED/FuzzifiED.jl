# This example calculates the spectrum of SO(5) DQCP on fuzzy sphere.
# This example reproduces Table II in Phys. Rev. X 14, 021044 (2024).
# On my portable computer, this calculation takes 9.852 s

using FuzzifiED
FuzzifiED.ElementType = Float64
≈(x, y) = abs(x - y) < eps(Float32)

nm = 7
nf = 4
no = nf * nm
ne = 2 * nm

qnd = [
    GetNeQNDiag(no), 
    GetLz2QNDiag(nm, nf), 
    GetFlavQNDiag(nm, nf, Dict(1=>1, 3=>-1)),
    GetFlavQNDiag(nm, nf, Dict(2=>1, 4=>-1), 2)
]
qnf = [
    GetParityQNOffd(nm, nf, [3,4,1,2], Dict(1=>-1, 2=>-1))
    GetFlavPermQNOffd(nm, nf, [[1,3]], Dict(1=>-1))
    GetFlavPermQNOffd(nm, nf, [[2,4]], Dict(2=>-1))
    GetFlavPermQNOffd(nm, nf, [2,1,4,3])
    GetRotyQNOffd(nm, nf)
]
cfs = Confs(no, [ne, 0, 0, 0], qnd)

ps_pot_u = [ 1.0]
ps_pot_v = [-0.9]
Ω = [ i - j == nf / 2 ? 1 : 0 for i = 1 : nf, j = 1 : nf ]
tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, nf, ps_pot_u) + 
    GetPairIntTerms(nm, nf, ps_pot_v, Ω)
)
tms_l2 = GetL2Terms(nm, nf)
tms_c2 = GetC2Terms(nm, nf, :Sp)

result = []
for P in (1,-1), (Z1, Z2, X) in ((1, 1, 1), (1, 1,-1), (1,-1, 0), (-1,-1, 1), (-1,-1,-1)), R in (1,-1)
    bs = Basis(cfs, [P, Z1, Z2, X, R], qnf)
    hmt = Operator(bs, tms_hmt)
    hmt_mat = OpMat(hmt)
    enrg, st = GetEigensystem(hmt_mat, 10)

    l2 = Operator(bs, tms_l2)
    l2_mat = OpMat(l2)
    l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)]
    
    c2 = Operator(bs, tms_c2)
    c2_mat = OpMat(c2)
    c2_val = [ st[:, i]' * c2_mat * st[:, i] for i in eachindex(enrg)]
    
    for i in eachindex(enrg)
        push!(result, round.([enrg[i], l2_val[i], c2_val[i], P], digits = 6))
    end
end

sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_T = filter(st -> st[2] ≈ 6 && st[3] ≈ 0 && st[4] ≈ 1, result)[1][1]
result_dim = [ [ 3 * (st[1] - enrg_0) / (enrg_T - enrg_0) ; st] for st in result ]
display(permutedims(hcat(result_dim...)))