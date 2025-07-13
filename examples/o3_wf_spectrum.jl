# This example calculates the spectrum of O(3) Wilson-Fisher CFT
# using the bilayer Heisenberg model.
# This example reproduces Table I and Figure 2 in Phys. Rev. B 110, 115113 (2024).
# On my table computer, this calculation takes 2.742 s

using FuzzifiED
const ⊗ = kron
const σ0 = [  1  0 ;  0  1 ]
const σ1 = [  1  0 ;  0  0 ]
const σ2 = [  0  0 ;  0  1 ]
const σx = [  0  1 ;  1  0 ]
const σy = [  0 im ;-im  0 ]
const σz = [  1  0 ;  0 -1 ]
FuzzifiED.ElementType = Float64
≈(x, y) = abs(x - y) < eps(Float32)

nm = 6
nf = 4
no = nf * nm
ne = 2 * nm

qnd = [
    GetNeQNDiag(no), 
    GetLz2QNDiag(nm, nf), 
    GetFlavQNDiag(nm, nf, [1, -1, 1, -1])
]
qnf = [
    GetParityQNOffd(nm, nf, [3,4,1,2], Dict(1=>-1, 2=>-1))
    GetFlavPermQNOffd(nm, nf, [[1,3],[2,4]])
    GetFlavPermQNOffd(nm, nf, [[1,2],[3,4]])
    GetRotyQNOffd(nm, nf)
]
cfs = Confs(no, [ne, 0, 0], qnd)

ps_pot_u0 = [ 1.0 ]
ps_pot_u1 = [ 0.55 ]
ps_pot_u2 = [ 0, 0.19 ]
tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, nf, ps_pot_u0) +
    GetDenIntTerms(nm, nf, ps_pot_u2, [σ1⊗σx,σ1⊗σy,σ1⊗σz],[σ2⊗σx,σ2⊗σy,σ2⊗σz]) -
    GetDenIntTerms(nm, nf, ps_pot_u1, [σ1⊗σx,σ1⊗σy,σ1⊗σz,σ2⊗σx,σ2⊗σy,σ2⊗σz]) -
    4 * π * 0.225 * GetPolTerms(nm, nf, σx⊗σ0)
)
tms_l2 = GetL2Terms(nm, nf)
tms_c2 = GetC2Terms(nm, nf, [σ0⊗σx,σ0⊗σy,σ0⊗σz]) * 2

result = []
for P in (1,-1), Z in (1,-1), X in (1,-1), R in (1,-1)
    bs = Basis(cfs, [P, Z, X, R], qnf)
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
        push!(result, round.([enrg[i], l2_val[i], c2_val[i], P, Z], digits = 6))
    end
end

sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_T = filter(st -> st[2] ≈ 6 && st[3] ≈ 0 && st[4] ≈ 1 && st[5] ≈ 1, result)[1][1]
result_dim = [ [ 3 * (st[1] - enrg_0) / (enrg_T - enrg_0) ; st] for st in result ]
display(permutedims(hcat(result_dim...)))