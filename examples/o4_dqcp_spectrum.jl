# This example calculates the spectrum of O(4) DQCP on fuzzy sphere.
# This example reproduces Figure 4 and Table I in arXiv : 2507.01322
# On my portable computer, this calculation takes 22.629 s

using FuzzifiED
const ⊗ = kron
const σ0 = [  1  0 ;  0  1 ]
const σx = [  0  1 ;  1  0 ]
const σy = [  0 im ;-im  0 ]
const σz = [  1  0 ;  0 -1 ]
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
    GetFlavQNDiag(nm, nf, Dict(2=>1, 4=>-1))
]
qnf = [
    GetParityQNOffd(nm, nf, [3,4,1,2], Dict(1=>-1, 2=>-1)),
    GetRotyQNOffd(nm, nf),
    GetFlavPermQNOffd(nm, nf, [[1,3]], Dict(1=>-1)),
    GetFlavPermQNOffd(nm, nf, [[2,4]], Dict(2=>-1)),
    GetFlavPermQNOffd(nm, nf, [[1,2],[3,4]])
]

ps_pot_u0 = [ 1.0 ]
ps_pot_u1 = [-0.6071 ]
ps_pot_u5 = [ 0.75 * 0.6071 ]
tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, nf, ps_pot_u0) +
    GetDenIntTerms(nm, nf, ps_pot_u1, [σ0⊗σx,σx⊗σy,σy⊗σy,σz⊗σy]) +
    GetDenIntTerms(nm, nf, ps_pot_u5, [σ0⊗σz])
)
tms_l2 = GetL2Terms(nm, nf)
tms_c2 = GetC2Terms(nm, nf, [σx⊗σz,σy⊗σz,σz⊗σz,σx⊗σ0,σy⊗σ0,σz⊗σ0])
result = []

cfs = Confs(no, [ne, 0, 0, 0], qnd)
for P in (1,-1), R in (1,-1), (Z1,Z2) in ((1,1),(1,-1),(-1,-1))
    bs = Basis(cfs, [P, R, Z1, Z2, 0], qnf)
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
        push!(result, [enrg[i], l2_val[i], c2_val[i], P])
    end
end

cfs = Confs(no, [ne, 0, 1, 1], qnd)
for P in (1,-1), R in (1,-1), X in (1,-1)
    bs = Basis(cfs, [P, R, 0, 0, X], qnf)
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
        push!(result, [enrg[i], l2_val[i], c2_val[i], P])
    end
end

sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_T = filter(st -> st[2] ≈ 6 && st[3] ≈ 0 && st[4] ≈ 1, result)[1][1]
spec = [ round.([ 3 * (st[1] - enrg_0) / (enrg_T - enrg_0) ; st ] .+ √eps(Float64), digits = 6) for st in result ]
display(permutedims(hcat(spec...)))
