# This example calculates the spectrum of magnetic line defect in 3d Ising model
# in lz = 0, P = ±1 and lz = 1 sectors, calibrated by bulk T.
# On my table computer, this calculation takes 2.246 s

using FuzzifiED
const σ1 = [  1  0 ;  0  0 ]
const σ2 = [  0  0 ;  0  1 ]
const σx = [  0  1 ;  1  0 ]
≊(x, y) = abs(x - y) < eps(Float32)

nm = 12
no = nm * 2
tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, 2 ; ps_pot = 2 .* [4.75, 1.], mat_a = σ1, mat_b = σ2) - 
    3.16 * GetPolTerms(nm, 2 ; mat = σx) )
    
qnd = [ 
    GetNeQNDiag(2 * nm), 
    GetLz2QNDiag(nm, 2) 
]
qnf = [ 
    GetParityQNOffd(nm, 2, [2, 1], [-1, 1]), 
    GetFlavPermQNOffd(nm, 2, [2, 1]), 
    GetRotyQNOffd(nm, 2) 
]
cfs = Confs(2 * nm, [nm, 0], qnd)
bs = Basis(cfs, [1, 1, 1], qnf)
hmt = Operator(bs, bs, tms_hmt ; red_q = 1, sym_q = 1)
hmt_mat = OpMat(hmt ; type = Float64)
enrg, st = GetEigensystem(hmt_mat, 6)
enrg_0 = enrg[1]
enrg_T = enrg[3]

result = []

qnd = [
    GetNeQNDiag(no),
    GetLz2QNDiag(nm, 2), 
    GetPinOrbQNDiag(no, [1, no - 1]), 
    GetPinOrbQNDiag(no, [2, no])
]
qnf = [
    GetParityQNOffd(nm, 2, [2, 1], [-1, 1]),
    GetRotyQNOffd(nm, 2)
]
qnf1 = [
    GetParityQNOffd(nm, 2, [2, 1], [1, -1]) * GetRotyQNOffd(nm, 2)
]
cfs = Confs(no, [nm, 0, 2, 0], qnd)
global enrg_d = 0
for P in (1, -1), R in (1, -1)
    bs = Basis(cfs, [P, R], qnf)
    hmt = Operator(bs, bs, tms_hmt ; red_q = 1, sym_q = 1)
    hmt_mat = OpMat(hmt ; type = Float64)
    enrg, st = GetEigensystem(hmt_mat, 10)
    if (P == 1 && R == 1) global enrg_d = enrg[1] end
    dim = (enrg .- enrg_d) ./ (enrg_T - enrg_0) * 3
    for i in eachindex(enrg)
        push!(result, round.([dim[i], enrg[i], 0, P, R], digits = 6))
    end
end

cfs = Confs(no, [nm, 2, 2, 0], qnd)
for PR in (1, -1)
    bs = Basis(cfs, [PR], qnf1)
    hmt = Operator(bs, bs, tms_hmt ; red_q = 1, sym_q = 1)
    hmt_mat = OpMat(hmt ; type = Float64)
    enrg, st = GetEigensystem(hmt_mat, 10)
    dim = (enrg .- enrg_d) ./ (enrg_T - enrg_0) * 3
    for i in eachindex(enrg)
        push!(result, round.([dim[i], enrg[i], 1, 0, PR], digits = 6))
    end
end

sort!(result, by = st -> real(st[1]))
for (lz, P) in ((0, 1), (0, -1), (1, 0))
    display(permutedims(hcat(
        filter(st -> st[3] ≊ lz && st[4] ≊ P, result)...
    )))
end