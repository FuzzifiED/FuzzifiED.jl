# This example calculates the spectrum of magnetic line defect in 3d Ising model
# in lz = 0, P = ±1 and lz = 1 sectors, calibrated by bulk T.
# This example reproduces Table I in Nat. Commun. 15, 3659 (2024)
# On my table computer, this calculation takes 0.962 s

using FuzzifiED
const σ1 = [  1  0 ;  0  0 ]
const σ2 = [  0  0 ;  0  1 ]
const σx = [  0  1 ;  1  0 ]
FuzzifiED.ElementType = Float64
≈(x, y) = abs(x - y) < eps(Float32)

let 

nm = 12
no = nm * 2
tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, 2, 2 .* [4.75, 1.], σ1, σ2) - 
    3.16 * GetPolTerms(nm, 2, σx) )
    
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
hmt = Operator(bs, tms_hmt)
hmt_mat = OpMat(hmt)
enrg, st = GetEigensystem(hmt_mat, 6)
enrg_0 = enrg[1]
enrg_T = enrg[3]

spec = []

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
enrg_d = 0
for P in (1, -1), R in (1, -1)
    bs = Basis(cfs, [P, R], qnf)
    hmt = Operator(bs, tms_hmt)
    hmt_mat = OpMat(hmt)
    enrg, st = GetEigensystem(hmt_mat, 10)
    if (P == 1 && R == 1) enrg_d = enrg[1] end
    dim = (enrg .- enrg_d) ./ (enrg_T - enrg_0) * 3
    for i in eachindex(enrg)
        push!(spec, round.([dim[i], enrg[i], 0, P, R] .+ √eps(Float64), digits = 6))
    end
end

cfs = Confs(no, [nm, 2, 2, 0], qnd)
for PR in (1, -1)
    bs = Basis(cfs, [PR], qnf1)
    hmt = Operator(bs, tms_hmt)
    hmt_mat = OpMat(hmt)
    enrg, st = GetEigensystem(hmt_mat, 10)
    dim = (enrg .- enrg_d) ./ (enrg_T - enrg_0) * 3
    for i in eachindex(enrg)
        push!(spec, round.([dim[i], enrg[i], 1, 0, PR] .+ √eps(Float64), digits = 6))
    end
end

sort!(spec, by = st -> real(st[1]))
display(permutedims(hcat(spec...)))

end