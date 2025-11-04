# This example calculates the spectrum of the defect creation and changing operators 
# of the magnetic line defect in 3d Ising model.
# This example reproduces Table 2 and Figure 5 in SciPost Phys. 17, 021 (2024).
# On my table computer, this calculation takes 1.733 s

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
enrg, st = GetEigensystem(hmt_mat, 5)
enrg_0 = enrg[1]
enrg_T = enrg[3]

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
cfs = Confs(no, [nm, 0, 2, 0], qnd)
bs = Basis(cfs, [1, 1], qnf)
hmt = Operator(bs, tms_hmt)
hmt_mat = OpMat(hmt)
enrg, st = GetEigensystem(hmt_mat, 3)
enrg_d = enrg[1]

qnd = [
    GetNeQNDiag(no),
    GetLz2QNDiag(nm, 2), 
    GetPinOrbQNDiag(no, [1]), 
    GetPinOrbQNDiag(no, [2])
]
qnf = [
    GetParityQNOffd(nm, 2, [2, 1], [-1, 1])
]
cfs = Confs(no, [nm, 0, 1, 0], qnd)
spec_p0 = []
for P in (1, -1)
    bs = Basis(cfs, [P], qnf)
    hmt = Operator(bs, tms_hmt)
    hmt_mat = OpMat(hmt)
    enrg, st = GetEigensystem(hmt_mat, 10)
    dim = (enrg .- (enrg_d + enrg_0) / 2) ./ (enrg_T - enrg_0) * 3
    for i in eachindex(enrg)
        push!(spec_p0, round.([dim[i], enrg[i], 0, P] .+ √eps(Float64), digits = 6))
    end
end
sort!(spec_p0, by = st -> real(st[1]))

qnd = [
    GetNeQNDiag(no),
    GetLz2QNDiag(nm, 2), 
    GetPinOrbQNDiag(no, [1, no]), 
    GetPinOrbQNDiag(no, [2, no - 1])
]
qnf = [
    GetParityQNOffd(nm, 2, [2,1],[-1,1]),    
    GetRotyQNOffd(nm, 2) * GetFlavPermQNOffd(nm, 2, [2,1])
]
cfs = Confs(no, [nm, 0, 2, 0], qnd)
spec_pm = []
for P in (1, -1), RZ in (1, -1)
    bs = Basis(cfs, [P, RZ], qnf)
    hmt = Operator(bs, tms_hmt)
    hmt_mat = OpMat(hmt)
    enrg, st = GetEigensystem(hmt_mat, 10)
    dim = (enrg .- enrg_d) ./ (enrg_T - enrg_0) * 3
    for i in eachindex(enrg)
        push!(spec_pm, round.([dim[i], enrg[i], 0, P, RZ] .+ √eps(Float64), digits = 6))
    end
end
sort!(spec_pm, by = st -> real(st[1]))

display(permutedims(hcat(spec_p0...)))
display(permutedims(hcat(spec_pm...)))

end