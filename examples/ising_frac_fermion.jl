# This example calculates the spectrum of 3d Ising model on fuzzy sphere
# for fermions at fractional filling ν = 1/3.
# This example reproduces Figure 10 in Phys. Rev. X 15, 031007 (2025).
# On my portable computer, this calculation takes 33.527 s.
# We acknowlege Cristian Voinea for his help in reproducing the results. 

using FuzzifiED
const σ1 = [ 1 0 ; 0 0 ]
const σ2 = [ 0 0 ; 0 1 ]
const σ0 = [ 1 0 ; 0 1 ]
const σx = [ 0 1 ; 1 0 ]
FuzzifiED.ElementType = Float64
≈(x, y) = abs(x - y) < eps(Float32)

ne = 7
nm = 3 * ne - 2
qnd = [ 
    GetNeQNDiag(2 * nm), 
    GetLz2QNDiag(nm, 2) ]
qnf = [ 
    GetFlavPermQNOffd(nm, 2, [2, 1]), 
    GetRotyQNOffd(nm, 2) ]

tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, 2, [1.0, 1.0])
    + GetDenIntTerms(nm, 2, 2 .* [0.0, 0.0, 0.49, 0.09],  σ1, σ2)
    - 0.183 * GetPolTerms(nm, 2, σx) 
)
tms_l2 = GetL2Terms(nm, 2)

cfs = Confs(2 * nm, [ne, 0], qnd)

result = []
for Z in [1, -1], R in [1, -1]
    bs = Basis(cfs, [Z, R], qnf)
    hmt = Operator(bs, tms_hmt)
    hmt_mat = OpMat(hmt)
    enrg, st = GetEigensystem(hmt_mat, 20)

    l2 = Operator(bs, tms_l2)
    l2_mat = OpMat(l2)
    l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)]

    for i in eachindex(enrg)
        push!(result, [enrg[i], l2_val[i], Z])
    end
end

sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_T = filter(st -> st[2] ≈ 6 && st[3] ≈ 1, result)[1][1]
spec = [ round.([ 3 * (st[1] - enrg_0) / (enrg_T - enrg_0) ; st] .+ √eps(Float64), digits = 6) for st in result ]
display(permutedims(hcat(spec...)))