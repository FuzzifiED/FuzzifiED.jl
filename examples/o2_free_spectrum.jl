# This example calculates the spectrum of O(2) free scalar CFT.
# On my table computer, this calculation takes 6.836 s

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

mat_0 = diagm(ComplexF64[0, 0, 1])
mat_ω = [
    [0 0  1 ; 0  0 1 ; -1 -1 0] * im / √2, 
    [0 0 -1 ; 0  0 1 ; -1  1 0] / √2,
    [1 0  0 ; 0 -1 0 ;  0  0 0]
]

FuzzifiED.ObsNormRadSq = nm 
obs_ne = GetDensityObs(nm, nf)
obs_ω = GetDensityObs.(nm, nf, mat_ω)
tms_l2 = GetL2Terms(nm, nf)

cfs = Dict{Int64, Confs}()
for s = 0 : 3 
    cfs[s] = Confs(no, [nm, 0, s], qnd)
end
tms_hmt = SimplifyTerms(
    GetIntegral(obs_ne' * obs_ne)
    + 0.148 * GetIntegral(obs_ω' * Laplacian.(obs_ω))
    - 0.070 / nm  * GetPolTerms(nm, nf, mat_0)
)

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
enrg_ϕ = filter(st -> st[2] ≈ 0 && st[3] ≈ 1, result)[1][1]
spec = [ round.([ 0.5 * (st[1] - enrg_0) / (enrg_ϕ - enrg_0) ; st] .+ √eps(Float64), digits = 6) for st in result ]
display(permutedims(hcat(spec...)))
