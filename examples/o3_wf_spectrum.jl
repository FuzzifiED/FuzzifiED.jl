# This example calculates the spectrum of O(3) Wilson-Fisher CFT.
# This example takes the model from 2510.09755 
# and reproduces partly Tables I and II, Figures 1 and 2.
# On my table computer, this calculation takes 6.184 s

using FuzzifiED
using LinearAlgebra
FuzzifiED.ElementType = Float64
≈(x, y) = abs(x - y) < √eps(Float64)

mat_0 = diagm([0, 0, 0, 1.0])
mat_A = [
    [ 0  1 0 0 ; 1 0  1 0 ; 0 1  0 0 ; 0 0 0 0 ] / √2, # Ax
    [ 0 -1 0 0 ; 1 0 -1 0 ; 0 1  0 0 ; 0 0 0 0 ] / √2 * im, # Ay
    [ 1  0 0 0 ; 0 0  0 0 ; 0 0 -1 0 ; 0 0 0 0 ] # Az
]
mat_V = [
    [ 0 0 0  1 ; 0 0 0 0 ; 0 0 0 -1 ; 1 0 -1 0 ] / √2, # Vx
    [ 0 0 0 -1 ; 0 0 0 0 ; 0 0 0 -1 ; 1 0  1 0 ] / √2 * im, # Vy
    [ 0 0 0  0 ; 0 0 0 1 ; 0 0 0  0 ; 0 1  0 0 ] # Vz
]

nm = 8
nf = 4
no = nm * nf 

qnd = [
    GetNeQNDiag(no),
    GetLz2QNDiag(nm, nf),
    GetFlavQNDiag(nm, nf, Dict([1 => 1, 3 => -1])), 
    GetFlavQNDiag(nm, nf, Dict([f => 1 for f = 1 : nf - 1]), 0, 2)
]
qnf = [
    GetRotyQNOffd(nm, nf),
    GetFlavPermQNOffd(nm, nf, Dict([1 => 3, 3 => 1]))
] ;

cfs = Dict{Int64, Confs}()
for Z = 0 : 1
    cfs[Z] = Confs(no, [nm, 0, 0, Z], qnd)
end 

h = 14.992 + 6.59 * nm ^ (-2.16 / 2)
tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, nf, [6.5, 1.0])
    - 1.4 * 2 * GetDenIntTerms(nm, nf, [6.5, 1.0], mat_V)
    - 2 * h * GetPolTerms(nm, nf, mat_0)
)
tms_l2 = GetL2Terms(nm, nf)
tms_c2 = 4 * GetC2Terms(nm, nf, mat_A) ;

result = []
for Z in [0, 1], X in [1, -1], R in [1, -1]
    bs = Basis(cfs[Z], [R, X], qnf)
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
        push!(result, [enrg[i], l2_val[i], c2_val[i], X])
    end
end

sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_ϕ = filter(st -> st[2] ≈ 0 && st[3] ≈ 2 && st[4] ≈ 1, result)[1][1]
spec = [ round.([ 0.518936 * (st[1] - enrg_0) / (enrg_ϕ - enrg_0) ; st] .+ √eps(Float64), digits = 6) for st in result ]
display(permutedims(hcat(spec...)))
