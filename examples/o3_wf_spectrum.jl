# This example calculates the spectrum of O(3) Wilson-Fisher CFT.
# On my table computer, this calculation takes 6.576 s

using FuzzifiED
using LinearAlgebra
FuzzifiED.ElementType = Float64
≈(x, y) = abs(x - y) < √eps(Float64)

nm = 8
nf = 4
no = nm * nf 

# The lines below are the same for all O(N) Wilson-Fisher and free scalars

qnd = [
    GetNeQNDiag(no) ;
    GetLz2QNDiag(nm, nf) ;
    [ GetFlavQNDiag(nm, nf, Dict([i => 1, i + 1 => -1])) for i = 1 : 2 : nf - 2 ] ; 
]
if (iseven(nf)) 
    push!(qnd, GetFlavQNDiag(nm, nf, Dict([f => 1 for f = 1 : nf - 1]), 0, 2))
end 
qnf = [
    GetRotyQNOffd(nm, nf) ; 
    [GetFlavPermQNOffd(nm, nf, Dict([i => i+1, i+1 => i])) for i = 1 : 2 : nf - 2] ; 
    [GetFlavPermQNOffd(nm, nf, Dict([i => i+2, i+1 => i+3, i+2 => i, i+3 => i+1])) for i = 1 : 2 : nf - 4 ]
] 

function mat_one(i, j)
    mat = zeros(nf, nf)
    mat[i, j] = 1.0
    return mat
end
mat_U = Matrix{ComplexF64}(I, nf, nf)
for i = 1 : 2 : nf - 2
    mat_U[i : i + 1, i : i + 1] = [1/√2 1/√2 ; im/√2 -im/√2]
end
mat_0 = mat_one(nf, nf)
mat_A = [mat_U' * (mat_one(i, j) - mat_one(j, i)) * mat_U for i = 1 : nf - 2 for j = i + 1 : nf - 1 ]
mat_Vp = [mat_U' * (mat_one(i, nf) + mat_one(nf, i)) * mat_U for i = 1 : nf - 1] 
mat_Vm = [mat_U' * (mat_one(i, nf) - mat_one(nf, i)) * mat_U for i = 1 : nf - 1] 

FuzzifiED.ObsNormRadSq = nm 
obs_n0 = GetDensityObs(nm, nf, mat_0)
obs_ne = GetDensityObs(nm, nf)
obs_nVp = GetDensityObs.(nm, nf, mat_Vp)
obs_nVm = GetDensityObs.(nm, nf, mat_Vm)
obs_nA = GetDensityObs.(nm, nf, mat_A)

tms_l2 = GetL2Terms(nm, nf)
tms_c2 = 4 * GetC2Terms(nm, nf, mat_A) 

# The lines above are the same for all O(N) Wilson-Fisher and free scalars

cfs = Dict{Int64, Confs}()
for Z = 0 : 1
    cfs[Z] = Confs(no, [nm, 0, 0, Z], qnd)
end 
tms_hmt = SimplifyTerms(
    GetIntegral(obs_ne' * obs_ne)
    - 0.264 * GetIntegral(obs_nVp' * Laplacian.(obs_nVp))
    - 0.096 * GetPolTerms(nm, nf, mat_0)
)

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
        push!(result, [enrg[i], l2_val[i], c2_val[i], (-1) ^ Z, X])
    end
end

sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_T = filter(st -> st[2] ≈ 6 && st[3] ≈ 0 && st[4] ≈ 1, result)[1][1]
spec = [ round.([ 3 * (st[1] - enrg_0) / (enrg_T - enrg_0) ; st] .+ √eps(Float64), digits = 6) for st in result ]
display(permutedims(hcat(spec...)))
