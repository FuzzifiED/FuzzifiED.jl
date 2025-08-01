# This example calculates the spectrum of free real scalar.
# This example reproduces Figure 2 in arXiv : 2506.14904. 
# On my portable computer, this calculation takes 5.288 s

using FuzzifiED
const σx = [  0  1 ;  1  0 ]
const σy = [  0 im ;-im  0 ]
const σz = [  1  0 ;  0 -1 ]
FuzzifiED.ElementType = Float64
≈(x, y) = abs(x - y) < eps(Float32)

nm = 12
nf = 2
no = nm * nf
qnd = [ 
    GetNeQNDiag(no), 
    GetLz2QNDiag(nm, nf),
    GetFlavQNDiag(nm, nf, [0, 1], 1, 2)
]
qnf = [ 
    GetParityQNOffd(nm, nf, [2, 1], [-1, 1]), 
    GetRotyQNOffd(nm, nf) 
]

FuzzifiED.ObsNormRadSq = nm
n0 = GetDensityObs(nm, nf)
ny = GetDensityObs(nm, nf, σy)
nz = GetDensityObs(nm, nf, σz)
tms_hmt = SimplifyTerms(
    GetIntegral(n0 * n0)
    + 0.18 * GetIntegral(ny * Laplacian(ny))
    - 0.42 / 4π / nm * GetIntegral(nz)
)
tms_l2 = GetL2Terms(nm, 2)

cfs = Dict{Int64, Confs}()
cfs[ 1] = Confs(2 * nm, [nm, 0, 0], qnd)
cfs[-1] = Confs(2 * nm, [nm, 0, 1], qnd)

result = []
for P in [1, -1], Z in [1, -1], R in [1, -1]
    bs = Basis(cfs[Z], [P, R], qnf)
    hmt = Operator(bs, tms_hmt)
    hmt_mat = OpMat(hmt)
    enrg, st = GetEigensystem(hmt_mat, 20)

    l2 = Operator(bs, tms_l2)
    l2_mat = OpMat(l2)
    l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)]

    for i in eachindex(enrg)
        push!(result, round.([enrg[i], l2_val[i], P, Z], digits = 6))
    end
end

sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_ϕ = result[2][1]
result_dim = [ [ 0.5 * (st[1] - enrg_0) / (enrg_ϕ - enrg_0) ; st] for st in result ]
display(permutedims(hcat(result_dim...)))