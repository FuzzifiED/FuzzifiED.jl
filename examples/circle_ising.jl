# This example calculates the spectrum of 2d Ising CFT on a fuzzy thin torus. 
# This example reproduces Figure 4 and Tables I--III in Phys. Rev. B 111, 085113 (2025)
# On my table computer, this calculation takes 3.498 s
# We acknowlege Wei Zhu for his help in reproducing the results. 

using FuzzifiED
using FuzzifiED.FuzzyManifolds
FuzzifiED.ElementType = Float64
const σ1 = [  1  0 ;  0  0 ]
const σ2 = [  0  0 ;  0  1 ]
const σx = [  0  1 ;  1  0 ]

nf = 2
nm = 12
no = nf * nm
qnd = [
    GetNeQNDiag(no),
    GetTorusLz2QNDiag(nm, nf)
]
qnf = [ 
    GetFlavPermQNOffd(nm, 2, [2, 1]),
    GetTorusTranslQNOffd(nm, 2)
]

cfs = Dict{Int64, Confs}()
for lz = 0 : 3
    cfs[lz] = Confs(no, [nm, 2 * lz], qnd)
end

ps_pot = [4, 1] * 2
h = 0.9216
lx = 3.25
tms_hmt = SimplifyTerms(
    GetTorusDenIntTerms(nm, nf, lx, ps_pot, σ1, σ2) - 
    h * GetPolTerms(nm, 2, σx) 
)

result = []
for lz = 0 : 3, Z in [1, -1]
    bs = Basis(cfs[lz], [Z, -(-1)^nm], qnf)
    hmt = Operator(bs, tms_hmt)
    hmt_mat = OpMat(hmt)
    enrg, st = GetEigensystem(hmt_mat, 20)

    for i in eachindex(enrg)
        push!(result, [enrg[i], lz, Z])
    end
end

sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_T = filter(st -> st[2] ≈ 2 && st[3] ≈ 1, result)[1][1]
spec = [ round.(real([ (st[1] - enrg_0) / (enrg_T - enrg_0) * 2 ; st]) .+ eps(Float32), digits = 6) for st in result ]
display(permutedims(hcat(spec...)))
