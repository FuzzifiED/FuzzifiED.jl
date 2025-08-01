# This example calculates the spectrum of the U(1)_2 coupled to 
# a complex scalar Chern-Simons matter CFT on the fuzzy sphere. 
# This example reproduces Figs. 2d and 3 in arXiv : 2507.19580.
# On my portable computer, this calculation takes 2.405 s

using FuzzifiED
using FuzzifiED.Fuzzifino
FuzzifiED.ElementType = Float64 
≈(x, y) = abs(x - y) < √(eps(Float64))

nmf = 6
nof = 2 * nmf 
s = (nmf - 1) / 2
nmb = nmf * 2 - 1 
nob = nmb 

qnd = [
    SQNDiag(GetNeQNDiag(nof), nob) + 2 * GetBosonNeSQNDiag(nof, nob), 
    SQNDiag(GetLz2QNDiag(nmf, 2), nob) + GetBosonLz2SQNDiag(nof, nmb, 1), 
    SQNDiag(GetFlavQNDiag(nmf, 2, [1, -1]), nob)
]

qnf = [
    SQNOffd(GetFlavPermQNOffd(nmf, 2, [2, 1], [1, -1]), nob)
    SQNOffd(GetRotyQNOffd(nmf, 2), nob) * GetBosonRotySQNOffd(nof, nmb, 1)
]

cfs = SConfs(nof, nob, nmf, [2 * nmf, 0, 0], qnd) 

tms_l2 = GetL2STerms(nmf, 2, nmb, 1)
tms_c2 = STerms(GetC2Terms(nmf, 2, :SU)) 

FuzzifiED.ObsNormRadSq = nmf
tms_hop = GetIntegral(GetBosonSObs(nmb, 1, 1)' * GetFermionSObs(nmf, 2, 1) * GetFermionSObs(nmf, 2, 2))
den_e = StoreComps(GetFerDensitySObs(nmf, 2) + 2 * GetBosDensitySObs(nmb, 1))
tms_int_e = GetIntegral(den_e * den_e)

tms_hmt = SimplifyTerms(
    tms_int_e
    - 0.5 * (tms_hop + tms_hop') 
    + 0.312 * STerms(GetPolTerms(nmf, 2))
)
result = []
for Z in [1, -1], R in [-1, 1]
    bs = SBasis(cfs, [R, Z], qnf)
    hmt = SOperator(bs, tms_hmt)
    hmt_mat = OpMat(hmt)
    enrg, st = GetEigensystem(hmt_mat, 20)

    l2 = SOperator(bs, tms_l2)
    l2_mat = OpMat(l2)
    l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)]

    c2 = SOperator(bs, tms_c2)
    c2_mat = OpMat(c2)
    c2_val = [ st[:, i]' * c2_mat * st[:, i] for i in eachindex(enrg)]
    @show c2_val

    for i in eachindex(enrg)
        (st[:, i]' * st[:, i] ≈ 1) || continue
        push!(result, round.([enrg[i], l2_val[i], c2_val[i]] .+ √eps(Float64), digits = 7))
    end
end

sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_1 = (filter(st -> st[2] ≈ 6 && st[3] ≈ 0, result)[1][1] - enrg_0) / 3
result_dim = [ [ (st[1] - enrg_0) / enrg_1 ; st] for st in result if st[2] < 12 ]
display(permutedims(hcat(result_dim...)))

# As a minimal example, in this code, we calibrate by setting Δ_T=3.
# To calibrate by optimising conformality, replace Lines 69 with
# enrg_cal = [
#     filter(st -> st[2] ≈ 2 && st[3] ≈ 0, result)[1][1] - filter(st -> st[2] ≈ 0 && st[3] ≈ 0, result)[2][1], # ∂S - S
#     filter(st -> st[2] ≈ 2 && st[3] ≈ 2, result)[1][1] - enrg_0, # J
#     filter(st -> st[2] ≈ 2 && st[3] ≈ 2, result)[2][1] - enrg_0, # ϵ∂J
#     filter(st -> st[2] ≈ 6 && st[3] ≈ 2, result)[1][1] - enrg_0, # ∂J
#     filter(st -> st[2] ≈ 6 && st[3] ≈ 0, result)[1][1] - enrg_0 # T
# ]
# dim_cal = Float64[1, 2, 3, 3, 3]
# enrg_1 = (enrg_cal' * dim_cal) / (dim_cal' * dim_cal)
