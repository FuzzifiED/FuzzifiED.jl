# input nmf before calling this

using FuzzifiED
using FuzzifiED.Fuzzifino
FuzzifiED.SilentStd = true
FuzzifiED.ElementType = Float64 
≈(x, y) = abs(x - y) < √(eps(Float64))
σx = [  0  1 ;  1  0 ]
σy = [  0 im ;-im  0 ]
σz = [  1  0 ;  0 -1 ]

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

tms_lz = [STerm(m - s - 1, [1, (m - 1) * 2 + f, 0, (m - 1) * 2 + f]) for m = 1 : nmf for f = 1 : 2] +
    [ STerm(m - 2s - 1, [1, -m, 0, -m]) for m = 1 : nmb ]
tms_lp = [ STerm(sqrt((nmf - m) * m), [1, m * 2 + f, 0, (m - 1) * 2 + f]) for m = 1 : nmf - 1 for f = 1 : 2 ] +
    [ STerm(sqrt((nmb - m) * m), [1, -(m + 1), 0, -m]) for m = 1 : nmb - 1 ]
tms_lm = tms_lp' 
tms_l2 = SimplifyTerms(tms_lz * tms_lz - tms_lz + tms_lp * tms_lm) 

tms_c2 = STerms(GetC2Terms(nmf, 2, [σx, σy, σz])) 

FuzzifiED.ObsNormRadSq = nmf
tms_hop = GetIntegral(GetBosonSObs(nmb, 1, 1)' * GetFermionSObs(nmf, 2, 1) * GetFermionSObs(nmf, 2, 2))
den_f = StoreComps(GetFerDensitySObs(nmf, 2))
den_b = StoreComps(GetBosDensitySObs(nmb, 1))
tms_int_0 = GetIntegral(den_f * den_f + 4 * den_b * den_b)
tms_int_1 = GetIntegral(4 * den_f * den_b)

tms_hmt = SimplifyTerms(
    tms_int_0 
    + tms_int_1 
    - 0.5 * (tms_hop + tms_hop') 
    + 0.312 * STerms(GetPolTerms(nmf, 2))
)
result = []
for Z in [1, -1], R in [-1, 1]
    bs = SBasis(cfs, [R, Z], qnf)
    hmt = SOperator(bs, tms_hmt + √eps(Float64) * tms_l2)
    hmt_mat = OpMat(hmt)
    enrg, st = GetEigensystem(hmt_mat, 20)

    l2 = SOperator(bs, tms_l2)
    l2_mat = OpMat(l2)
    l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)]

    c2 = SOperator(bs, tms_c2)
    c2_mat = OpMat(c2)
    c2_val = [ st[:, i]' * c2_mat * st[:, i] for i in eachindex(enrg)]

    for i in eachindex(enrg)
        (st[:, i]' * st[:, i] ≈ 1) || continue
        push!(result, round.([enrg[i], l2_val[i], c2_val[i]] .+ √eps(Float64), digits = 7))
    end
end

sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_T = filter(st -> st[2] ≈ 6 && st[3] ≈ 0, result)[1][1]
result_dim = [ [ 3 * (st[1] - enrg_0) / (enrg_T - enrg_0) ; st] for st in result if st[2] < 12 ]
display(permutedims(hcat(result_dim...)))
