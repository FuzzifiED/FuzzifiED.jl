# This example calculates the spectrum of super-Ising SCFT. 
# This example reproduces Figure 12--14 in arXiv : 2509.08038
# On my portable computer, this calculation takes 67.520 s

using FuzzifiED
using FuzzifiED.Fuzzifino
FuzzifiED.ElementType = Float64
≈(x, y) = abs(x - y) < √eps(Float64)

nmf = 9
FuzzifiED.ObsNormRadSq = nmf
nof = 2 * nmf
nmb = nmf - 1
nob = nmb
qnd = [
    GetNeSQNDiag(nof, nob),
    GetBosonLz2SQNDiag(nof, nmb, 1) + SQNDiag(GetLz2QNDiag(nmf, 2), nob)
]
tms_l2 = GetL2STerms(nmf, 2, nmb, 1) 
cfs = Dict{Int64, SConfs}()
for lz = 0 : 1 
    cfs[lz] = SConfs(nof, nob, nmf, [nmf, lz], qnd)
end 

den_e = StoreComps(GetFerDensitySObs(nmf, 2) + GetBosDensitySObs(nmb, 1))
tms_int_e = SimplifyTerms(GetIntegral(den_e * den_e))

obs_bf = StoreComps(GetFermionSObs(nmf, 2, 1)' * GetBosonSObs(nmb, 1, 1))
tms_hop = SimplifyTerms(GetIntegral(obs_bf * DPlus(obs_bf))) ; 

den_σ = StoreComps(GetFerDensitySObs(nmf, 2, [0 1 ; 1 0]))
tms_int_σ = SimplifyTerms(GetIntegral(den_σ * Laplacian(den_σ))) 

den_χχ = StoreComps(GetBosDensitySObs(nmb, 1))
tms_int_yuk = SimplifyTerms(GetIntegral(den_σ * den_χχ + den_χχ * den_σ) / 2) 

tms_pol_σ = SimplifyTerms(GetIntegral(den_σ))
tms_pol_χχ= SimplifyTerms(GetIntegral(GetBosDensitySObs(nmb, 1)))
tms_pol_ϵ = SimplifyTerms(GetIntegral(GetFerDensitySObs(nmf, 2, [0 0 ; 0 1]))) 

filter!(tm -> length(tm.cstr) > 4, tms_int_e)
filter!(tm -> length(tm.cstr) > 4, tms_hop)
filter!(tm -> length(tm.cstr) > 4, tms_int_σ)
filter!(tm -> length(tm.cstr) > 4, tms_int_yuk)

GetHmtTerms(t, U, g, hx, hz, μ) = SimplifyTerms(
    tms_int_e - t * (tms_hop + tms_hop')
    + U * tms_int_σ + g * tms_int_yuk
    - hx * tms_pol_σ - hz * tms_pol_ϵ - μ * tms_pol_χχ
) 

tms_hmt = GetHmtTerms(1.5, 0.25, 1.0, 0.0698080753241176, 0.0773934248387392, 0.07224789112627583)
result = []
nst = Dict(0 => 20, 1 => 20)
for lz = 0 : 1
    bs = SBasis(cfs[lz])
    hmt = SOperator(bs, tms_hmt)
    hmt_mat = OpMat(hmt)
    enrg, st = GetEigensystem(hmt_mat, nst[lz])

    l2 = SOperator(bs, tms_l2)
    l2_mat = OpMat(l2)
    l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)] ; 

    nf = SOperator(bs, STerms(GetPolTerms(nof, 2, [0 0 ; 0 1])) / nof)
    nf_mat = OpMat(nf)
    nf_val = [ st[:, i]' * nf_mat * st[:, i] for i in eachindex(enrg)]

    for i in eachindex(enrg)
        push!(result, round.([enrg[i], l2_val[i], nf_val[i]], digits = 6))
    end
end
sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_σ = filter(st -> st[2] ≈ 0, result)[2][1]
spec = [ [ 0.58444 * 1.0061416719249328 * (st[1] - enrg_0) / (enrg_σ - enrg_0) ; st] for st in result ] 
display(permutedims(stack(spec)))