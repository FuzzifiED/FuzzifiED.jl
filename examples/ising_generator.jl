# This example examines the quality of conformal symmetry at nm = 12
# by examining the matrix elements of conformal generators P^z + K^z
# and compare the states (P^z + K^z)|Φ⟩ with the CFT expectations. 
# This example reproduces Figure 7 in SciPost Phys. 18, 086 (2025).
# On my portable computer, this calculation takes 2.490 s

using FuzzifiED
using LinearAlgebra
const σ1 = [  1  0 ;  0  0 ]
const σp = [  0  1 ;  0  0 ]
const σ2 = [  0  0 ;  0  1 ]
const σx = [  0  1 ;  1  0 ]
const Δϵ = 1.412625
const Δσ = 0.5181489
FuzzifiED.ElementType = Float64
≈(x, y) = abs(x - y) < eps(Float32)

compare_st(st0, st1) = abs(st0' * st1) ^ 2 / ((st0' * st0) * (st1' * st1)) ;

nm = 12
qnd = [ 
    GetNeQNDiag(2 * nm), 
    GetLz2QNDiag(nm, 2) ]
qnf = [ 
    GetParityQNOffd(nm, 2, [2, 1], [-1, 1]), 
    GetFlavPermQNOffd(nm, 2, [2, 1]), 
    GetRotyQNOffd(nm, 2) ]

tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, 2, 2 .* [4.75, 1.], σ1, σ2) - 
    3.16 * GetPolTerms(nm, 2, σx) )
tms_l2 = GetL2Terms(nm, 2)

cfs = Confs(2 * nm, [nm, 0], qnd)

result = []
bss = Dict{Vector{Int64}, Basis}()
l2_mats = Dict{Vector{Int64}, OpMat}()
for P in [1], Z in [1, -1], R in [1, -1]
    bs = Basis(cfs, [P, Z, R], qnf)
    hmt = Operator(bs, tms_hmt)
    hmt_mat = OpMat(hmt)
    enrg, st = GetEigensystem(hmt_mat, 5)

    l2 = Operator(bs, tms_l2)
    l2_mat = OpMat(l2)
    l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)]

    bss[[P, Z, R]] = bs
    l2_mats[[P, Z, R]] = l2_mat

    for i in eachindex(enrg)
        push!(result, ([enrg[i], l2_val[i], P, Z, st[:, i]]))
    end
end

sort!(result, by = st -> st[1])
st_0  = result[1][end]
st_ϵ  = filter(st -> st[2] ≈ 0 && st[3] ≈ 1 && st[4] ≈ 1, result)[2][end]
st_∂ϵ = filter(st -> st[2] ≈ 2 && st[3] ≈ 1 && st[4] ≈ 1, result)[1][end]
st_∂∂ϵ= filter(st -> st[2] ≈ 6 && st[3] ≈ 1 && st[4] ≈ 1, result)[2][end]
st_□ϵ = filter(st -> st[2] ≈ 0 && st[3] ≈ 1 && st[4] ≈ 1, result)[3][end]
st_σ  = filter(st -> st[2] ≈ 0 && st[3] ≈ 1 && st[4] ≈-1, result)[1][end]
st_∂σ = filter(st -> st[2] ≈ 2 && st[3] ≈ 1 && st[4] ≈-1, result)[1][end]
st_∂∂σ= filter(st -> st[2] ≈ 6 && st[3] ≈ 1 && st[4] ≈-1, result)[1][end]
st_□σ = filter(st -> st[2] ≈ 0 && st[3] ≈ 1 && st[4] ≈-1, result)[2][end]

pair_pm = StoreComps(GetPairingMod(nm, 2, σp))
den_x = GetDensityObs(nm, 2, σx)

tms_cand = [ 
    [GetComponent(den_x, 1, 0)] ; 
    [GetComponent(FilterL2(pair_pm, l)' * FilterL2(pair_pm, l), 1, 0) for l in nm - 2 : nm - 1]
]
opst = [ Operator(bss[[1,1,1]], bss[[1,1,-1]], tms ; red_q = 1) * st_0 for tms in tms_cand]

mat = [sti' * stj for sti in opst, stj in opst]
eigval, eigvec = eigen(mat);
@show real(eigval)

tms_pk = SimplifyTerms(eigvec[:, 1]' * tms_cand)
norm_pk = st_∂ϵ' * Operator(bss[[1,1,1]], bss[[1,1,-1]], tms_pk ; red_q = 1) * st_ϵ / √(2Δϵ)
tms_pk /= norm_pk

pk_mats = Dict{Vector{Int64}, OpMat}()
pk_mats[[1, 1, 1]] = OpMat(Operator(bss[[1, 1, 1]], bss[[1, 1,-1]], tms_pk ; red_q = 1))
pk_mats[[1,-1, 1]] = OpMat(Operator(bss[[1,-1, 1]], bss[[1,-1,-1]], tms_pk ; red_q = 1))
pk_mats[[1, 1,-1]] = OpMat(Operator(bss[[1, 1,-1]], bss[[1, 1, 1]], tms_pk ; red_q = 1))
pk_mats[[1,-1,-1]] = OpMat(Operator(bss[[1,-1,-1]], bss[[1,-1, 1]], tms_pk ; red_q = 1))
st_pk∂ϵ = pk_mats[[1, 1,-1]] * st_∂ϵ
st_pk∂σ = pk_mats[[1,-1,-1]] * st_∂σ
st_pkϵ  = pk_mats[[1, 1, 1]] * st_ϵ
st_pkσ  = pk_mats[[1,-1, 1]] * st_σ
st_pk∂ϵ_2 = l2_mats[[1, 1, 1]] * st_pk∂ϵ ./ 6
st_pk∂ϵ_0 = st_pk∂ϵ .- st_pk∂ϵ_2
st_pk∂σ_2 = l2_mats[[1,-1, 1]] * st_pk∂σ ./ 6
st_pk∂σ_0 = st_pk∂σ .- st_pk∂σ_2

ovl_fzs_ϵ_∂ϵ   = abs(st_pk∂ϵ' * st_ϵ  )
ovl_fzs_∂ϵ_∂∂ϵ = abs(st_pk∂ϵ' * st_∂∂ϵ)
ovl_fzs_∂ϵ_□ϵ  = abs(st_pk∂ϵ' * st_□ϵ )
ovl_fzs_σ_∂σ   = abs(st_pk∂σ' * st_σ  )
ovl_fzs_∂σ_∂∂σ = abs(st_pk∂σ' * st_∂∂σ)
ovl_fzs_∂σ_□σ  = abs(st_pk∂σ' * st_□σ )
ovl_cft_ϵ_∂ϵ   = √(2Δϵ)
ovl_cft_∂ϵ_∂∂ϵ = √(8(Δϵ+1)/3)
ovl_cft_∂ϵ_□ϵ  = √(2(2Δϵ-1)/3)
ovl_cft_σ_∂σ   = √(2Δσ)
ovl_cft_∂σ_∂∂σ = √(8(Δσ+1)/3)
ovl_cft_∂σ_□σ  = √(2(2Δσ-1)/3)
@show ovl_fzs_ϵ_∂ϵ  , ovl_cft_ϵ_∂ϵ   
@show ovl_fzs_∂ϵ_∂∂ϵ, ovl_cft_∂ϵ_∂∂ϵ 
@show ovl_fzs_∂ϵ_□ϵ , ovl_cft_∂ϵ_□ϵ  
@show ovl_fzs_σ_∂σ  , ovl_cft_σ_∂σ   
@show ovl_fzs_∂σ_∂∂σ, ovl_cft_∂σ_∂∂σ 
@show ovl_fzs_∂σ_□σ , ovl_cft_∂σ_□σ  

@show compare_st(st_∂σ , st_pkσ   )
@show compare_st(st_□σ , st_pk∂σ_0)
@show compare_st(st_σ  , st_pk∂σ_0)
@show compare_st(st_∂∂σ, st_pk∂σ_2)
@show compare_st(st_∂ϵ , st_pkϵ   )
@show compare_st(st_□ϵ , st_pk∂ϵ_0)
@show compare_st(st_ϵ  , st_pk∂ϵ_0)
@show compare_st(st_∂∂ϵ, st_pk∂ϵ_2)