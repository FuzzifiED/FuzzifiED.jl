# This example defines a cost function as the square sum of the deviations of 
# ∂^nσ, ∂^nϵ and T to evaluate the conformal symmetry for Ising model 
# and minimises this cost function to find the best parameter.
# On my table computer, this calculation takes 5.061 s
# We acknowlege Wenhan Guo for his help in optimising this example. 

using FuzzifiED
using Optim
FuzzifiED.SilentStd = true
FuzzifiED.ElementType = Float64
const σ1 = [  1  0 ;  0  0 ]
const σ2 = [  0  0 ;  0  1 ]
const σx = [  0  1 ;  1  0 ]
≈(x, y) = abs(x - y) < eps(Float32)

nm = 10
qnd = [ 
    GetNeQNDiag(2 * nm), 
    GetLz2QNDiag(nm, 2) ]
qnf = [ 
    GetParityQNOffd(nm, 2, [2, 1], [-1, 1]), 
    GetFlavPermQNOffd(nm, 2, [2, 1]), 
    GetRotyQNOffd(nm, 2) ]

tms_l2 = GetL2Terms(nm, 2)

cfs = Confs(2 * nm, [nm, 0], qnd)

secf_lst = [
    (1, 1, 1, 5),
    (1, 1,-1, 3),
    (1,-1, 1, 5),
    (1,-1,-1, 3)
]

function cost(ps_pot)
    tms_hmt = SimplifyTerms(
        GetDenIntTerms(nm, 2, ps_pot, σ1, σ2) - 
        GetPolTerms(nm, 2, σx) 
    )
    result = []
    for (P, Z, R, nst) in secf_lst
        bs = Basis(cfs, [P, Z, R], qnf)
        hmt = Operator(bs, tms_hmt)
        hmt_mat = OpMat(hmt)
        enrg, st = GetEigensystem(hmt_mat, nst)

        l2 = Operator(bs, tms_l2)
        l2_mat = OpMat(l2)
        l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)]

        for i in eachindex(enrg)
            push!(result, [enrg[i], l2_val[i], P, Z])
        end
    end

    sort!(result, by = st -> real(st[1]))
    try
        enrg_0 = result[1][1]
        enrg_T = filter(st -> st[2] ≈ 6 && st[3] ≈ 1 && st[4] ≈ 1, result)[1][1]
        spec = [ [ 3 * (st[1] - enrg_0) / (enrg_T - enrg_0) ; st] for st in result ]
        Δσ   = filter(st -> st[3] ≈ 0 && st[5] ≈ -1, spec)[1][1]
        Δ∂σ  = filter(st -> st[3] ≈ 2 && st[5] ≈ -1, spec)[1][1]
        Δ∂∂σ = filter(st -> st[3] ≈ 6 && st[5] ≈ -1, spec)[1][1]
        Δ□σ  = filter(st -> st[3] ≈ 0 && st[5] ≈ -1, spec)[2][1]
        Δϵ   = filter(st -> st[3] ≈ 0 && st[5] ≈  1, spec)[2][1]
        Δ∂ϵ  = filter(st -> st[3] ≈ 2 && st[5] ≈  1, spec)[1][1]

        dim_1 = [Δ∂σ - Δσ, Δ∂∂σ - Δσ, Δ□σ - Δσ, Δ∂ϵ - Δϵ, 3.0]
        dim_0 = Float64[1, 2, 2, 1, 3]
        coeff_2 = sum(abs.(dim_1) .^ 2)
        coeff_1 = sum(conj(dim_1) .* dim_0)
        coeff_0 = sum(dim_0 .^ 2)
        renorm = coeff_1 / coeff_2 
        cost = coeff_0 - abs(coeff_1) ^ 2 / coeff_2
        lin_cost = sqrt(cost / 5)
        @show ps_pot
        @show renorm, lin_cost
        return cost
    catch
        return 100.0
    end
end

optimize(cost, 
    [3.0, 0.6], 
    NelderMead(), 
    Optim.Options(g_tol = 1e-6)
)