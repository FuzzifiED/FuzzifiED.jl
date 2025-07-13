# This example calculates the entanglement entropy of the Ising ground state
# along the real space cut of θ = 0.500π and 0.499π respectively,
# and use these two data to extract finite size F-function without sustracting the IQHE contribution. 
# This example reproduces Figures 3 in Phys. Rev. B 111, 155151 (2025).
# On my portable computer, this calculation takes 4.197 s

using FuzzifiED
using SpecialFunctions
const σ1 = [  1  0 ;  0  0 ]
const σ2 = [  0  0 ;  0  1 ]
const σx = [  0  1 ;  1  0 ]
FuzzifiED.ElementType = Float64
≈(x, y) = abs(x - y) < eps(Float32)

let

nm = 10
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

cfs = Confs(2 * nm, [nm, 0], qnd)
bs = Basis(cfs, [1, 1, 1], qnf)
hmt = Operator(bs, tms_hmt)
hmt_mat = OpMat(hmt)
enrg, st = GetEigensystem(hmt_mat, 3)
st_g = st[:, 1]

secf_lst = [ [[1], [1]], [[-1], [-1]] ]
secd_lst = Vector{Vector{Int64}}[]
for nea = 0 : nm 
    neb = nm - nea 
    for lza = -min(nea, neb) * (nm - 1) : 2 : min(nea, neb) * (nm - 1)
        lzb = -lza 
        push!(secd_lst, [[nea, lza], [neb, lzb]])
    end
end

ent_entropy = Dict{Float64, Float64}()
θ1 = 0.500 * π
θ2 = 0.499 * π
qnd_a = qnd
qnf_a = [GetFlavPermQNOffd(nm, 2, [2, 1])]
for θ in (θ1, θ2)
    amp_oa = [ sqrt(beta_inc(m, nm - m + 1, (cos(θ) + 1) / 2)[1]) for m = 1 : nm for f = 1 : 2]
    ent_spec = GetEntSpec(st_g, bs, secd_lst, secf_lst ; qnd_a, qnf_a, amp_oa)

    eig_rho = vcat(values(ent_spec)...)
    tr_rho = sum(eig_rho)
    @show θ, tr_rho
    ent_entropy[θ] = -sum(eig_rho .* log.(eig_rho))
    @show θ, ent_entropy[θ]
end
F = (ent_entropy[θ1] * sin(θ2) - ent_entropy[θ2] * sin(θ1)) / (sin(θ1) - sin(θ2))
@show F

end