# This example calculates the entanglement entropy of the Ising ground state
# along the orbital space cut at m = 0, and also the lowest entanglement spectrum 
# in the half-filled sectors
# On my table computer, this calculation takes 0.506 s

using FuzzifiED
const σ1 = [  1  0 ;  0  0 ]
const σ2 = [  0  0 ;  0  1 ]
const σx = [  0  1 ;  1  0 ]
FuzzifiED.ElementType = Float64
≈(x, y) = abs(x - y) < eps(Float32)

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

cfs = Confs(2 * nm, [nm, 0], qnd)
bs = Basis(cfs, [1, 1, 1], qnf)
hmt = Operator(bs, tms_hmt)
hmt_mat = OpMat(hmt)
enrg, st = GetEigensystem(hmt_mat, 3)
st_g = st[:, 1]

secf_lst = [ [[1], [1]], [[-1], [-1]] ]
secd_lst = Vector{Vector{Int64}}[]
qnd_a = [ qnd ; GetPinOrbQNDiag(2 * nm, collect(nm + 1 : 2 * nm)) ]
qnd_b = [ qnd ; GetPinOrbQNDiag(2 * nm, collect(1 : nm)) ]
amp_oa = [ o <= nm ? 1 : 0 for o = 1 : 2 * nm]

for nea = 0 : nm 
    neb = nm - nea 
    for lza = -min(nea, neb) * (nm - 1) : 2 : 0
        lzb = -lza 
        push!(secd_lst, [[nea, lza, 0], [neb, lzb, 0]])
    end
end
ent_spec = GetEntSpec(st_g, bs, secd_lst, secf_lst ; qnd_a, qnd_b, qnf_a = [GetFlavPermQNOffd(nm, 2, [2, 1])], amp_oa)

eig_rho = vcat(values(ent_spec)...)
tr_rho = sum(eig_rho)
@show tr_rho
ent_entropy = -sum(eig_rho .* log.(eig_rho))
@show ent_entropy

ent_spec_p = hcat(sort([ [sec.secd_a[2] / 2, -log(val)] 
    for (sec, vals) in ent_spec 
    if sec.secd_a[1] == nm ÷ 2 && sec.secf_a[1] == 1
    for val in vals 
], by = p -> p[2])...)
ent_spec_m = hcat(sort([ [sec.secd_a[2] / 2, -log(val)] 
    for (sec, vals) in ent_spec 
    if sec.secd_a[1] == nm ÷ 2 && sec.secf_a[1] == -1
    for val in vals 
], by = p -> p[2])...)
display(permutedims(ent_spec_p[:, 1:20]))
display(permutedims(ent_spec_m[:, 1:20]))
# To plot the entanglement spectrum, uncomment the following lines 
# using Plots 
# scatter(ent_spec_p[1, :], ent_spec_p[2, :])
# scatter!(ent_spec_m[1, :], ent_spec_m[2, :])