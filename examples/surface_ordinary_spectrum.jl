# This example calculates the spectrum of ordinary surface CFT in 3d Ising model
# calibrated by surface displacement operator D in the orbital boundary scheme.
# This example reproduces Figures 3 and 4 in  SciPost Phys. 18, 031 (2025).
# On my portable computer, this calculation takes 2.196 s

using FuzzifiED
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
qnf = [ GetFlavPermQNOffd(nm, 2, [2, 1]) ]

ps_pot = 2 .* [4.75, 1.]
fld_h = 3.16
int_mat = GetIntMatrix(2 * nm, ps_pot)
fld_pin = [ sum([ int_mat[m1, m2, m2] for m2 = nm + 1 : 2 * nm]) for m1 = 1 : nm ]

tms_hmt = SimplifyTerms(
    RemoveOrbs(GetDenIntTerms(nm * 2, 2, 2 .* [4.75, 1.], σ1, σ2), collect(2 * nm + 1 : 4 * nm))
    + GetPolTerms(nm, 2 ; fld_m = fld_pin) / 2
    - 3.16 * GetPolTerms(nm, 2, σx) 
)

enrg_0 = 0
enrg_D = 0
spec = []
for lz = 0 : 4, Z in (1, -1)
    cfs = Confs(2 * nm, [nm, 2 * lz], qnd)
    bs = Basis(cfs, [Z], qnf)
    hmt = Operator(bs, tms_hmt)
    hmt_mat = OpMat(hmt)
    enrg, st = GetEigensystem(hmt_mat, 20)
    if (lz == 0 && Z == 1)
        enrg_0 = enrg[1]
        enrg_D = enrg[2]
    end
    dim = 3 .* (enrg .- enrg_0) ./ (enrg_D - enrg_0)
    for i in eachindex(enrg)
        push!(spec, round.([dim[i], enrg[i], lz, Z] .+ √eps(Float64), digits = 6))
    end
end

sort!(spec, by = st -> real(st[1]))
display(permutedims(hcat(spec...)))

end