# This example calculates the scaling dimension of the cusp of the magnetic line defect 
# in 3d Ising model as a function of the angle θ.
# This example reproduces Table 2, upper panel in JHEP 11 (2024) 061
# On my portable computer, this calculation takes 33.841 s

using FuzzifiED
const σ1 = [  1  0 ;  0  0 ]
const σ2 = [  0  0 ;  0  1 ]
const σx = [  0  1 ;  1  0 ]
const σz = [  1  0 ;  0 -1 ]
FuzzifiED.ElementType = Float64
≈(x, y) = abs(x - y) < eps(Float32)

let

nm = 10
no = nm * 2
    
qnd = [ 
    GetNeQNDiag(2 * nm), 
    GetLz2QNDiag(nm, 2) 
]
qnf = [ 
    GetParityQNOffd(nm, 2, [2, 1], [-1, 1]), 
    GetFlavPermQNOffd(nm, 2, [2, 1]), 
    GetRotyQNOffd(nm, 2) 
]
cfs = Confs(2 * nm, [nm, 0], qnd)
bs = Basis(cfs, [1, 1, 1], qnf)
tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, 2, 2 .* [4.75, 1.], σ1, σ2) - 
    3.16 * GetPolTerms(nm, 2, σx) )
hmt = Operator(bs, tms_hmt)
hmt_mat = OpMat(hmt)
enrg, st = GetEigensystem(hmt_mat, 6)
enrg_0 = enrg[1]
enrg_T = enrg[3]

qnd = [ GetNeQNDiag(no) ]
cfs = Confs(no, [nm], qnd)
bs = Basis(cfs)
θs = collect(π : -0.1 * π : 0)
dims = []
enrg_d = 0
hd = 1000 
nz = StoreComps(GetDensityObs(nm, 2, σz))

for θ in θs
    tms_hmt = SimplifyTerms(
        GetDenIntTerms(nm, 2, 2 .* [4.75, 1.], σ1, σ2)
        - 3.16 * GetPolTerms(nm, 2, σx) 
        - hd / nm * GetPointValue(nz, 0.0, 0.0)
        - hd / nm * GetPointValue(nz, θ,   0.0)
    )
    hmt = Operator(bs, tms_hmt)
    hmt_mat = OpMat(hmt)
    enrg, st = GetEigensystem(hmt_mat, 3)

    if (θ ≈ π) enrg_d = enrg[1] end
    dim = 3 * (enrg[1] - enrg_d) / (enrg_T - enrg_0)
    append!(dims, dim)
end
display(hcat([θs, dims]...))

end