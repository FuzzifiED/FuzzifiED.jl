# This example calculates the correlator of n^i for ϕ, ϕ^2 and π of free real scalar.
# This example reproduces Figure 3 in arXiv : 2506.14904. 
# On my portable computer, this calculation takes 1.183 s

using FuzzifiED
using LegendrePolynomials
const σx = [  0  1 ;  1  0 ]
const σy = [  0 im ;-im  0 ]
const σz = [  1  0 ;  0 -1 ]
≈(x, y) = abs(x - y) < eps(Float32)

let ovl_IzI, ovl_ϕ2zI, ovl_ϕxI, ovl_ϕyI 

nm = 12
nf = 2
no = nm * nf
qnd = [ 
    GetNeQNDiag(no), 
    GetLz2QNDiag(nm, nf),
    GetFlavQNDiag(nm, nf, [0, 1], 1, 2)
]
qnf = [ 
    GetParityQNOffd(nm, nf, [2, 1], [-1, 1]), 
    GetRotyQNOffd(nm, nf) 
]

FuzzifiED.ObsNormRadSq = nm
n0 = GetDensityObs(nm, nf)
ny = GetDensityObs(nm, nf, σy)
nz = GetDensityObs(nm, nf, σz)
nx = GetDensityObs(nm, nf, σx)
tms_hmt = SimplifyTerms(
    GetIntegral(n0 * n0)
    + 0.18 * GetIntegral(ny * Laplacian(ny))
    - 0.42 / 4π / nm * GetIntegral(nz)
)
tms_l2 = GetL2Terms(nm, 2)

cfs = Dict{Int64, Confs}()
cfs[ 1] = Confs(2 * nm, [nm, 0, 0], qnd)
cfs[-1] = Confs(2 * nm, [nm, 0, 1], qnd)

bss = Dict{Tuple{Int64, Int64, Int64}, Basis}()
sts = Dict{Tuple{Int64, Int64, Int64}, Matrix{ComplexF64}}()
for P in [1], Z in [1, -1], R in [1, -1]
    bs = Basis(cfs[Z], [P, R], qnf)
    bss[P, Z, R] = bs
    
    if (R == -1) continue end 
    hmt = Operator(bs, tms_hmt)
    hmt_mat = OpMat(hmt)
    enrg, st = GetEigensystem(hmt_mat, 3)
    sts[P, Z, R] = st
end

stI = sts[1, 1, 1][:, 1]
stϕ = sts[1,-1, 1][:, 1]
stϕ2= sts[1, 1, 1][:, 2]

cor_z = Dict{Int64, Float64}()
for l = 0 : nm - 1
    nz_l = Operator(bss[1, 1, 1], bss[1, 1, (-1)^l], GetComponent(nz, l, 0.0))
    st_zI = nz_l * stI
    if (l == 0) 
        ovl_IzI = stI' * st_zI
        ovl_ϕ2zI = stϕ2' * st_zI
    end
    cor_l = st_zI' * st_zI / abs(ovl_ϕ2zI) ^ 2 
    if (l == 0)
        cor_l -= abs(ovl_IzI) ^ 2 / abs(ovl_ϕ2zI) ^ 2
    end 
    cor_z[l] = cor_l 
end
cor_z_fn(θ :: Float64) = sum([(2 * l + 1) * Pl(cos(θ), l) * cor for (l, cor) in cor_z])

cor_x = Dict{Int64, Float64}()
for l = 0 : nm - 1
    nx_l = Operator(bss[1, 1, 1], bss[1,-1, (-1)^l], GetComponent(nx, l, 0.0))
    st_xI = nx_l * stI
    if (l == 0) 
        ovl_ϕxI = stϕ' * st_xI
    end
    cor_l = st_xI' * st_xI / abs(ovl_ϕxI) ^ 2 
    cor_x[l] = cor_l 
end
cor_x_fn(θ :: Float64) = sum([(2 * l + 1) * Pl(cos(θ), l) * cor for (l, cor) in cor_x])

cor_y = Dict{Int64, Float64}()
for l = 0 : nm - 1
    ny_l = Operator(bss[1, 1, 1], bss[1,-1, (-1)^l], GetComponent(ny, l, 0.0))
    st_yI = ny_l * stI
    if (l == 0) 
        ovl_ϕyI = stϕ' * st_yI
    end
    cor_l = st_yI' * st_yI / abs(ovl_ϕyI) ^ 2 / 4
    cor_y[l] = cor_l 
end
cor_y_fn(θ :: Float64) = sum([(2 * l + 1) * Pl(cos(θ), l) * cor for (l, cor) in cor_y])

display(permutedims(hcat([[θ / π, cor_z_fn(θ), cor_x_fn(θ), cor_y_fn(θ)] for θ = 0 : π / 10 : 2 * π]...)))

end 