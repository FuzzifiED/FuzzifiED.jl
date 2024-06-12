# This example calculates the full spectrum of 3d Ising model on fuzzy sphere
# at nm = 10 for sector (P,Z,R) = (1,1,1)
# On my table computer, this calculation takes 0.789 s

using FuzzifiED
using LinearAlgebra
const σ1 = [  1  0 ;  0  0 ]
const σ2 = [  0  0 ;  0  1 ]
const σx = [  0  1 ;  1  0 ]

nm = 10
qnd = [ 
    GetNeQNDiag(2 * nm), 
    GetLz2QNDiag(nm, 2) ]
qnf = [ 
    GetParityQNOffd(nm, 2, [2, 1], [-1, 1]), 
    GetFlavPermQNOffd(nm, 2, [2, 1]), 
    GetRotyQNOffd(nm, 2) ]
tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, 2 ; ps_pot = 2 .* [4.75, 1.], mat_a = σ1, mat_b = σ2) - 
    3.16 * GetPolTerms(nm, 2 ; mat = σx) )

cfs = Confs(2 * nm, [nm, 0], qnd)
bs = Basis(cfs, [1, 1, 1], qnf)
hmt = Operator(bs, bs, tms_hmt ; red_q = 1, sym_q = 1)
hmt_mat = OpMat(hmt ; type = Float64)
hmt_mat_full = MatrixFromOpMat(hmt_mat)
enrg, st = eigen(hmt_mat_full)
@show enrg