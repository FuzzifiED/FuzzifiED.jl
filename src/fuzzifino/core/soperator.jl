import FuzzifiED.AbstractOperator
export SOperator


"""
    SOperator
    
The mutable type `SOperator` records the sum of terms together with information about its symmetry and the basis of the state it acts on and the basis of the resulting state.

# Fields
* `bsd :: SBasis` is the basis of the initial state.
* `bsf :: SBasis` is the basis of the final state.
* `red_q :: Int64` is a flag that records whether or not the conversion to a sparse martrix can be simplified : if `bsd` and `bsf` have exactly the same set of quantum numbers, and the operator fully respects the symmetries, and all the elements in `bsd.cffac` and `bsf.cffac` has the same absolute value, then `red_q = 1` ; otherwise `red_q = 0`.
* `sym_q :: Int64` records the symmetry of the operator : if the matrix is Hermitian, then `sym_q = 1` ; if it is symmetric, then `sym_q = 2` ; otherwise `sym_q = 0`.
* `ntm :: Int64` is the number of terms.
* `nc :: Int64` is the maximum number of operators in an operator string
* `cstrs :: Matrix{Int64}` is a matrix recording the operator string of each term. Each column corresponds to a term and is padded to the maximum length with `-1`'s.
* `coeffs :: Vector{ComplexF64}` corresponds to the coefficients in each term.
"""
mutable struct SOperator <: AbstractOperator
    bsd :: SBasis 
    bsf :: SBasis
    red_q :: Int64
    sym_q :: Int64
    ntm :: Int64
    nc :: Int64
    cstrs :: Matrix{Int64}
    coeffs :: Vector{ComplexF64}
end 


"""
    SOperator(bsd :: SBasis[, bsf :: SBasis], terms :: STerms ; red_q :: Int64, sym_q :: Int64, num_th :: Int64, disp_std :: Bool) :: SOperator

generates an operator object from a series of terms. 

# Arguments
* `bsd :: SBasis` is the basis of the initial state.
* `bsf :: SBasis` is the basis of the final state. Facultative, the same as `bsd` by default. 
* `terms :: STerms` records the terms.
* `red_q :: Int64` is a flag that records whether or not the conversion to a sparse martrix can be simplified : if `bsd` and `bsf` have exactly the same set of quantum numbers, and the operator fully respects the symmetries, and all the elements in `bsd.cffac` and `bsf.cffac` has the same absolute value, then `red_q = 1` ; otherwise `red_q = 0` ; Facultative, if `bsf` is not given, 1 by default, otherwise 0 by default.
* `sym_q :: Int64` records the symmetry of the operator : if the matrix is Hermitian, then `sym_q = 1` ; if it is symmetric, then `sym_q = 2` ; otherwise `sym_q = 0`. Facultative, if `bsf` is not given, 1 by default, otherwise 0 by default.
* `num_th :: Int64`, the number of threads. Facultative, `NumThreads` by default.
* `disp_std :: Bool`, whether or not the log shall be displayed. Facultative, `!SilentStd` by default. 
"""
function SOperator(bsd :: SBasis, bsf :: SBasis, terms :: STerms ; red_q :: Int64 = 0, sym_q :: Int64 = 0)
    ntm = length(terms)
    nc = div(maximum([length(tm.cstr) for tm in terms]), 2)
    cstrs_vec = [ [tm.cstr ; fill(-1, 2 * nc - length(tm.cstr))] for tm in terms]
    cstrs = reduce(hcat, cstrs_vec)
    coeffs = [ tm.coeff for tm in terms ]
    return SOperator(bsd, bsf, red_q, sym_q, ntm, nc, cstrs, coeffs)
end
SOperator(bsd :: SBasis, terms :: STerms ; red_q :: Int64 = 1, sym_q :: Int64 = 1) = SOperator(bsd, bsd, terms ; red_q, sym_q)


"""
    *(op :: SOperator, st_d :: Vector{ComplexF64} ; num_th :: Int64, disp_std :: Bool) :: Vector{ComplexF64}
    *(op :: SOperator, st_d :: Vector{Float64} ; num_th :: Int64, disp_std :: Bool) :: Vector{Float64}

Measure the action of an operator on a state. `st_d` must be of length `op.bsd.dim`. Returns a vector of length `op.bsf.dim` that represents the final state.

# Facultative arguments

* `num_th :: Int64`, the number of threads. Facultative, `NumThreads` by default.
* `disp_std :: Bool`, whether or not the log shall be displayed. Facultative, `!SilentStd` by default. 
"""
function Base.:*(op :: SOperator, st_d :: Vector{ComplexF64} ; num_th = NumThreads, disp_std = !SilentStd)
    st_f = Vector{ComplexF64}(undef, op.bsf.dim)
    binom_d = [ binomial(i + j, i) for i = 0 : op.bsd.cfs.nebm, j = 0 : op.bsd.cfs.nob]
    binom_f = [ binomial(i + j, i) for i = 0 : op.bsf.cfs.nebm, j = 0 : op.bsf.cfs.nob]
    @ccall Libpathino.__sop_MOD_action_sop(op.bsd.cfs.nof :: Ref{Int64}, op.bsd.cfs.nob :: Ref{Int64}, op.bsd.cfs.norf :: Ref{Int64}, op.bsd.cfs.norb :: Ref{Int64}, 
    op.bsd.cfs.nebm :: Ref{Int64}, op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conff :: Ref{Int64}, op.bsd.cfs.confb :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, binom_d :: Ref{Int64}, 
    op.bsf.cfs.nebm :: Ref{Int64}, op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conff :: Ref{Int64}, op.bsf.cfs.confb :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, binom_f :: Ref{Int64}, 
    op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstrs :: Ref{Int64}, op.coeffs :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, ComplexF64.(st_d) :: Ref{ComplexF64}, st_f :: Ref{ComplexF64}, num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}) :: Nothing 
    return st_f
end
function Base.:*(op :: SOperator, st_d :: Vector{Float64} ; num_th = NumThreads, disp_std = !SilentStd)
    st_f = Vector{ComplexF64}(undef, op.bsf.dim)
    binom_d = [ binomial(i + j, i) for i = 0 : op.bsd.cfs.nebm, j = 0 : op.bsd.cfs.nob]
    binom_f = [ binomial(i + j, i) for i = 0 : op.bsf.cfs.nebm, j = 0 : op.bsf.cfs.nob]
    @ccall Libpathino.__sop_MOD_action_sop(op.bsd.cfs.nof :: Ref{Int64}, op.bsd.cfs.nob :: Ref{Int64}, op.bsd.cfs.norf :: Ref{Int64}, op.bsd.cfs.norb :: Ref{Int64}, 
    op.bsd.cfs.nebm :: Ref{Int64}, op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conff :: Ref{Int64}, op.bsd.cfs.confb :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, binom_d :: Ref{Int64}, 
    op.bsf.cfs.nebm :: Ref{Int64}, op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conff :: Ref{Int64}, op.bsf.cfs.confb :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, binom_f :: Ref{Int64}, 
    op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstrs :: Ref{Int64}, op.coeffs :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, ComplexF64.(st_d) :: Ref{ComplexF64}, st_f :: Ref{ComplexF64}, num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}) :: Nothing 
    return real.(st_f)
end


"""
    *(st_fp :: LinearAlgebra.Adjoint{ComplexF64, Vector{ComplexF64}}, op :: SOperator, st_d :: Vector{ComplexF64} ; num_th :: Int64, disp_std :: Bool) :: ComplexF64
    *(st_fp :: LinearAlgebra.Adjoint{Float64, Vector{Float64}}, op :: SOperator, st_d :: Vector{Float64} ; num_th :: Int64, disp_std :: Bool) :: Float64

Measuring the inner product between two states and an operator. `st_d` must be of length `op.bsd.dim` and `st_fp` must be of length `op.bsf.dim`, and `st_fp` must be an adjoint. 

# Facultative arguments

* `num_th :: Int64`, the number of threads. Facultative, `NumThreads` by default.
* `disp_std :: Bool`, whether or not the log shall be displayed. Facultative, `!SilentStd` by default. 
"""
function Base.:*(st_fp :: LinearAlgebra.Adjoint{ComplexF64, Vector{ComplexF64}}, op :: SOperator, st_d :: Vector{ComplexF64} ; num_th = NumThreads, disp_std = !SilentStd)
    ovl_ref = Ref{ComplexF64}(0)
    binom_d = [ binomial(i + j, i) for i = 0 : op.bsd.cfs.nebm, j = 0 : op.bsd.cfs.nob]
    binom_f = [ binomial(i + j, i) for i = 0 : op.bsf.cfs.nebm, j = 0 : op.bsf.cfs.nob]
    @ccall Libpathino.__sop_MOD_overlap_sop(op.bsd.cfs.nof :: Ref{Int64}, op.bsd.cfs.nob :: Ref{Int64}, op.bsd.cfs.norf :: Ref{Int64}, op.bsd.cfs.norb :: Ref{Int64}, 
    op.bsd.cfs.nebm :: Ref{Int64}, op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conff :: Ref{Int64}, op.bsd.cfs.confb :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, binom_d :: Ref{Int64}, 
    op.bsf.cfs.nebm :: Ref{Int64}, op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conff :: Ref{Int64}, op.bsf.cfs.confb :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, binom_f :: Ref{Int64}, 
    op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstrs :: Ref{Int64}, op.coeffs :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, ComplexF64.(st_d) :: Ref{ComplexF64}, ComplexF64.(st_fp') :: Ref{ComplexF64}, ovl_ref :: Ref{ComplexF64}, num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}) :: Nothing 
    return ovl_ref[]
end
function Base.:*(st_fp :: LinearAlgebra.Adjoint{Float64, Vector{Float64}}, op :: SOperator, st_d :: Vector{Float64} ; num_th = NumThreads, disp_std = !SilentStd)
    ovl_ref = Ref{ComplexF64}(0)
    binom_d = [ binomial(i + j, i) for i = 0 : op.bsd.cfs.nebm, j = 0 : op.bsd.cfs.nob]
    binom_f = [ binomial(i + j, i) for i = 0 : op.bsf.cfs.nebm, j = 0 : op.bsf.cfs.nob]
    @ccall Libpathino.__sop_MOD_overlap_sop(op.bsd.cfs.nof :: Ref{Int64}, op.bsd.cfs.nob :: Ref{Int64}, op.bsd.cfs.norf :: Ref{Int64}, op.bsd.cfs.norb :: Ref{Int64}, 
    op.bsd.cfs.nebm :: Ref{Int64}, op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conff :: Ref{Int64}, op.bsd.cfs.confb :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, binom_d :: Ref{Int64}, 
    op.bsf.cfs.nebm :: Ref{Int64}, op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conff :: Ref{Int64}, op.bsf.cfs.confb :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, binom_f :: Ref{Int64}, 
    op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstrs :: Ref{Int64}, op.coeffs :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, ComplexF64.(st_d) :: Ref{ComplexF64}, ComplexF64.(st_fp') :: Ref{ComplexF64}, ovl_ref :: Ref{ComplexF64}, num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}) :: Nothing 
    return real(ovl_ref[])
end
