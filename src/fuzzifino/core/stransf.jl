export STransf


"""
    STransf 

The mutable type `STransf` records a transformation in the same form as a [SQNOffd](@ref)
```math
𝒵:\\ c_o↦ α_o^* c^{(p_o)}_{π_o},  c_o^†↦ α_o c^{(1-p_o)}_{π_o}
```
together with information about its symmetry and the basis of the state it acts on and the basis of the resulting state.

# Fields
* `bsd :: SBasis` is the basis of the initial state.
* `bsf :: SBasis` is the basis of the final state.
* `permf :: Vector{Int64}`, permb :: Vector{Int64}`, `phf :: Vector{Int64}` and `facf :: Vector{ComplexF64}`, `facb :: Vector{ComplexF64}` records the transformation in the same form as a [SQNOffd](@ref).
"""
mutable struct STransf 
    bsd :: SBasis 
    bsf :: SBasis
    permf :: Vector{Int64}
    permb :: Vector{Int64}
    phf :: Vector{Int64}
    facf :: Vector{ComplexF64}
    facb :: Vector{ComplexF64}
end


"""
    STransf(bsd :: SBasis, bsf :: SBasis, qnf :: SQNOffd)

generates a transformation object from a SQNOffd. 

# Arguments
* `bsd :: SBasis` is the basis of the initial state.
* `bsf :: SBasis` is the basis of the final state. Facultative, the same as `bsd` by default. 
* `qnf :: SQNOffd` records the transformation ; 
"""
function STransf(bsd :: SBasis, bsf :: SBasis, qnf :: SQNOffd)
    return STransf(bsd, bsf, qnf.permf, qnf.permb, qnf.phf, qnf.facf, qnf.facb)
end
STransf(bsd :: SBasis, qnf :: SQNOffd) = STransf(bsd, bsd, qnf)


"""
    *(trs :: STransf, st_d :: Vector{ComplexF64} ; num_th = NumThreads) :: Vector{ComplexF64}
    *(trs :: STransf, st_d :: Vector{Float64} ; num_th = NumThreads) :: Vector{Float64}

Act a transformation on a state. `st_d` must be of length `trs.bsd.dim`. Returns a vector of length `trs.bsf.dim` that represents the final state.

# Facultative arguments

* `num_th :: Int64`, the number of threads. Facultative, `NumThreads` by default. 
"""
function Base.:*(trs :: STransf, st_d :: Vector{ComplexF64} ; num_th = NumThreads)
    st_f = Vector{ComplexF64}(undef, trs.bsf.dim)
    binom_d = [ binomial(i + j, i) for i = 0 : trs.bsd.cfs.nebm, j = 0 : trs.bsd.cfs.nob]
    binom_f = [ binomial(i + j, i) for i = 0 : trs.bsf.cfs.nebm, j = 0 : trs.bsf.cfs.nob]
    @ccall Libpathino.__sbs_MOD_action_strs(trs.bsd.cfs.nof :: Ref{Int64}, trs.bsd.cfs.nob :: Ref{Int64}, trs.bsd.cfs.norf :: Ref{Int64}, trs.bsd.cfs.norb :: Ref{Int64}, 
    trs.bsd.cfs.nebm :: Ref{Int64}, trs.bsd.cfs.ncf :: Ref{Int64}, trs.bsd.dim :: Ref{Int64}, trs.bsd.cfs.conff :: Ref{Int64}, trs.bsd.cfs.confb :: Ref{Int64}, trs.bsd.cfs.lid :: Ref{Int64}, trs.bsd.cfs.rid :: Ref{Int64}, trs.bsd.szz :: Ref{Int64}, trs.bsd.cfgr :: Ref{Int64}, trs.bsd.cffac :: Ref{ComplexF64}, trs.bsd.grel :: Ref{Int64}, trs.bsd.grsz :: Ref{Int64}, binom_d :: Ref{Int64}, 
    trs.bsf.cfs.nebm :: Ref{Int64}, trs.bsf.cfs.ncf :: Ref{Int64}, trs.bsf.dim :: Ref{Int64}, trs.bsf.cfs.conff :: Ref{Int64}, trs.bsf.cfs.confb :: Ref{Int64}, trs.bsf.cfs.lid :: Ref{Int64}, trs.bsf.cfs.rid :: Ref{Int64}, trs.bsf.szz :: Ref{Int64}, trs.bsf.cfgr :: Ref{Int64}, trs.bsf.cffac :: Ref{ComplexF64}, trs.bsf.grel :: Ref{Int64}, trs.bsf.grsz :: Ref{Int64}, binom_f :: Ref{Int64}, 
    trs.permf :: Ref{Int64}, trs.permb :: Ref{Int64}, trs.phf :: Ref{Int64}, trs.facf :: Ref{ComplexF64}, trs.facb :: Ref{ComplexF64}, ComplexF64.(st_d) :: Ref{ComplexF64}, st_f :: Ref{ComplexF64}, num_th :: Ref{Int64}) :: Nothing 
    return st_f
end
function Base.:*(trs :: STransf, st_d :: Vector{Float64} ; num_th = NumThreads)
    st_f = Vector{ComplexF64}(undef, trs.bsf.dim)
    binom_d = [ binomial(i + j, i) for i = 0 : trs.bsd.cfs.nebm, j = 0 : trs.bsd.cfs.nob]
    binom_f = [ binomial(i + j, i) for i = 0 : trs.bsf.cfs.nebm, j = 0 : trs.bsf.cfs.nob]
    @ccall Libpathino.__sbs_MOD_action_strs(trs.bsd.cfs.nof :: Ref{Int64}, trs.bsd.cfs.nob :: Ref{Int64}, trs.bsd.cfs.norf :: Ref{Int64}, trs.bsd.cfs.norb :: Ref{Int64}, 
    trs.bsd.cfs.nebm :: Ref{Int64}, trs.bsd.cfs.ncf :: Ref{Int64}, trs.bsd.dim :: Ref{Int64}, trs.bsd.cfs.conff :: Ref{Int64}, trs.bsd.cfs.confb :: Ref{Int64}, trs.bsd.cfs.lid :: Ref{Int64}, trs.bsd.cfs.rid :: Ref{Int64}, trs.bsd.szz :: Ref{Int64}, trs.bsd.cfgr :: Ref{Int64}, trs.bsd.cffac :: Ref{ComplexF64}, trs.bsd.grel :: Ref{Int64}, trs.bsd.grsz :: Ref{Int64}, binom_d :: Ref{Int64}, 
    trs.bsf.cfs.nebm :: Ref{Int64}, trs.bsf.cfs.ncf :: Ref{Int64}, trs.bsf.dim :: Ref{Int64}, trs.bsf.cfs.conff :: Ref{Int64}, trs.bsf.cfs.confb :: Ref{Int64}, trs.bsf.cfs.lid :: Ref{Int64}, trs.bsf.cfs.rid :: Ref{Int64}, trs.bsf.szz :: Ref{Int64}, trs.bsf.cfgr :: Ref{Int64}, trs.bsf.cffac :: Ref{ComplexF64}, trs.bsf.grel :: Ref{Int64}, trs.bsf.grsz :: Ref{Int64}, binom_f :: Ref{Int64}, 
    trs.permf :: Ref{Int64}, trs.permb :: Ref{Int64}, trs.phf :: Ref{Int64}, trs.facf :: Ref{ComplexF64}, trs.facb :: Ref{ComplexF64}, ComplexF64.(st_d) :: Ref{ComplexF64}, st_f :: Ref{ComplexF64}, num_th :: Ref{Int64}) :: Nothing 
    return real.(st_f)
end
