export SBasis


"""
    SBasis
    
The mutable type `SBasis` stores the information of the SBasis that respects both diagonal and off-diagonal quantum numbers. The states in the SBasis is in the form 
```math
|I⟩=λ_{i_{I1}}|i_{I1}⟩+λ_{i_{I2}}|i_{I2}⟩+⋯+λ_{i_{Im_I}}|i_{Im_I}⟩
```
where ``|i⟩`` is a direct product state, _i. e._, the configurations ``|i_{Ik}⟩`` are grouped into a state ``|I⟩``. 

# Fields
* `cfs :: SConfs` stores the configurations that respect the QNDiags.
* `dim :: Int64` is the dimension of the SBasis.
* `szz :: Int64` records the maximum size ``\\max m_g`` of groups.
* `cfgr :: Vector{Int64}` is a vector of length `cfs.ncf` and records which group ``|I⟩`` each configuration ``|i⟩`` belong to.
* `cffac :: Vector{ComplexF64}` is a vector of length `cfs.ncf` and records the coefficients ``λ_i`` of each configuration.
* `grel :: Matrix{Int64}` is a `szz`×`dim` matrix that records the configurations in each group ``|i_{Ik}⟩ (k = 1,…,m_I)``
* `grsz :: Vector{Int64}` is a vector of length `dim` that records the size ``m_I`` of each group.
"""
mutable struct SBasis
    cfs :: SConfs
    dim :: Int64
    szz :: Int64
    cfgr :: Vector{Int64}
    cffac :: Vector{ComplexF64}
    grel :: Matrix{Int64}
    grsz :: Vector{Int64}
end 


"""
    SBasis(cfs :: SConfs, secf :: Vector{ComplexF64}, qnf :: Vector{SQNOffd} ; num_th :: Int64, disp_std :: Bool)

generates the SBasis that respects the off-diagonal ``ℤ_p`` quantum numbers (secfQNOffd)

# Arguments 

* `cfs :: SConfs` is the diagonal QN–preserving configurations.
* `secf :: Vector{ComplexF64}` is a vector of length the same as the number of discrete symmetries that records the eigenvalue of each transformation in the sector.
* `qnf :: Vector{SQNOffd}` is a vector of off-diagonal quantum numbers.
* `num_th :: Int64`, the number of threads. Facultative, `NumThreads` by default.
* `disp_std :: Bool`, whether or not the log shall be displayed. Facultative, `!SilentStd` by default. 

# Output

* `bs :: SBasis` is the resulting [SBasis](@ref SBasis) object.
"""
function SBasis(cfs :: SConfs, secf :: Vector{<:Number}, qnf :: Vector{SQNOffd} ; num_th = NumThreads, disp_std = !SilentStd)
    if (length(secf) == 0) return SBasis(cfs) end
    nqnf = length(secf)
    cyc = [ qnfi.cyc for qnfi in qnf ]
    permf_o_mat = reduce(hcat, [ qnfi.permf for qnfi in qnf ])
    permb_o_mat = reduce(hcat, [ qnfi.permb for qnfi in qnf ])
    phf_o_mat = reduce(hcat, [ qnfi.phf for qnfi in qnf ])
    facf_o_mat = reduce(hcat, [ qnfi.facf for qnfi in qnf ])
    facb_o_mat = reduce(hcat, [ qnfi.facb for qnfi in qnf ])
    dim_ref = Ref{Int64}(0)
    cfgr = Vector{Int64}(undef, cfs.ncf)
    cffac = Vector{ComplexF64}(undef, cfs.ncf)
    szz = prod([ abs(secf[i]) < 1E-8 ? 1 : cyc[i] for i = 1 : nqnf ])
    binom = [ binomial(i + j, i) for i = 0 : cfs.nebm, j = 0 : cfs.nob]
    @ccall Libpathino.__sbs_MOD_generate_sbs_cfgr(
        cfs.nof :: Ref{Int64}, cfs.nob :: Ref{Int64}, cfs.norf :: Ref{Int64}, cfs.norb :: Ref{Int64}, cfs.nebm :: Ref{Int64},
        cfs.ncf :: Ref{Int64}, cfs.lid :: Ref{Int64}, cfs.rid :: Ref{Int64}, cfs.conff :: Ref{Int64}, cfs.confb :: Ref{Int64}, 
        nqnf :: Ref{Int64}, ComplexF64.(secf) :: Ref{ComplexF64}, 
        cyc :: Ref{Int64}, permf_o_mat :: Ref{Int64}, permb_o_mat :: Ref{Int64}, phf_o_mat :: Ref{Int64}, facf_o_mat :: Ref{ComplexF64}, facb_o_mat :: Ref{ComplexF64}, 
        szz :: Ref{Int64}, dim_ref :: Ref{Int64}, cfgr :: Ref{Int64}, cffac :: Ref{ComplexF64}, binom :: Ref{Int64}, 
        num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}
    ) :: Nothing
    dim = dim_ref[]
    grel = Matrix{Int64}(undef, szz, dim)
    grsz = Vector{Int64}(undef, dim)
    @ccall Libpathino.__sbs_MOD_generate_sbs_grel(
        cfs.ncf :: Ref{Int64}, szz :: Ref{Int64}, dim :: Ref{Int64}, 
        cfgr :: Ref{Int64}, grel :: Ref{Int64}, grsz :: Ref{Int64}, 
        (disp_std ? 1 : 0) :: Ref{Int64}
    ) :: Nothing
    return SBasis(cfs, dim, szz, cfgr, cffac, grel, grsz)
end 


"""
    SBasis(cfs :: SConfs)

Generate a SBasis from the configurations without off-diagonal ``ℤ_n`` symmetries.

# Arguments 

* `cfs :: SConfs` is the diagonal QN–preserving configurations.

# Output

* `bs :: SBasis` is the resulting `SBasis` object.
"""
function SBasis(cfs :: SConfs)
    dim = cfs.ncf
    szz = 1 
    cfgr = collect(1 : dim)
    cffac = fill(ComplexF64(1), dim)
    grel = reshape(collect(1 : dim), 1, :)
    grsz = fill(1, dim)
    return SBasis(cfs, dim, szz, cfgr, cffac, grel, grsz)
end 
