export Basis, GetConfWeight


"""
    Basis

The mutable type `Basis` stores the information of the basis that respects both diagonal and off-diagonal quantum numbers. The states in the basis is in the form 
```math
|I⟩=λ_{i_{I1}}|i_{I1}⟩+λ_{i_{I2}}|i_{I2}⟩+⋯+λ_{i_{Im_I}}|i_{Im_I}⟩
```
where ``|i⟩`` is a direct product state, _i. e._, the configurations ``|i_{Ik}⟩`` are grouped into a state ``|I⟩``. 

# Fields
* `cfs :: Confs` stores the configurations that respect the QNDiags.
* `dim :: Int64` is the dimension of the basis.
* `szz :: Int64` records the maximum size ``\\max m_I`` of groups.
* `cfgr :: Vector{Int64}` is a vector of length `cfs.ncf` and records which group ``|I⟩`` each configuration ``|i⟩`` belong to.
* `cffac :: Vector{ComplexF64}` is a vector of length `cfs.ncf` and records the coefficients ``λ_i`` of each configuration.
* `grel :: Matrix{Int64}` is a `szz`×`dim` matrix that records the configurations in each group ``|i_{Ik}⟩ (k = 1,…,m_I)``.
* `grsz :: Vector{Int64}` is a vector of length `dim` that records the size ``m_I`` of each group.
"""
mutable struct Basis
    cfs :: Confs
    dim :: Int64
    szz :: Int64
    cfgr :: Vector{Int64}
    cffac :: Vector{ComplexF64}
    grel :: Matrix{Int64}
    grsz :: Vector{Int64}
end 


"""
    Basis(cfs :: Confs, secf :: Vector{ComplexF64}, qnf :: Vector{QNOffd} ; num_th :: Int64, disp_std :: Bool)

generates the basis that respects the off-diagonal ``ℤ_p`` quantum numbers (QNOffd)

# Arguments 

* `cfs :: Confs` is the diagonal QN–preserving configurations.
* `secf :: Vector{ComplexF64}` is a vector of length the same as the number of discrete symmetries that records the eigenvalue of each transformation in the sector.
* `qnf :: Vector{QNOffd}` is a vector of off-diagonal quantum numbers.
* `num_th :: Int64`, the number of threads. Facultative, `NumThreads` by default.
* `disp_std :: Bool`, whether or not the log shall be displayed. Facultative, `!SilentStd` by default. 

# Output

* `bs :: Basis` is the resulting [Basis](@ref Basis) object.
"""
function Basis(cfs :: Confs, secf :: Vector{<:Number}, qnf :: Vector{QNOffd} ; num_th = NumThreads, disp_std = !SilentStd)
    if (length(secf) == 0) return Basis(cfs) end
    nqnf = length(secf)
    cyc = [ qnfi.cyc for qnfi in qnf ]
    perm_o_mat = reduce(hcat, [ qnfi.perm for qnfi in qnf ])
    ph_o_mat = reduce(hcat, [ qnfi.ph for qnfi in qnf ])
    fac_o_mat = reduce(hcat, [ qnfi.fac for qnfi in qnf ])
    dim_ref = Ref{Int64}(0)
    cfgr = Vector{Int64}(undef, cfs.ncf)
    cffac = Vector{ComplexF64}(undef, cfs.ncf)
    szz = prod([ abs(secf[i]) < 1E-8 ? 1 : cyc[i] for i = 1 : nqnf ])
    @ccall Libpath.__bs_MOD_generate_bs_cfgr(
        cfs.no :: Ref{Int64}, cfs.nor :: Ref{Int64}, cfs.ncf :: Ref{Int64}, cfs.lid :: Ref{Int64}, cfs.rid :: Ref{Int64}, cfs.conf :: Ref{Int64}, 
        nqnf :: Ref{Int64}, ComplexF64.(secf) :: Ref{ComplexF64}, 
        cyc :: Ref{Int64}, perm_o_mat :: Ref{Int64}, ph_o_mat :: Ref{Int64}, fac_o_mat :: Ref{ComplexF64}, 
        szz :: Ref{Int64}, dim_ref :: Ref{Int64}, cfgr :: Ref{Int64}, cffac :: Ref{ComplexF64}, 
        num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}
    ) :: Nothing
    dim = dim_ref[]
    grel = Matrix{Int64}(undef, szz, dim)
    grsz = Vector{Int64}(undef, dim)
    @ccall Libpath.__bs_MOD_generate_bs_grel(
        cfs.ncf :: Ref{Int64}, szz :: Ref{Int64}, dim :: Ref{Int64}, 
        cfgr :: Ref{Int64}, grel :: Ref{Int64}, grsz :: Ref{Int64}, 
        (disp_std ? 1 : 0) :: Ref{Int64}
    ) :: Nothing
    return Basis(cfs, dim, szz, cfgr, cffac, grel, grsz)
end 


"""
    Basis(cfs :: Confs)

Generate a basis from the configurations without off-diagonal ``ℤ_n`` symmetries.

# Arguments 

* `cfs :: Confs` is the diagonal QN–preserving configurations.

# Output

* `bs :: Basis` is the resulting `Basis` object.
"""
function Basis(cfs :: Confs)
    dim = cfs.ncf
    szz = 1 
    cfgr = collect(1 : dim)
    cffac = fill(ComplexF64(1), dim)
    grel = reshape(collect(1 : dim), 1, :)
    grsz = fill(1, dim)
    return Basis(cfs, dim, szz, cfgr, cffac, grel, grsz)
end 


"""
    GetConfWeight(bs :: Basis, st :: Vector{<:Number}, cf :: Int64) :: ComplexF64

looks up a the weight of a configuration in a state. 

# Arguments 

* `bs :: Basis` is the basis of the state.
* `st :: Vector{ComplexF64}` or `st :: Vector{Float64}` is a vector of length `bs.dim` that stores the state.
* `cf :: Int64` stores the configuration to be looked-up expressed in a binary number. If the `o-1`-th bit of `conf[i]` is 1, then the `o`-th site in the `i`-th configuration is occupied ; if the bit is 0, then the site is empty. 

# Output

* The weight of the configuration in the state.
"""
function GetConfWeight(bs :: Basis, st :: Union{Vector{ComplexF64}, Vector{Float64}}, cf :: Int64)
    id = GetConfId(bs.cfs, cf)
    (bs.cfgr[id] == -1) && return zero(st[1])
    return st[bs.cfgr[id]] * bs.cffac[id]
end
