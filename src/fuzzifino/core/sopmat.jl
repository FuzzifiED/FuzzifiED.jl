import FuzzifiED: OpMat

function OpMat{T}(op :: SOperator ; num_th = NumThreads, disp_std = !SilentStd) where T <: ComplexF64
    colptr = Vector{Int64}(undef, op.bsd.dim + 1)
    nel_ref = Ref{Int64}(0)
    binom_d = [ binomial(i + j, i) for i = 0 : op.bsd.cfs.nebm, j = 0 : op.bsd.cfs.nob]
    binom_f = [ binomial(i + j, i) for i = 0 : op.bsf.cfs.nebm, j = 0 : op.bsf.cfs.nob]
    @ccall Libpathino.__sop_MOD_count_sop(op.bsd.cfs.nof :: Ref{Int64}, op.bsd.cfs.nob :: Ref{Int64}, op.bsd.cfs.norf :: Ref{Int64}, op.bsd.cfs.norb :: Ref{Int64}, 
    op.bsd.cfs.nebm :: Ref{Int64}, op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conff :: Ref{Int64}, op.bsd.cfs.confb :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, binom_d :: Ref{Int64}, 
    op.bsf.cfs.nebm :: Ref{Int64}, op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conff :: Ref{Int64}, op.bsf.cfs.confb :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, binom_f :: Ref{Int64}, 
        op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstrs :: Ref{Int64}, op.coeffs :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, op.sym_q :: Ref{Int64}, nel_ref :: Ref{Int64}, colptr :: Ref{Int64}, num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}) :: Nothing
    nel = nel_ref[]
    rowid = Vector{Int64}(undef, nel)
    elval = Vector{ComplexF64}(undef, nel)
    @ccall Libpathino.__sop_MOD_generate_sop(op.bsd.cfs.nof :: Ref{Int64}, op.bsd.cfs.nob :: Ref{Int64}, op.bsd.cfs.norf :: Ref{Int64}, op.bsd.cfs.norb :: Ref{Int64}, 
    op.bsd.cfs.nebm :: Ref{Int64}, op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conff :: Ref{Int64}, op.bsd.cfs.confb :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, binom_d :: Ref{Int64}, 
    op.bsf.cfs.nebm :: Ref{Int64}, op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conff :: Ref{Int64}, op.bsf.cfs.confb :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, binom_f :: Ref{Int64}, 
        op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstrs :: Ref{Int64}, op.coeffs :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, op.sym_q :: Ref{Int64}, nel_ref :: Ref{Int64}, colptr :: Ref{Int64}, rowid :: Ref{Int64}, elval :: Ref{ComplexF64}, num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}) :: Nothing
    return OpMat{T}(op.bsd.dim, op.bsf.dim, op.sym_q, nel, colptr, rowid, elval)
end

function OpMat{T}(op :: SOperator ; num_th = NumThreads, disp_std = !SilentStd) where T <: Float64
    colptr = Vector{Int64}(undef, op.bsd.dim + 1)
    nel_ref = Ref{Int64}(0)
    binom_d = [ binomial(i + j, i) for i = 0 : op.bsd.cfs.nebm, j = 0 : op.bsd.cfs.nob]
    binom_f = [ binomial(i + j, i) for i = 0 : op.bsf.cfs.nebm, j = 0 : op.bsf.cfs.nob]
    @ccall Libpathino.__sop_MOD_count_sop(op.bsd.cfs.nof :: Ref{Int64}, op.bsd.cfs.nob :: Ref{Int64}, op.bsd.cfs.norf :: Ref{Int64}, op.bsd.cfs.norb :: Ref{Int64}, 
    op.bsd.cfs.nebm :: Ref{Int64}, op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conff :: Ref{Int64}, op.bsd.cfs.confb :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, binom_d :: Ref{Int64}, 
    op.bsf.cfs.nebm :: Ref{Int64}, op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conff :: Ref{Int64}, op.bsf.cfs.confb :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, binom_f :: Ref{Int64}, 
        op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstrs :: Ref{Int64}, op.coeffs :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, op.sym_q :: Ref{Int64}, nel_ref :: Ref{Int64}, colptr :: Ref{Int64}, num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}) :: Nothing
    nel = nel_ref[]
    rowid = Vector{Int64}(undef, nel)
    elval = Vector{Float64}(undef, nel)
    @ccall Libpathino.__sop_MOD_generate_sop_re(op.bsd.cfs.nof :: Ref{Int64}, op.bsd.cfs.nob :: Ref{Int64}, op.bsd.cfs.norf :: Ref{Int64}, op.bsd.cfs.norb :: Ref{Int64}, 
    op.bsd.cfs.nebm :: Ref{Int64}, op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conff :: Ref{Int64}, op.bsd.cfs.confb :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, binom_d :: Ref{Int64}, 
    op.bsf.cfs.nebm :: Ref{Int64}, op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conff :: Ref{Int64}, op.bsf.cfs.confb :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, binom_f :: Ref{Int64}, 
        op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstrs :: Ref{Int64}, op.coeffs :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, op.sym_q :: Ref{Int64}, nel_ref :: Ref{Int64}, colptr :: Ref{Int64}, rowid :: Ref{Int64}, elval :: Ref{Float64}, num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}) :: Nothing
    return OpMat{T}(op.bsd.dim, op.bsf.dim, op.sym_q, nel, colptr, rowid, elval)
end


"""
    OpMat[{type}](op :: SOperator ; num_th :: Int64, disp_std :: Bool) :: OpMat{type}

Generates the sparse matrix from the operator. The parameter `type` is either `Float64` or `ComplexF64` ; it is facultative, given by `ElementType` by default. 

# Arguments 

* `op :: SOperator` is the operator.
* `type :: DataType` specifies the type of the matrix. It can either be `ComplexF64` or `Float64`. Facultative, the same as `ElementType` by default
* `num_th :: Int64`, the number of threads. Facultative, `NumThreads` by default.
* `disp_std :: Bool`, whether or not the log shall be displayed. Facultative, `!SilentStd` by default. 
"""
OpMat(op :: SOperator ; type :: DataType = ElementType, num_th = NumThreads, disp_std = !SilentStd) = OpMat{type}(op ; num_th, disp_std)
