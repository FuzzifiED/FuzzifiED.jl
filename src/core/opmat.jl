export OpMat, GetEigensystem


"""
    OpMat{ComplexF64}
    OpMat{Float64}

The mutable type `OpMat{T}` where the element type `T` can be Float64 and ComplexF64 stores a sparse matrix in the same form as `SparseMatrixCSC` in `SparseArrays`. If the matrix is Hermitian or symmetric, only the lower triangle is stored. 

# Fields
* `dimd :: Int64` and `dimf :: Int64` are the number of columns and rows of the matrix.
* `symq :: Int64` records whether or not the matrix is Hermitian or symmetric.
* `nel :: Int64` records the number of elements.
* `colptr :: Vector{Int64}`, `rowid :: Vector{Int64}` and `elval :: Vector{ComplexF64}` records the elements of the sparse matrix as in the `SparseMatrixCSC` elements of Julia. 
"""
mutable struct OpMat{T <: Union{Float64, ComplexF64}}
    dimd :: Int64
    dimf :: Int64
    sym_q :: Int64
    nel :: Int64 
    colptr :: Vector{Int64}
    rowid :: Vector{Int64}
    elval :: Vector{T}

    function OpMat{T}(dimd :: Int64, dimf :: Int64, sym_q :: Int64, nel :: Int64, colptr :: Vector{Int64}, rowid :: Vector{Int64}, elval :: Vector{T}) where T <: Union{Float64, ComplexF64}
        return new(dimd, dimf, sym_q, nel, colptr, rowid, elval)
    end
    function OpMat{ComplexF64}(op :: Operator ; num_th = NumThreads, disp_std = !SilentStd)
        colptr = Vector{Int64}(undef, op.bsd.dim + 1)
        nel_ref = Ref{Int64}(0)
        @ccall Libpath.__op_MOD_count_op(op.bsd.cfs.no :: Ref{Int64}, op.bsd.cfs.nor :: Ref{Int64}, 
            op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conf :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, 
            op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conf :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, 
            op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstrs :: Ref{Int64}, op.coeffs :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, op.sym_q :: Ref{Int64}, nel_ref :: Ref{Int64}, colptr :: Ref{Int64}, num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}) :: Nothing
        nel = nel_ref[]
        rowid = Vector{Int64}(undef, nel)
        elval = Vector{ComplexF64}(undef, nel)
        @ccall Libpath.__op_MOD_generate_op(op.bsd.cfs.no :: Ref{Int64}, op.bsd.cfs.nor :: Ref{Int64}, 
            op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conf :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, 
            op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conf :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, 
            op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstrs :: Ref{Int64}, op.coeffs :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, op.sym_q :: Ref{Int64}, nel_ref :: Ref{Int64}, colptr :: Ref{Int64}, rowid :: Ref{Int64}, elval :: Ref{ComplexF64}, num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}) :: Nothing
        return new(op.bsd.dim, op.bsf.dim, op.sym_q, nel, colptr, rowid, elval)
    end
    function OpMat{Float64}(op :: Operator ; num_th = NumThreads, disp_std = !SilentStd)
        colptr = Vector{Int64}(undef, op.bsd.dim + 1)
        nel_ref = Ref{Int64}(0)
        @ccall Libpath.__op_MOD_count_op(op.bsd.cfs.no :: Ref{Int64}, op.bsd.cfs.nor :: Ref{Int64}, 
            op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conf :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, 
            op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conf :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, 
            op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstrs :: Ref{Int64}, op.coeffs :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, op.sym_q :: Ref{Int64}, nel_ref :: Ref{Int64}, colptr :: Ref{Int64}, num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}) :: Nothing
        nel = nel_ref[]
        rowid = Vector{Int64}(undef, nel)
        elval = Vector{Float64}(undef, nel)
        @ccall Libpath.__op_MOD_generate_op_re(op.bsd.cfs.no :: Ref{Int64}, op.bsd.cfs.nor :: Ref{Int64}, 
            op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conf :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, 
            op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conf :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, 
            op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstrs :: Ref{Int64}, op.coeffs :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, op.sym_q :: Ref{Int64}, nel_ref :: Ref{Int64}, colptr :: Ref{Int64}, rowid :: Ref{Int64}, elval :: Ref{Float64}, num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}) :: Nothing
        return new(op.bsd.dim, op.bsf.dim, op.sym_q, nel, colptr, rowid, elval)
    end
end 
"""
    OpMat[{T}](op :: Operator ; num_th :: Int64, disp_std :: Bool)

Generates the sparse matrix from the operator. The parametric type `T` is either `Float64` or `ComplexF64` ; it is facultative, given by `ElementType` by default. 

# Arguments 

* `op :: Operator` is the operator.
* `T :: DataType` specifies the type of the matrix. It can either be `ComplexF64` or `Float64`. Facultative, the same as `ElementType` by default.
* `num_th :: Int64`, the number of threads. Facultative, `NumThreads` by default.
* `disp_std :: Bool`, whether or not the log shall be displayed. Facultative, `!SilentStd` by default. 
"""
OpMat(op :: Operator ; type :: DataType = ElementType, num_th = NumThreads, disp_std = !SilentStd) = OpMat{type}(op ; num_th, disp_std)


"""
    GetEigensystem(mat :: OpMat{ComplexF64}, nst :: Int64 ; tol :: Float64, ncv :: Int64, initvec :: Vector{ComplexF64}, num_th :: Int64, disp_std :: Bool) :: Tuple{Vector{ComplexF64}, Matrix{ComplexF64}}
    GetEigensystem(mat :: OpMat{Float64}, nst :: Int64 ; tol :: Float64, ncv :: Int64, initvec :: Vector{Float64}, num_th :: Int64, disp_std :: Bool) :: Tuple{Vector{Float64}, Matrix{Float64}}

calls the Arpack package to calculate the lowest eigenstates of sparse matrix. 

# Arguments 

* `mat :: OpMat{ComplexF64}` or `mat :: OpMat{Float64}` is the matrix.
* `nst :: Int64` is the number of eigenstates to be calculated.
* `tol :: Float64` is the tolerence for the Arpack process. The default value is `1E-8`.
* `ncv :: Int64` is an auxiliary parameter needed in the Arpack process. The default value is `max(2 * nst, nst + 10)`.
* `initvec :: Vector{ComplexF64}` or `initvec :: Vector{Float64}` is the initial vector. If empty, a random initialisation shall be used. Facultative, empty by default.
* `num_th :: Int64`, the number of threads. Facultative, `NumThreads` by default.
* `disp_std :: Bool`, whether or not the log shall be displayed. Facultative, `!SilentStd` by default. 

# Output

* A length-`nst` array that has the same type as `mat` recording the eigenvalues, and 
* A `dimd`×`nst` matrix that has the same type as `mat` where every column records an eigenstate. 
"""
function GetEigensystem(mat :: OpMat{ComplexF64}, nst :: Int64 ; tol :: Float64 = 1E-8, ncv :: Int64 = min(max(2 * nst, nst + 10), mat.dimd), initvec :: Vector{ComplexF64} = ComplexF64[], num_th = NumThreads, disp_std = !SilentStd)
    eigval = Vector{ComplexF64}(undef, nst + 1)
    eigvec = Matrix{ComplexF64}(undef, mat.dimd, nst)
    if isempty(initvec)
        @ccall Libpath.__diag_MOD_diagonalisation(mat.dimd :: Ref{Int64}, mat.sym_q :: Ref{Int64}, mat.nel :: Ref{Int64}, mat.colptr :: Ref{Int64}, mat.rowid :: Ref{Int64}, mat.elval :: Ref{ComplexF64}, nst :: Ref{Int64}, tol :: Ref{Float64}, ncv :: Ref{Int64}, eigval :: Ref{ComplexF64}, eigvec :: Ref{ComplexF64}, num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}) :: Nothing
    else
        @ccall Libpath.__diag_MOD_diagonalisation_init(mat.dimd :: Ref{Int64}, mat.sym_q :: Ref{Int64}, mat.nel :: Ref{Int64}, mat.colptr :: Ref{Int64}, mat.rowid :: Ref{Int64}, mat.elval :: Ref{ComplexF64}, nst :: Ref{Int64}, tol :: Ref{Float64}, ncv :: Ref{Int64}, initvec :: Ref{ComplexF64}, eigval :: Ref{ComplexF64}, eigvec :: Ref{ComplexF64}, num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}) :: Nothing
    end
    return eigval[1 : end - 1], eigvec
end 
function GetEigensystem(mat :: OpMat{Float64}, nst :: Int64 ; tol :: Float64 = 1E-8, ncv :: Int64 = min(max(2 * nst, nst + 10), mat.dimd), initvec :: Vector{Float64} = Float64[], num_th = NumThreads, disp_std = !SilentStd)
    eigval = Vector{Float64}(undef, nst + 1)
    eigvec = Matrix{Float64}(undef, mat.dimd, nst)
    if isempty(initvec)
        @ccall Libpath.__diag_re_MOD_diagonalisation_re(mat.dimd :: Ref{Int64}, mat.sym_q :: Ref{Int64}, mat.nel :: Ref{Int64}, mat.colptr :: Ref{Int64}, mat.rowid :: Ref{Int64}, mat.elval :: Ref{Float64}, nst :: Ref{Int64}, tol :: Ref{Float64}, ncv :: Ref{Int64}, eigval :: Ref{Float64}, eigvec :: Ref{Float64}, num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}) :: Nothing
    else
        @ccall Libpath.__diag_re_MOD_diagonalisation_init_re(mat.dimd :: Ref{Int64}, mat.sym_q :: Ref{Int64}, mat.nel :: Ref{Int64}, mat.colptr :: Ref{Int64}, mat.rowid :: Ref{Int64}, mat.elval :: Ref{Float64}, nst :: Ref{Int64}, tol :: Ref{Float64}, ncv :: Ref{Int64}, initvec :: Ref{Float64}, eigval :: Ref{Float64}, eigvec :: Ref{Float64}, num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}) :: Nothing
    end
    return eigval[1 : end - 1], eigvec
end 

"""
    *(mat :: OpMat{ComplexF64}, st_d :: Vector{ComplexF64} ; num_th :: Int64) :: Vector{ComplexF64}
    *(mat :: OpMat{Float64}, st_d :: Vector{Float64} ; num_th :: Int64) :: Vector{Float64}

Measure the action of a sparse matrix on a state. `st_d` must be of length `mat.dimd`. Returns a vector of length `mat.dimf` that represents the final state. 

# Facultative argument

* `num_th :: Int64`, the number of threads. Facultative, `NumThreads` by default. 
"""
function Base.:*(mat :: OpMat{ComplexF64}, st_d :: Vector{ComplexF64} ; num_th = NumThreads)
    st_f = Vector{ComplexF64}(undef, mat.dimf)
    @ccall Libpath.__diag_MOD_vec_prod(mat.dimd :: Ref{Int64}, mat.dimf :: Ref{Int64}, mat.sym_q :: Ref{Int64}, mat.nel :: Ref{Int64}, mat.colptr :: Ref{Int64}, mat.rowid :: Ref{Int64}, mat.elval :: Ref{ComplexF64}, st_d :: Ref{ComplexF64}, st_f :: Ref{ComplexF64}, num_th :: Ref{Int64}) :: Nothing
    return st_f
end
function Base.:*(mat :: OpMat{Float64}, st_d :: Vector{Float64} ; num_th = NumThreads)
    st_f = Vector{Float64}(undef, mat.dimf)
    @ccall Libpath.__diag_re_MOD_vec_prod_re(mat.dimd :: Ref{Int64}, mat.dimf :: Ref{Int64}, mat.sym_q :: Ref{Int64}, mat.nel :: Ref{Int64}, mat.colptr :: Ref{Int64}, mat.rowid :: Ref{Int64}, mat.elval :: Ref{Float64}, real.(st_d) :: Ref{Float64}, st_f :: Ref{Float64}, num_th :: Ref{Int64}) :: Nothing
    return st_f
end


"""
    *(st_fp :: LinearAlgebra.Adjoint{ComplexF64, Vector{ComplexF64}}, mat :: OpMat{ComplexF64}, st_d :: Vector{ComplexF64} ; num_th :: Int64) :: ComplexF64
    *(st_fp :: LinearAlgebra.Adjoint{Float64, Vector{Float64}}, mat :: OpMat{Float64}, st_d :: Vector{Float64} ; num_th :: Int64) :: Float64

Measuring the inner product between two states and a sparse matrix. `st_d` must be of length `mat.dimd` and `st_fp` must be of length `mat.dimf`, and `st_fp` must be an adjoint. 

# Facultative argument

* `num_th :: Int64`, the number of threads. Facultative, `NumThreads` by default. 
"""
function Base.:*(st_fp :: LinearAlgebra.Adjoint{ComplexF64, Vector{ComplexF64}}, mat :: OpMat{ComplexF64}, st_d :: Vector{ComplexF64} ; num_th = NumThreads)
    ovl_ref = Ref{ComplexF64}(0)
    @ccall Libpath.__diag_MOD_scal_prod(mat.dimd :: Ref{Int64}, mat.dimf :: Ref{Int64}, mat.sym_q :: Ref{Int64}, mat.nel :: Ref{Int64}, mat.colptr :: Ref{Int64}, mat.rowid :: Ref{Int64}, mat.elval :: Ref{ComplexF64}, st_d :: Ref{ComplexF64}, st_fp' :: Ref{ComplexF64}, ovl_ref :: Ref{ComplexF64}, num_th :: Ref{Int64}) :: Nothing
    return ovl_ref[]
end
function Base.:*(st_fp :: LinearAlgebra.Adjoint{Float64, Vector{Float64}}, mat :: OpMat{Float64}, st_d :: Vector{Float64} ; num_th = NumThreads)
    ovl_ref = Ref{Float64}(0)
    @ccall Libpath.__diag_re_MOD_scal_prod_re(mat.dimd :: Ref{Int64}, mat.dimf :: Ref{Int64}, mat.sym_q :: Ref{Int64}, mat.nel :: Ref{Int64}, mat.colptr :: Ref{Int64}, mat.rowid :: Ref{Int64}, mat.elval :: Ref{Float64}, st_d :: Ref{Float64}, st_fp' :: Ref{Float64}, ovl_ref :: Ref{Float64}, num_th :: Ref{Int64}) :: Nothing
    return ovl_ref[]
end
