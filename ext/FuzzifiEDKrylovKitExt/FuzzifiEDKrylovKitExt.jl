module FuzzifiEDKrylovKitExt

using KrylovKit
using FuzzifiED
import FuzzifiED: GetEigensystemKrylov
import FuzzifiED: NumThreads, SilentStd, ElementType, AbstractOperator


"""
    GetEigensystemKrylov(mat :: OpMat{ComplexF64}, nst :: Int64 ; initvec :: Vector{ComplexF64}, num_th :: Int64, disp_std :: Bool, kwargs...) :: Tuple{Vector{ComplexF64}, Matrix{ComplexF64}}
    GetEigensystemKrylov(mat :: OpMat{Float64}, nst :: Int64 ; initvec :: Vector{Float64}, num_th :: Int64, disp_std :: Bool, kwargs...) :: Tuple{Vector{Float64}, Matrix{Float64}}

This method calls the `eigsolve` from Julia `KrylovKit.jl` package instead of Arpack from Fortran to calculate the lowest eigenstates of sparse matrix. The performance should be similar. For an example, refer to [`ising_spectrum_krylov.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/ising_spectrum_krylov.jl).

# Arguments 

* `mat :: OpMat{ComplexF64}` or `mat :: OpMat{Float64}` is the matrix.
* `nst :: Int64` is the number of eigenstates to be calculated.
* `tol :: Float64` is the tolerence for the KrylovKit process. The default value is `1E-8`.
* `ncv :: Int64` is the maximum dimension of the Krylov subspace. The default value is `max(2 * nst, nst + 10)`. If `krylovdim` is also given, `ncv` will not be used.
* `initvec :: Vector{ComplexF64}` or `initvec :: Vector{Float64}` is the initial vector. Facultative, a random initialisation by default.
* `num_th :: Int64`, the number of threads. Facultative, `NumThreads` by default.
* `disp_std :: Bool`, whether or not the log shall be displayed. Facultative, `!SilentStd` by default.
* `kwargs...` is the options that will directly sent into `eigsolve`, see its [documentation](https://jutho.github.io/KrylovKit.jl/stable/man/eig/#KrylovKit.eigsolve) for detail.

# Output

* A length-`nst` array that has the same type as `mat` recording the eigenvalues, and 
* A `dimd`×`nst` matrix that has the same type as `mat` where every column records an eigenstate. 
"""
function GetEigensystemKrylov(mat :: OpMat{T}, nst :: Int64 ; tol :: Float64 = 1E-8, ncv :: Int64 = max(2 * nst, nst + 10), initvec = rand(T, mat.dimd), num_th = NumThreads, disp_std = !SilentStd, kwargs...) where T <: Union{ComplexF64,Float64}
    kwargs1 = haskey(kwargs, :krylovdim) ? kwargs : (kwargs..., krylovdim = ncv)
    eigval, eigvec, info = eigsolve(x -> *(mat, x ; num_th), initvec, nst, :SR ; tol, kwargs1...)
    if (disp_std) print(info) end
    return Vector{T}(eigval), Matrix{T}(hcat(eigvec...))
end

function GetEigensystemKrylov(op :: AbstractOperator, nst :: Int64, type :: DataType = ElementType ; tol :: Float64 = 1E-8, ncv :: Int64 = max(2 * nst, nst + 10), initvec = rand(type, op.bsd.dim), num_th = NumThreads, disp_std = !SilentStd, kwargs...)
    kwargs1 = haskey(kwargs, :krylovdim) ? kwargs : (kwargs..., krylovdim = ncv)
    eigval, eigvec, info = eigsolve(x -> *(op, x ; num_th), initvec, nst, :SR ; tol, kwargs1...)
    if (disp_std) print(info) end
    return Vector{type}(eigval), Matrix{type}(hcat(eigvec...))
end

end