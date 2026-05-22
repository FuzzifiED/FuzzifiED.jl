module FuzzifiEDSparseArraysExt

using SparseArrays
using LinearAlgebra
using FuzzifiED

import FuzzifiED: OpMat

"""
    canonicalize_csc!(colptr::Vector{Int}, rowval::Vector{Int}, nzval)

Sorts the row indices and corresponding non-zero values in-place for each column of a sparse matrix in compressed sparse column (CSC) format. 
This function ensures that the row indices within each column are in ascending order, which is expected by the `SparseMatrixCSC` constructor in the `SparseArrays` package. 
The function takes three arguments: `colptr`, which is a vector of column pointers; `rowval`, which is a vector of row indices; and `nzval`, which is a vector of non-zero values corresponding to the row indices. 
The sorting is performed for each column, and the function modifies the input vectors in-place.
"""

function canonicalize_csc!(
    colptr::Vector{Int},
    rowval::Vector{Int},
    nzval
)

    n = length(colptr) - 1

    for col in 1:n

        r = colptr[col]:(colptr[col+1]-1)

        # skip small columns
        length(r) <= 1 && continue

        p = sortperm(view(rowval, r))

        rowval[r] = rowval[r][p]
        nzval[r]  = nzval[r][p]

    end
end


"""
    SparseMatrixCSC(mat :: OpMat{ComplexF64}) :: SparseMatrixCSC{Int64,ComplexF64}
    SparseMatrixCSC(mat :: OpMat{Float64}) :: SparseMatrixCSC{Int64,Float64}

converts the `OpMat` objects to a `SparseMatrixCSC` object in the `SparseArrays` package.
"""
function SparseArrays.SparseMatrixCSC(mat :: OpMat)
    rowval = copy(mat.rowid)
    nzval = copy(mat.elval)
    canonicalize_csc!(mat.colptr, rowval, nzval)
    matcsc1 = SparseMatrixCSC(mat.dimf, mat.dimd, mat.colptr, rowval, nzval)
    if (mat.sym_q == 0) 
        return matcsc1
    elseif (mat.sym_q == 1)
        return matcsc1 + adjoint(matcsc1) - spdiagm(diag(matcsc1))
    elseif (mat.sym_q == 2)
        return matcsc1 + transpose(matcsc1) - spdiagm(diag(matcsc1))
    end
end


"""
    OpMat(matcsc :: SparseMatrixCSC{Int64,ComplexF64}) :: OpMat{ComplexF64}
    OpMat(matcsc :: SparseMatrixCSC{Int64,Float64}) :: OpMat{Float64}

converts the `SparseMatrixCSC` object in the `SparseArrays` package to an `OpMat` objects.
"""
function OpMat(matcsc :: SparseMatrixCSC)
    return OpMat{typeof(matcsc.nzval[1])}(matcsc.n, matcsc.m, 0, length(matcsc.rowval), matcsc.colptr, matcsc.rowval, matcsc.nzval)
end


"""
    Matrix(mat :: OpMat{ComplexF64}) :: Matrix{ComplexF64}
    Matrix(mat :: OpMat{Float64}) :: Matrix{Float64}

converts the `OpMat` objects to a full matrix.
"""
function LinearAlgebra.Matrix(mat :: OpMat)
    return Matrix(SparseMatrixCSC(mat))
end


end