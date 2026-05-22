module FuzzifiEDSparseArraysExt

using SparseArrays
using LinearAlgebra
using FuzzifiED

import FuzzifiED: OpMat


"""
    SparseMatrixCSC(mat :: OpMat{ComplexF64}) :: SparseMatrixCSC{Int64,ComplexF64}
    SparseMatrixCSC(mat :: OpMat{Float64}) :: SparseMatrixCSC{Int64,Float64}

converts the `OpMat` objects to a `SparseMatrixCSC` object in the `SparseArrays` package.
"""
function SparseArrays.SparseMatrixCSC(mat :: OpMat)
    I = Int[]
    J = Int[]
    V = eltype(mat.elval)[]
    for col in 1:mat.dimd
        for ptr in mat.colptr[col]:(mat.colptr[col+1]-1)
            push!(I, mat.rowid[ptr])
            push!(J, col)
            push!(V, mat.elval[ptr])
        end
    end
    matcsc1 = sparse(I, J, V, mat.dimf, mat.dimd)
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