# Other Extensions

Apart from ITensor extension, FuzzifiED also provides other extensions, _viz._ HDF5 extension, KrylovKit extension, CUDA extension and SparseArrays extension. 

## HDF5 Extension 

The HDF5 extension supports writing the types `Confs`, `Basis`, `Terms`, `Operator`, `OpMat{ComplexF64}` and `OpMat{Float64}` into HDF5 files and reading them from groups and subgroups in HDF5 format. This extension requires the packages `HDF5`. To use this extension, include at the heading 
```julia
using HDF5
```
A typical file operation process looks like
```julia
f = h5open(file_name, "cw")
# include the file name as a string 
# Modes : "cw" for write and "r" for read
...
close(f)
```

To write, include in the middle 
```julia
write(f, group_name :: String, cfs :: Confs)
write(f, group_name :: String, bs  :: Basis)
write(f, group_name :: String, tms :: Terms)
write(f, group_name :: String, op  :: Operator)
write(f, group_name :: String, mat :: OpMat{ComplexF64})
write(f, group_name :: String, mat :: OpMat{Float64})
```
To read, include in the middle 
```julia
cfs = read(f, group_name :: String, Confs)
bs  = read(f, group_name :: String, Basis)
tms = read(f, group_name :: String, Terms)
op  = read(f, group_name :: String, Operator)
mat = read(f, group_name :: String, OpMat{ComplexF64})
mat = read(f, group_name :: String, OpMat{Float64})
```

## SparseArrays Extension 

The SparseArrays extension supports the conversion between `OpMat` in FuzzifiED and the `SparseMatrixCSC` format. This extension requires the packages `SparseArrays`. To use this extension, include at the heading 
```julia
using SparseArrays
```
```@docs
SparseMatrixCSC(mat :: OpMat)
OpMat(matcsc :: SparseMatrixCSC)
```

## KrylovKit Extension

The KrylovKit extension supports an interface diagonalising the sparse matrix using the `KrylovKit` pakage in Julia. This extension requires the packages `KrylovKit`. To use this extension, include at the heading 
```julia
using KrylovKit
```
Besides Arpack in Fortran, we also provide an interface calling KrylovKit in Julia.
```@docs
GetEigensystemKrylov(mat :: OpMat{ComplexF64}, nst :: Int64)
```

### Modified Version of KrylovKit

We have forked `KrylovKit` and made some modifications to better suit our need. To install the modified packages, please use 

```julia
using Pkg 
Pkg.add(url="https://github.com/FuzzifiED/KrylovKit.jl.git")
```

We add an algorithm `LockedLanczos`. It is a hard-locking (deflated) variant of the Lanczos algorithm that is useful when large number of eigenvectors (_e. g._ in the order of 1000) are required, for use in `eigsolve` with a real symmetric or complex Hermitian linear operator. Whenever a Ritz pair has converged to within a tolerence, it is locked, _i. e._ removed from the active Krylov subspace, but retained in the set against which every subsequent Krylov vector is orthogonalized.

## CUDA Extension

The CUDA extension supports the conversion between `OpMat` in FuzzifiED and the `CuSparseMatrixCSC` in `CUDA.CUSPARSE` as well as the acceleration of diagonalisation on GPU. This extension requires the packages `CUDA`, `KrylovKit` and `SparseArrays`. To use this extension, include at the heading 
```julia
using CUDA, KrylovKit, SparseArrays
```
```@docs
CUSPARSE.CuSparseMatrixCSC(mat :: OpMat{ComplexF64})
GetEigensystemCuda(mat :: OpMat{ComplexF64}, nst :: Int64)
```