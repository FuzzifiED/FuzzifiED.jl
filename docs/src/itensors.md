# ITensor extension

FuzzifiED also supports an ITensor extension, including convertion with the `Site` and `OpSum` objects from ITensor library and management of DMRG sweeping process. To use the extension, make sure the packages `ITensors`, `ITensorMPS` are properly installed, and include
```julia
using ITensors, ITensorMPS
```
at the heading of the Julia script. Also note that the ITensor package supports Julia at 1.10 the best.

During intialisation, the optimal configuration for parallelisation is already automatically set.
```julia
BLAS.set_num_threads(1);
NDTensors.Strided.disable_threads();
ITensors.enable_threaded_blocksparse();
```

## Format conversion

FuzzyfiED defines a new SityType "FuzzyFermion" that can be initialised from QNDiags to avoid overwriting the original "Fermion" type.
```@docs
ITensors.space( :: SiteType"FuzzyFermion"; o :: Int, qnd :: Vector{QNDiag})
```

The `Sites` objects in ITensor can be converted to a QNDiags and Confs with the QNs extracted. 
```@docs
GetQNDiags(sites :: Vector{Index{Vector{Pair{QN, Int64}}}})
Confs(sites :: Vector{Index{Vector{Pair{QN, Int64}}}}, sec_qn :: QN)
```
Conversely, the  `Sites` objects in ITensors can also be generated from a set of diagonal quantum numbers 
```@docs
GetSites(qnd :: Vector{QNDiag})
```
If the number of quantum numbers are too many, it can be truncated by 
```@docs
TruncateQNDiag(qnd :: Vector{QNDiag} ; trunc_lth :: Int64 = 3, trunc_wt :: Vector{Int64} = [ 10 ^ (i - trunc_lth) for i = trunc_lth : length(qnd)]) 
```

The `OpSum` objects in ITensor can be converted with the collection of `Term`'s
```@docs
Terms(opsum :: Sum{Scaled{ComplexF64, Prod{Op}}})
OpSum(tms :: Terms)
```

## Easy sweep

This tool kit facilitates the management of DMRG process. It automatically records the intermediate results and recover these results if a job is stopped and run again on HPC. It also manages the gradual increase of maximal bond dimensions and the determination of convergence by the criteria of energy. This extension required the packages `HDF5`. We also recommand using the package [`ITensorMPOConstruction`](https://github.com/ITensor/ITensorMPOConstruction.jl) for the generation of Hamiltonian MPO, which can be installed by 
```julia
using Pkg ; Pkg.add(url="https://github.com/ITensor/ITensorMPOConstruction.jl.git")
```

```@docs
EasySweep(id :: String, hmt :: MPO, st00 :: MPS ; path :: String = "./", dim_list :: Vector{Int64} = [1000,2000,3000,4000,5000,6000], proj :: Vector{String} = String[], e_tol1 :: Float64 = 1E-6, e_tol :: Float64 = 1E-7, cutoff :: Vector{Float64} = [1E-9], maxdim0 :: Vector{Int64} = [10,20,50,100,200,500], noise0 :: Vector{Float64} = [1E-4,3E-5,1E-5,3E-6,1E-6,3E-7], noise :: Vector{Float64} = [1E-6,2E-7,5E-8,1E-8,0], nsweeps :: Int64 = 10, weight :: Float64 = 100.0, observer :: AbstractObserver = EasySweepObserver(e_tol1))
SweepOne(id :: String, hmt :: MPO, st0 :: MPS, dim1 :: Int64 ; path :: String = "./", cutoff :: Vector{Float64} = [1E-9], maxdim :: Vector{Int64} = [dim1], nsweeps :: Int64 = 10, noise :: Vector{Float64} = [1E-6,1E-7,0], proj :: Vector{String} = String[], e_tol :: Float64 = 1E-6, weight :: Float64 = 100.0, observer :: AbstractObserver = EasySweepObserver(e_tol))
GetMPOSites(id :: String, tms :: Union{Terms, Sum{Scaled{ComplexF64, Prod{Op}}}}, qnd :: Vector{QNDiag} ; path :: String = "./", mpo_method :: Function = MPO)
GetMPO(id :: String, tms :: Union{Terms, Sum{Scaled{ComplexF64, Prod{Op}}}}, sites :: Vector{<:Index} ; path :: String = "./", mpo_method :: Function = MPO)
```

## Modified version of ITensor

We have forked `ITensors` and made some modifications to better suit our need. To install the modified packages, please use 

```julia
using Pkg 
Pkg.add(url="https://github.com/FuzzifiED/ITensors.jl.git")
Pkg.add(url="https://github.com/FuzzifiED/ITensorMPS.jl.git")
Pkg.add(url="https://github.com/FuzzifiED/ITensorMPOConstruction.jl.git")
```

The modifications made include

* Modify `ITensors` to allow up to 10 quantum numbers ; 
* Modify `ITensorMPS` to allow `write_when_maxdim_exceeds` in DMRG in the presence of projection matrices ;
* Modify `ITensorMPOConstruction` to allow building MPOs that has symmetry flux ;
* Modify `ITensorMPOConstruction` to allow compatibility with the newest version of `ITensors`.

Please be warned that the robustness and the backward compatibility of these modifications are not warranted. 