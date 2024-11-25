# Fuzzifino

Fuzzifino is a module for exact diagonalisation (ED) calculation on the fuzzy sphere for systems with both bosons and fermions and similar number of boson and fermion orbitals. The usage is similar to FuzzifiED, with new types `SQNDiag`, `SQNOffd`, `SConf`, `SBasis`, `STerm` and `SOperator` defined. To use the module, include also at the start of your Julia script
```julia
using FuzzifiED.Fuzzifino
```

## Environment parameter

```@docs
Libpathino
```

## Quantum numbers

The diagonal and off-diagonal quantum numbers are implemented as
```@docs
SQNDiag
SQNOffd
```

## Configurations
```@docs
SConfs
```
It can be generated from the QNDiags.
```@docs
SConfs(nof :: Int64, nob :: Int64, nebm :: Int64, secd :: Vector{Int64}, qnd :: Vector{SQNDiag} ; num_th :: Int64 = NumThreads, disp_std :: Bool = !SilentStd)
```

## Basis
```@docs
SBasis
```
It can be generated by the following methods.
```@docs
SBasis(cfs :: SConfs, secf :: Vector{<:Number}, qnf :: Vector{SQNOffd} ; num_th = NumThreads, disp_std = !SilentStd)
SBasis(cfs :: SConfs)
```

## Term

```@docs
STerm
```
The product of terms with a number, the sum and product of terms, adjoint and particle-hole transformation are defined
```@docs
*(fac :: Number, tms :: Vector{STerm})
+(tms1 :: Vector{STerm}, tms2 :: Vector{STerm})
*(tms1 :: Vector{STerm}, tms2 :: Vector{STerm})
adjoint(tms :: Vector{STerm})
```
The terms can be simplified by 
```@docs
NormalOrder(tm :: STerm)
SimplifyTerms(tms :: Vector{STerm})
```

## Operator

```@docs
SOperator
```
It can be generated by the following methods.
```@docs
SOperator(bsd :: SBasis, bsf :: SBasis, terms :: Vector{STerm} ; red_q :: Int64 = 0, sym_q :: Int64 = 0)
```
The product of an operator on a state and the inner product of a final state, an operator and an initial state can be calculated by
```@docs
*(op :: SOperator, st_d :: Vector{ComplexF64} ; num_th = NumThreads, disp_std = !SilentStd)
*(st_fp :: LinearAlgebra.Adjoint{ComplexF64, Vector{ComplexF64}}, op :: SOperator, st_d :: Vector{ComplexF64} ; num_th = NumThreads, disp_std = !SilentStd)
```

## Sparse matrix

The OpMat can be generated from `SOperator` by the following methods.
```@docs
OpMat(op :: SOperator)
```
After the generation of sparse matrix, the diagonalisation can be condicted with FuzzifiED. 

## Test

The module is tested by a simple example

* [`test_boson.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/test_boson.jl) tests the nearest-neighbour tight-binding model $H=\sum_i(b^\dagger_ib_{i+1}+f^\dagger_if_{i+1}+\mathrm{h.c.})$. The example diagonalises the sector with the number of bosons and fermions both $N_o/2$, and even under the reflection with respect to a bond center $i\mapsto N_o+1-i$, and measures the total particle number squared $\left[\sum_i(b_i^\dagger b_i+f^\dagger_if_i)\right]^2$