# Built-in models

## Diagonal quantum numbers on fuzzy sphere

The following diagonal quantum numbers (symmetry charges) on fuzzy sphere are built in in FuzzifiED.

```@docs
GetNeQNDiag(no :: Int64)
GetLz2QNDiag(nm :: Int64, nf :: Int64)
GetFlavQNDiag(nm :: Int64, nf :: Int64, qf :: Union{Dict{Int64, Int64}, Vector{Int64}}, id :: Int64 = 1, modul :: Int64 = 1)
GetZnfChargeQNDiag(nm :: Int64, nf :: Int64)
GetPinOrbQNDiag(no :: Int64, pin_o :: Vector{Int64}, id :: Int64 = 1) 
```

## Off-diagonal quantum numbers on fuzzy sphere

The following off-diagonal quantum numbers (transformations) on fuzzy sphere are built in in FuzzifiED.

```@docs
GetParityQNOffd(nm :: Int64, nf :: Int64, permf :: Union{Dict{Int64, Int64}, Vector{Vector{Int64}}, Vector{Int64}} = Dict{Int64, Int64}(), fac :: Union{Dict{Int64, <: Number}, Vector{<: Number}} = Dict{Int64, ComplexF64}()) 
GetFlavPermQNOffd(nm :: Int64, nf :: Int64, permf :: Union{Dict{Int64, Int64}, Vector{Vector{Int64}}, Vector{Int64}}, fac :: Union{Dict{Int64, <: Number}, Vector{<: Number}} = Dict{Int64, ComplexF64}())
GetRotyQNOffd(nm :: Int64, nf :: Int64)
```

## Operators on fuzzy sphere

```@docs 
GetIntMatrix(nm :: Int64, ps_pot :: Vector{<:Number})
GetDenIntTerms(nm :: Int64, nf :: Int64, ps_pot :: Vector{<:Number}, mat_a :: Matrix{<:Number} = Matrix{Float64}(I, nf, nf), mat_b :: Matrix{<:Number} = Matrix(mat_a'))
GetDenIntTerms(nm :: Int64, nf :: Int64, ps_pot :: Vector{<:Number}, mats_a :: Vector{<:AbstractMatrix{<:Number}}, mats_b :: Vector{<:AbstractMatrix{<:Number}} = [Matrix(mat_a') for mat_a in mats_a])
GetPairIntTerms(nm :: Int64, nf :: Int64, ps_pot :: Vector{<:Number}, mat_a :: Matrix{<:Number}, mat_b :: Matrix{<:Number} = Matrix(mat_a'))
GetPolTerms(nm :: Int64, nf :: Int64, mat :: Matrix{<:Number})
GetL2Terms(nm :: Int64, nf :: Int64)
GetC2Terms(nm :: Int64, nf :: Int64, mat_gen :: Vector{<:AbstractMatrix{<:Number}})
```

## Spherical observables

FuzzifiED supports local observables on sphere that can be decomposed into angular components ``\Phi(\Omega)=\sum_{lm}\Phi_{lm}Y^{(s)}_{lm}``
```@docs
SphereObs
```
It can be initialised with the following methods 
```@docs
SphereObs(s2 :: Int64, l2m :: Int64, get_comp :: Function)
SphereObs(s2 :: Int64, l2m :: Int64, cmps :: Dict{Tuple{Int64, Int64}, Terms})
```
The following methods explicitly calculates and stores each component
```@docs
StoreComps!(obs :: SphereObs)
StoreComps(obs :: SphereObs)
```
The multiplication, addition, conjugate and Laplacian operation of an observable is supported 
```@docs
*(fac :: Number, obs :: SphereObs) 
+(obs1 :: SphereObs, obs2 :: SphereObs) 
adjoint(obs :: SphereObs)
*(obs1 :: SphereObs, obs2 :: SphereObs)
Laplacian(obs :: SphereObs)
```
Here, a parameter `ObsNormRadSq` is used to specify the normalisation of the observables
```@docs
FuzzifiED.ObsNormRadSq
```
One can evaluate an angular component, the value at one point, the uniform spatial integral of an observable, or filter a specific set of components.
```@docs
GetComponent(obs :: SphereObs, l :: Number, m :: Number)
GetPointValue(obs :: SphereObs, θ :: Float64, ϕ :: Float64)
GetIntegral(obs :: SphereObs)
FilterComponent(obs :: SphereObs, flt) 
```
Three types of operators, _viz._ electrons and density operators, and pairing operators are built-in.
```@docs
GetElectronObs(nm :: Int64, nf :: Int64, f :: Int64)
GetDensityObs(nm :: Int64, nf :: Int64, mat :: Matrix{<:Number})
GetPairingObs(nm :: Int64, nf :: Int64, mat :: Matrix{<:Number})
```

## Angular modes

An angular modes object is similar to spherical observables except that it superposes in the rule of Clebsch-Gordan coefficients and does not have the notion of locality.
```@docs
AngModes
```
It can be initialised with the following methods 
```@docs
AngModes(l2m :: Int64, get_comp :: Function)
AngModes(l2m :: Int64, cmps :: Dict{Tuple{Int64, Int64}, Terms})
```
The following methods explicitly calculates and stores each component
```@docs
StoreComps!(amd :: AngModes)
StoreComps(amd :: AngModes)
```
The multiplication, addition and conjugate of an observable is supported 
```@docs
*(fac :: Number, amd :: AngModes) 
+(obs1 :: AngModes, obs2 :: AngModes) 
adjoint(amd :: AngModes)
*(obs1 :: AngModes, obs2 :: AngModes)
```
One can take out either one or a set of components
```@docs
GetComponent(amd :: AngModes, l :: Number, m :: Number)
FilterL2(amd :: AngModes, l :: Number) 
FilterComponent(amd :: AngModes, flt) 
```
Three types of operators, _viz._ electrons and density operators, and pairing operators are built-in.
```@docs
GetElectronMod(nm :: Int64, nf :: Int64, f :: Int64)
GetPairingMod(nm :: Int64, nf :: Int64, mat :: Matrix{<:Number})
GetDensityMod(nm :: Int64, nf :: Int64, mat :: Matrix{<:Number})
```