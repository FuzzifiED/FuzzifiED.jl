# Fuzzy Manifolds

FuzzyManifolds is a module to support calculations on other geometry regularised by lowest Landau level besides sphere. The supported geometries include torus $T^2$. Other geometries (_e. g._, disk) can be added upon request. To use the module, include also at the start of your Julia script
```julia
using FuzzifiED.FuzzyManifolds
```

## Torus 

```@docs
GetTorusLz2QNDiag(nm :: Int64, nf :: Int64) 
GetTorusTranslQNOffd(nm :: Int64, nf :: Int64) 
GetTorusIntMatrix(nm :: Int64, lx :: Number, ps_pot :: Vector{<:Number})
GetTorusDenIntTerms(nm :: Int64, nf :: Int64, lx :: Number, ps_pot :: Vector{<:Number}, mat_a :: Matrix{<:Number} = Matrix{Float64}(I, nf, nf), mat_b :: Matrix{<:Number} = Matrix(mat_a'))
GetTorusPairIntTerms(nm :: Int64, nf :: Int64, lx :: Number, ps_pot :: Vector{<:Number}, mat_a :: Matrix{<:Number}, mat_b :: Matrix{<:Number} = Matrix(mat_a'))
```

## Related Examples 

* [`circle_ising.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/circle_ising.jl) calculates the spectrum of 2d Ising CFT on a fuzzy thin torus. This example reproduces Figure 4 and Tables I—III in [Han 2025](@ref Han2025)
