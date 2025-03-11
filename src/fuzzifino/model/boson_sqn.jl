export GetNeSQNDiag, GetBosonLz2SQNDiag, GetBosonFlavSQNDiag
export GetBosonFlavPermSQNOffd, GetBosonRotySQNOffd


""" 
    GetNeSQNDiag(nof :: Int64, nob :: Int64) :: SQNDiag 

Return the SQNDiag of the number of particles, implemented as 
```julia
SQNDiag("Ne", fill(1, nof), fill(1, nob))
```
"""
GetNeSQNDiag(nof :: Int64, nob :: Int64) = SQNDiag("Ne", fill(1, nof), fill(1, nob))


""" 
    GetBosonLz2SQNDiag(nof :: Int64, nm :: Int64, nf :: Int64) :: SQNDiag 

Return the SQNDiag of twice the angular momentum ``2L_z`` for pure bosons, implemented as 
```julia
SQNDiag("Lz", fill(0, nof), collect(0 : nm * nf - 1) .÷ nf .* 2 .- (nm - 1))
```
"""
GetBosonLz2SQNDiag(nof :: Int64, nm :: Int64, nf :: Int64) = SQNDiag("Lz", fill(0, nof), collect(0 : nm * nf - 1) .÷ nf .* 2 .- (nm - 1))

""" 
    GetBosonFlavSQNDiag(nof :: Int64, nm :: Int64, nf :: Int64, qf :: Dict{Int64, Int64}[, id :: Int64 = 1, modul :: Int64 = 1]) :: QNDiag 
    GetBosonFlavSQNDiag(nof :: Int64, nm :: Int64, nf :: Int64, qf :: Vector{Int64}[, id :: Int64 = 1, modul :: Int64 = 1]) :: QNDiag 

Return the QNDiag of linear combination of number of particles in each flavour for pure bosons, 
```math
    Q = ∑_{f}q_fn_f
```
the factor ``q_f`` can either be given by a length-``N_f`` vector or a dictionary containing non-zero terms. _E.g._, for ``Q=n_{f=1}-n_{f=3}`` in a 4-flavour system, `qf = [1, 0, -1, 0]` or `qf = Dict(1 => 1, 3 => -3)`. `id` is an index to be put in the name to distinguish. For `qf` given as vector, the function is implemented as 
```julia
SQNDiag("Sz\$id", fill(0, nof), qf[collect(0 : nm * nf - 1) .% nf .+ 1], modul)
```
"""
GetBosonFlavSQNDiag(nof :: Int64, nm :: Int64, nf :: Int64, qf :: Union{Dict{Int64, Int64}, Vector{Int64}}, id :: Int64 = 1, modul :: Int64 = 1) = SQNDiag("Sz$id", fill(0, nof), FuzzifiED.DictOrVectorInt(qf, nf)[collect(0 : nm * nf - 1) .% nf .+ 1], modul)


"""
    GetBosonFlavPermSQNOffd(nof :: Int64, nm :: Int64, nf :: Int64, permf[, fac][, cyc :: Int64])
    
Return the flavour permutaiton transformation for pure bosons
```math
    𝒵: b^†_{mf}↦α_fb^†_{mπ_f}
```
# Arguments 

* `nm :: Int64` and `nf :: Int64` are the number of orbitals and the flavours.
* `permf :: Dict{Int64, Int64}`, `permf :: Vector{Vector{Int64}}` or `Vector{Int64}` gives the flavour permutation ``π_f``. It is either a vector of the cycles, a vector of the target flavours, or a dictionary of the changed elements. _E.g._, a permutation ``1↦4,2↦5,3↦3,4↦1,5↦6,6↦2`` can be expressed as `[4,5,3,1,6,2]`, `[[1,4],[2,5,6]]` or `Dict(1=>4,2=>5,4=>1,5=>6,6=>2)`. Facultative, identity by default. 
* `fac :: Dict{Int64, <: Number}` or `Vector{<: Number}` gives the factor ``α_f``. It is either a vector of all vectors, or a dictionary of all non-unity elements. Facultative, all unity by default. 
* `cyc :: Int64` is the period of the permutation. 
"""
function GetBosonFlavPermSQNOffd(
    nof :: Int64, nm :: Int64, nf :: Int64, 
    permf :: Union{Dict{Int64, Int64}, Vector{Vector{Int64}}, Vector{Int64}}, 
    fac :: Union{Dict{Int64, <: Number}, Vector{<: Number}} = Dict{Int64, ComplexF64}(), 
    cyc :: Int64 = 2
)
    permf1 = FuzzifiED.PermDictOrVector(permf, nf)
    fac1 = ComplexF64.(FuzzifiED.DictOrVectorPhase(fac, nf))
    return SQNOffd(
        collect(1 : nof),
        vcat([permf1[f] + (m - 1) * nf for f = 1 : nf, m = 1 : nm]...),
        fill(ComplexF64(1), nof),
        vcat([fac1[f] for f = 1 : nf, m = 1 : nm]...),
        cyc
    )
end
GetBosonFlavPermSQNOffd(nof :: Int64, nm :: Int64, nf :: Int64, permf :: Union{Dict{Int64, Int64}, Vector{Vector{Int64}}, Vector{Int64}}, cyc :: Int64) = GetBosonFlavPermSQNOffd(nof, nm, nf, permf, Dict{Int64, Int64}(), cyc)


"""
    GetBosonRotySQNOffd(nof :: Int64, nm :: Int64, nf :: Int64)
    
Return the ``π``-rotation with respect to the ``y``-axis for pure bosons. 
```math
    ℛ_y: b^†_{mf}↦(-)^{m+s}b^†_{-mf}
```
"""
GetBosonRotySQNOffd(nof :: Int64, nm :: Int64, nf :: Int64) = SQNOffd(
    collect(1 : nof),
    vcat([f + (nm - m) * nf for f = 1 : nf, m = 1 : nm]...), 
    fill(ComplexF64(1), nof),
    ComplexF64(-1) .^ (collect(0 : nm * nf - 1) .÷ nf)
)