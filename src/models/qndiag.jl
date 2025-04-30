export GetNeQNDiag, GetLz2QNDiag, GetFlavQNDiag, GetZnfChargeQNDiag, GetPinOrbQNDiag


""" 
    GetNeQNDiag(no :: Int64[, modul :: Int64 = 1]) :: QNDiag 

Return the QNDiag of the number of electrons, implemented as 
```julia
QNDiag("Ne", fill(1, no), modul)
```
If `modul = 2`, parity is returned.
"""
GetNeQNDiag(no :: Int64, modul :: Int64 = 1) = QNDiag("Ne", fill(1, no), modul)


""" 
    GetLz2QNDiag(nm :: Int64, nf :: Int64) :: QNDiag 

Return the QNDiag of twice the angular momentum ``2L_z``, implemented as 
```julia
QNDiag("Lz", collect(0 : nm * nf - 1) .÷ nf .* 2 .- (nm - 1))
```
"""
GetLz2QNDiag(nm :: Int64, nf :: Int64) = QNDiag("Lz", collect(0 : nm * nf - 1) .÷ nf .* 2 .- (nm - 1))

function DictOrVectorInt(qf_dict :: Dict{Int64, Int64}, nf :: Int64)
    qf = zeros(Int64, nf)
    for (f, q) in qf_dict 
        qf[f] = q
    end
    return qf
end
DictOrVectorInt(qf :: Vector{Int64}, nf :: Int64) = qf


""" 
    GetFlavQNDiag(nm :: Int64, nf :: Int64, qf :: Dict{Int64, Int64}[, id :: Int64 = 1, modul :: Int64 = 1]) :: QNDiag 
    GetFlavQNDiag(nm :: Int64, nf :: Int64, qf :: Vector{Int64}[, id :: Int64 = 1, modul :: Int64 = 1]) :: QNDiag 

Return the QNDiag of linear combination of number of electrons in each flavour, 
```math
    Q = ∑_{f}q_fn_f
```
the factor ``q_f`` can either be given by a length-``N_f`` vector or a dictionary containing non-zero terms. _E.g._, for ``Q=n_{f=1}-n_{f=3}`` in a 4-flavour system, `qf = [1, 0, -1, 0]` or `qf = Dict(1 => 1, 3 => -3)`. `id` is an index to be put in the name to distinguish. For `qf` given as vector, the function is implemented as 
```julia
QNDiag("Sz\$id", qf[collect(0 : nm * nf - 1) .% nf .+ 1], modul)
```
"""
GetFlavQNDiag(nm :: Int64, nf :: Int64, qf :: Union{Dict{Int64, Int64}, Vector{Int64}}, id :: Int64 = 1, modul :: Int64 = 1) = QNDiag("Sz$id", DictOrVectorInt(qf, nf)[collect(0 : nm * nf - 1) .% nf .+ 1], modul)


""" 
    GetZnfChargeQNDiag(nm :: Int64, nf :: Int64) :: QNDiag 

Return the QNDiag of a ``ℤ_{N_f}``-charge, 
```math
    Q = ∑_{f=0}^{N_f-1}fn_f\\mod N_f
```
implemented as 
```julia
QNDiag("Q_Z\$nf", collect(0 : nm * nf - 1) .% nf, nf)
```
"""
GetZnfChargeQNDiag(nm :: Int64, nf :: Int64) = QNDiag("QZ$nf", collect(0 : nm * nf - 1) .% nf, nf)


""" 
    GetPinOrbQNDiag(no :: Int64, pin_o :: Vector{Int64}[, id :: Int64 = 1]) :: QNDiag 

Return the QNDiag of the number of electrons in the subset of sites `pin_o`, implemented as
```julia
QNDiag("Npin\$i", [ o in pin_o ? 1 : 0 for o = 1 : no])
```
This QNDiag is useful in pinning defects, where certain subset of sites need to be set empty or filled. To empty the sites, set this QNDiag to 0 ; to fill the sites, set this QNDiag to `length(pin_o)`. `id` is an index to be put in the name to distinguish. 
"""
GetPinOrbQNDiag(no :: Int64, pin_o :: Vector{Int64}, id :: Int64 = 1) = QNDiag("Npin$id", [ o in pin_o ? 1 : 0 for o = 1 : no])