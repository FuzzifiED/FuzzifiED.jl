export AngModes, GetElectronMod, GetPairingMod, GetDensityMod, FilterL2, ContractMod


"""
    AngModes
    
The mutable type `AngModes` stores angular momentum components of an operator on the sphere ``Φ_{lm}`` and superposes in the rule of Clebsch-Gordan coefficients. The usage is similar to the spherical observables, except that SphereObs superposes in the rule of spherical harmonics and has the notion of locality

# Fields

* `l2m :: Int64` is twice the maximal angular momentum ``2l_{\\max}`` of the components of the modes object. 
* `get_comp :: Function` is a function `get_comp(l2 :: Int64, m2 :: Int64) :: Terms` that sends the component specified by a tuple of integers ``(2l,2m)`` where ``|s|\\leq l\\leq l_{\\max}, -l\\leq m\\leq l`` to a list of terms that specifies the expression of the component. 
* `stored_q :: Bool` is a boolean that specifies whether or not each component of the modes object is stored.
* `comps :: Dict{Tuple{Int64, Int64}, Terms}` stores each component of the modes object in the format of a dictionary whose keys are the tuples of integers ``(2l,2m)`` and values are the lists of terms that specifies the expression of the component. 
"""
mutable struct AngModes
    l2m :: Int64
    get_comp :: Function
    stored_q :: Bool 
    comps :: Dict{Tuple{Int64, Int64}, Terms}
end


"""
    AngModes(l2m :: Int64, get_comp :: Function) :: AngModes

initialises the modes object from ``2l_{\\max}`` and the function ``(l,m)↦\\Phi_{lm}``

# Arguments

* `l2m :: Int64` is twice the maximal angular momentum ``2l_{\\max}`` of the components of the modes object. 
* `get_comp :: Function` is a function `get_comp(l2 :: Int64, m2 :: Int64) :: Terms` that sends the component specified by a tuple of integers ``(2l,2m)`` where ``|s|\\leq s\\leq l_{\\max}, -l\\leq m\\leq l`` to a list of terms that specifies the expression of the component. 
"""
function AngModes(l2m :: Int64, get_comp :: Function)
    return AngModes(l2m, get_comp, false, Dict{Tuple{Int64, Int64}, Terms}())
end


"""
    AngModes(l2m :: Int64, get_comp :: Function) :: AngModes

initialises the modes object from ``2l_{\\max}`` and a list of ``\\Phi_{lm}`` specified by a dictionary. 

# Arguments

* `l2m :: Int64` is twice the maximal angular momentum ``2l_{\\max}`` of the components of the modes object. 
* `comps :: Dict{Tuple{Int64, Int64}, Terms}` stores each component of the modes object in the format of a dictionary whose keys are the tuples of integers ``(2l,2m)`` and values are the lists of terms that specifies the expression of the component. 
"""
function AngModes(l2m :: Int64, cmps :: Dict{Tuple{Int64, Int64}, Terms})
    return AngModes(l2m, (l2, m2) -> (l2 ≤ l2m && abs(m2) ≤ l2 && haskey(cmps, (l2, m2))) ? cmps[(l2, m2)] : Term[], true, cmps)
end


"""
    StoreComps!(amd :: AngModes)
    
calculates and stores each component of the modes object `amd` and replace the function in `amd` by the list of calculated components. 
"""
function StoreComps!(amd :: AngModes)
    if (amd.stored_q) return end
    cmps = Dict{Tuple{Int64, Int64}, Terms}()
    l2m = amd.l2m
    for l2 = l2m % 2 : 2 : l2m
        for m2 = -l2 : 2 : l2 
            tms = SimplifyTerms(amd.get_comp(l2, m2))
            if (length(tms) ≠ 0) cmps[(l2, m2)] = tms end
        end
    end
    amd.stored_q = true 
    amd.comps = cmps
    amd.get_comp = (l2, m2) -> (l2 ≤ l2m && abs(m2) ≤ l2 && haskey(amd.comps, (l2, m2))) ? amd.comps[(l2, m2)] : Term[]
end


"""
    StoreComps(amd :: AngModes) :: AngModes
    
calculates and stores each component of the modes object `amd` and return a new modes object with the list of calculated components. 
"""
function StoreComps(amd :: AngModes)
    if (amd.stored_q) return amd end
    cmps = Dict{Tuple{Int64, Int64}, Terms}()
    l2m = amd.l2m
    for l2 = l2m % 2 : 2 : l2m
        for m2 = -l2 : 2 : l2 
            tms = SimplifyTerms(amd.get_comp(l2, m2))
            if (length(tms) ≠ 0) cmps[(l2, m2)] = tms end
        end
    end
    amd.stored_q = true 
    amd.comps = cmps
    return AngModes(l2m, cmps)
end


"""
    *(fac :: Number, amd :: AngModes) :: AngModes
    *(amd :: AngModes, fac :: Number) :: AngModes
    /(amd :: AngModes, fac :: Number) :: AngModes
    -(amd :: AngModes) :: AngModes
    
enables the multiplication of a mode with a number.
"""
function Base.:*(fac :: Number, amd :: AngModes) 
    return AngModes(amd.l2m, (l2, m2) -> fac * amd.get_comp(l2, m2))
end
function Base.:*(amd :: AngModes, fac :: Number) 
    return fac * amd
end
function Base.:/(amd :: AngModes, fac :: Number) 
    return (1 / fac) * amd
end
function Base.:-(amd :: AngModes) 
    return (-1) * amd
end


"""
    +(amd1 :: AngModes, amd2 :: AngModes) :: AngModes
    -(amd1 :: AngModes, amd2 :: AngModes) :: AngModes
    
enables the addition of two modes.
"""
function Base.:+(amd1 :: AngModes, amd2 :: AngModes) 
    l2m = max(amd1.l2m, amd2.l2m)
    return AngModes(l2m, (l2, m2) -> amd1.get_comp(l2, m2) + amd2.get_comp(l2, m2))
end
function Base.:-(amd1 :: AngModes, amd2 :: AngModes)
    return amd1 + (-1) * amd2
end


"""
    adjoint(amd :: AngModes) :: AngModes
    
enables the Hermitian conjugate of a spherical mode.
```math
\\begin{aligned}
    (Φ^†)_{lm}&=(-1)^{l+m}(Φ_{l,-m})^†
\\end{aligned}
```
"""
function Base.adjoint(amd :: AngModes)
    l2m = amd.l2m
    amd1 = AngModes(l2m, (l2, m2) -> (iseven((l2 + m2) ÷ 2) ? 1 : -1) * amd.get_comp(l2, -m2)')
    return amd1
end


"""
    *(amd1 :: AngModes, amd2 :: AngModes) :: AngModes
    
enables the multiplication of two modes in the rule of CG coefficients. 
```math 
    Φ_{lm}=∑_{l_1l_2m_1m_2}δ_{m,m_1+m_2}⟨l_1m_1,l_2m_2|lm⟩Φ_{l_1m_1}Φ_{l_2m_2}
```
"""
function Base.:*(amd1 :: AngModes, amd2 :: AngModes)
    l2m1 = amd1.l2m
    l2m2 = amd2.l2m
    l2m = l2m1 + l2m2
    gc = ((l2, m2) -> sum(Terms[sum(Terms[sum(Terms[
            (iseven((-l21 + l22 + m2) ÷ 2) ? 1 : -1) *
            sqrt(l2 + 1) *
            wigner3j(l21/2, l22/2, l2/2, m21/2, (m2-m21)/2, -m2/2) *
            amd1.get_comp(l21, m21) * amd2.get_comp(l22, m2 - m21)
        for m21 = max(-l21, -l22 + m2) : 2 : min(l21, l22 + m2)])
        for l21 = max(l2m1 % 2, abs(l2 - l22)) : 2 : min(l2m1, l2 + l22)])
        for l22 = l2m2 % 2 : 2 : l2m2]))
    return AngModes(l2m, gc)
end


"""
    GetComponent(amd :: AngModes, l :: Number, m :: Number) :: Terms

returns an angular component ``Φ_{lm}`` of a modes object in the format of a list of terms.
"""
function GetComponent(amd :: AngModes, l :: Number, m :: Number)
    return amd.get_comp(Int64(2 * l), Int64(2 * m))
end


"""
    ContractMod(amd1 :: AngModes, amd2 :: amd2, comps :: Dict) :: Terms
    ContractMod(amd1 :: AngModes, amd2 :: AngModes, l0 :: Number) :: Terms

Return the contraction of two angular modes 
```math 
    ∑_{l}U_l∑_{m=-l}^l(-1)^{l+m}Φ_{1,l}Φ_{2,l(-m)}
```
where the list of ``U_l`` is given by a dictionary, or
```math 
    U_{l₀}∑_{m=-l₀}^{l₀}(-1)^{l₀+m}Φ_{1,l₀m}Φ_{2,l₀(-m)}.
```
"""
function ContractMod(amd1 :: AngModes, amd2 :: AngModes, comps :: Dict)
    return SimplifyTerms(sum([ U * amd1.get_comp(Int64(2 * l), m2) * amd2.get_comp(Int64(2 * l), -m2) * (iseven((Int64(2 * l) + m2) ÷ 2) ? 1 : -1) for (l, U) in comps for m2 = -Int64(2 * l) : 2 : Int64(2 * l)]))
end
ContractMod(amd1 :: AngModes, amd2 :: AngModes, l0 :: Number) = ContractMod(amd1, amd2, Dict([l0 => 1]))

"""
    FilterComponent(amd :: AngModes, flt) :: AngModes 

returns an angular modes object with certain modes filtered out.

# Arguments 

* `amd :: AngModes` is the original angular modes
* `flt` is the filter function whose input is the pair ``(l,m)`` and output is a logical that indicates whether this mode is chosen. E.g., if one wants to filter out the modes with angular momentum `l0`, one should put `(l, m) -> l == l0`.
"""
function FilterComponent(amd :: AngModes, flt) 
    l2m = amd.l2m
    return AngModes(l2m, (l2, m2) -> flt(l2 / 2, m2 / 2) ? amd.get_comp(l2, m2) : Term[])
end


"""
    FilterL2(amd :: AngModes, l :: Number) :: AngModes 

returns an angular modes object with modes of a certain total angular momentum filtered out.
"""
function FilterL2(amd :: AngModes, l :: Number) 
    l2m = amd.l2m
    return AngModes(l2m, (l2, m2) -> l2 == Int64(l * 2) ? amd.get_comp(l2, m2) : Term[])
end


"""
    GetElectronMod(nm :: Int64, nf :: Int64, f :: Int64) :: AngModes

returns the modes of electron annihilation operator ``c_m``, with angular momentum ``s=(N_m-1)/2``

# Arguments

* `nf :: Int64` is the number of flavours.
* `nm :: Int64` is the number of orbitals.
* `f :: Int64` is the index of the flavour to be taken.
"""
function GetElectronMod(nm :: Int64, nf :: Int64, f :: Int64)
    gc = (l2, m2) -> (l2 == nm - 1) ? Terms(1.0, [0, f + nf * ((m2 + nm - 1) ÷ 2)]) : Term[]
    return AngModes(nm - 1, gc)
end


"""
    GetPairingMod(nm :: Int64, nf :: Int64, mat :: Matrix{<:Number}) :: AngModes

returns the modes of two electrons superposed in the rule of CG coefficients. 
```math 
    Δ_{lm}=∑_{m_1m_2}δ_{m,m_1+m_2}⟨sm_1,sm_2|lm⟩c_{a,m_1}M_{ab}c_{b,m_2}
```

# Arguments

* `nf :: Int64` is the number of flavours.
* `nm :: Int64` is the number of orbitals.
* `mat :: Int64` is the matrix ``M_{ff'}``.
"""
function GetPairingMod(nm :: Int64, nf :: Int64, mat :: Matrix{<:Number})
    el = [ StoreComps(GetElectronMod(nm, nf, f)) for f = 1 : nf ]
    amd = AngModes(2 * (nm - 1), Dict{Tuple{Int64, Int64}, Terms}())
    for f1 = 1 : nf, f2 = 1 : nf
        if abs(mat[f1, f2]) < 1E-13 continue end 
        amd += mat[f1, f2] * el[f1] * el[f2]
    end
    return amd
end


"""
    GetDensityMod(nm :: Int64, nf :: Int64, mat :: Matrix{<:Number}) :: AngModes

returns the modes of electron creation and annihilation superposed in the rule of CG coefficients. 
```math 
    n_{lm}=∑_{m_1m_2}δ_{m,-m_1+m_2}(-1)^{s+m_1}⟨s(-m_1),sm_2|lm⟩c^†_{m_1}c_{m_2}
```

# Arguments

* `nf :: Int64` is the number of flavours.
* `nm :: Int64` is the number of orbitals.
* `mat :: Int64` is the matrix ``M_{ff'}``. Facultative, identity matrix ``\\mathbb{I}`` by default.
"""
function GetDensityMod(nm :: Int64, nf :: Int64, mat :: Matrix{<:Number})
    el = [ StoreComps(GetElectronMod(nm, nf, f)) for f = 1 : nf ]
    amd = AngModes(2 * (nm - 1), Dict{Tuple{Int64, Int64}, Terms}())
    for f1 = 1 : nf, f2 = 1 : nf
        if abs(mat[f1, f2]) < 1E-13 continue end 
        amd += mat[f1, f2] * el'[f1] * el[f2]
    end
    return amd
end
