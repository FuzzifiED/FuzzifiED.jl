import FuzzifiED: StoreComps!, StoreComps, Laplacian, GetComponent, GetPointValue, GetIntegral, FilterComponent
export SSphereObs, GetFermionSObs, GetBosonSObs, GetFerDensitySObs, GetBosDensitySObs


"""
    SSphereObs

The mutable type `SSphereObs` stores the information of a local observable (or local operator) ``𝒪`` that can be decomposed into angular components.
```math 
    𝒪(\\Omega)=∑_{lm}𝒪_{lm}Y^{(s)}_{lm}
```

# Fields

* `s2 :: Int64` is twice the spin ``2s`` of the observable.
* `l2m :: Int64` is twice the maximal angular momentum ``2l_{\\max}`` of the components of the observable. 
* `get_comp :: Function` is a function `get_comp(l2 :: Int64, m2 :: Int64) :: STerms` that sends the component specified by a tuple of integers ``(2l,2m)`` where ``|s|\\leq l\\leq l_{\\max}, -l\\leq m\\leq l`` to a list of terms that specifies the expression of the component. 
* `stored_q :: Bool` is a boolean that specifies whether or not each component of the observable is stored.
* `comps :: Dict{Tuple{Int64, Int64}, STerms}` stores each component of the observable in the format of a dictionary whose keys are the tuples of integers ``(2l,2m)`` and values are the lists of terms that specifies the expression of the component. 
"""
mutable struct SSphereObs
    s2 :: Int64
    l2m :: Int64
    get_comp :: Function
    stored_q :: Bool 
    comps :: Dict{Tuple{Int64, Int64}, STerms}
end


"""
    SSphereObs(s2 :: Int64, l2m :: Int64, get_comp :: Function) :: SSphereObs

initialises the observable from ``2s``, ``2l_{\\max}`` and the function ``(l,m)↦𝒪_{lm}``.

# Arguments

* `s2 :: Int64` is twice the spin ``2s`` of the observable.
* `l2m :: Int64` is twice the maximal angular momentum ``2l_{\\max}`` of the components of the observable. 
* `get_comp :: Function` is a function `get_comp(l2 :: Int64, m2 :: Int64) :: STerms` that sends the component specified by a tuple of integers ``(2l,2m)`` where ``|s|\\leq s\\leq l_{\\max}, -l\\leq m\\leq l`` to a list of terms that specifies the expression of the component. 
"""
function SSphereObs(s2 :: Int64, l2m :: Int64, get_comp :: Function)
    return SSphereObs(s2, l2m, get_comp, false, Dict{Tuple{Int64, Int64}, STerms}())
end


"""
    SSphereObs(s2 :: Int64, l2m :: Int64, comps :: Dict{Tuple{Int64, Int64}, STerms}) :: SSphereObs

initialises the observable from ``2s``, ``2l_{\\max}`` and a list of ``𝒪_{lm}`` specified by a dictionary. 

# Arguments

* `s2 :: Int64` is twice the spin ``2s`` of the observable.
* `l2m :: Int64` is twice the maximal angular momentum ``2l_{\\max}`` of the components of the observable. 
* `comps :: Dict{Tuple{Int64, Int64}, STerms}` stores each component of the observable in the format of a dictionary whose keys are the tuples of integers ``(2l,2m)`` and values are the lists of terms that specifies the expression of the component. 
"""
function SSphereObs(s2 :: Int64, l2m :: Int64, cmps :: Dict{Tuple{Int64, Int64}, STerms})
    return SSphereObs(s2, l2m, (l2, m2) -> (l2 ≤ l2m && l2 ≥ abs(s2) && abs(m2) ≤ l2 && haskey(cmps, (l2, m2))) ? cmps[(l2, m2)] : STerm[], true, cmps)
end


"""
    StoreComps!(obs :: SSphereObs)
    
calculates and stores each component of the observable `obs` and replace the function in `obs` by the list of calculated components. 
"""
function StoreComps!(obs :: SSphereObs)
    if (obs.stored_q) return end
    cmps = Dict{Tuple{Int64, Int64}, STerms}()
    s2 = obs.s2
    l2m = obs.l2m
    for l2 = abs(s2) : 2 : l2m
        for m2 = -l2 : 2 : l2 
            cmps[(l2, m2)] = SimplifyTerms(obs.get_comp(l2, m2))
        end
    end
    obs.stored_q = true 
    obs.comps = cmps
    obs.get_comp = (l2, m2) -> (l2 ≤ l2m && l2 ≥ abs(s2) && abs(m2) ≤ l2 && haskey(obs.comps, (l2, m2))) ? obs.comps[(l2, m2)] : STerm[]
end


"""
    StoreComps(obs :: SSphereObs) :: SSphereObs
    
calculates and stores each component of the observable `obs` and return a new observable with the list of calculated components. 
"""
function StoreComps(obs :: SSphereObs)
    if (obs.stored_q) return obs end
    cmps = Dict{Tuple{Int64, Int64}, STerms}()
    s2 = obs.s2
    l2m = obs.l2m
    for l2 = abs(s2) : 2 : l2m
        for m2 = -l2 : 2 : l2 
            cmps[(l2, m2)] = SimplifyTerms(obs.get_comp(l2, m2))
        end
    end
    obs.stored_q = true 
    obs.comps = cmps
    return SSphereObs(s2, l2m, cmps)
end


"""
    *(fac :: Number, obs :: SSphereObs) :: SSphereObs
    *(obs :: SSphereObs, fac :: Number) :: SSphereObs
    /(obs :: SSphereObs, fac :: Number) :: SSphereObs
    -(obs :: SSphereObs) :: SSphereObs
    
enables the multiplication of an observable with a number.
"""
function Base.:*(fac :: Number, obs :: SSphereObs) 
    return SSphereObs(obs.s2, obs.l2m, (l2, m2) -> fac * obs.get_comp(l2, m2))
end
function Base.:*(obs :: SSphereObs, fac :: Number) 
    return fac * obs
end
function Base.:/(obs :: SSphereObs, fac :: Number) 
    return (1 / fac) * obs
end
function Base.:-(obs :: SSphereObs) 
    return (-1) * obs
end


"""
    +(obs1 :: SSphereObs, obs2 :: SSphereObs) :: SSphereObs
    -(obs1 :: SSphereObs, obs2 :: SSphereObs) :: SSphereObs
    
enables the addition of two observables.
"""
function Base.:+(obs1 :: SSphereObs, obs2 :: SSphereObs) 
    if (obs1.s2 ≠ obs2.s2) 
        print("Additions must have equal S")
        return 
    end
    s2 = obs1.s2
    l2m = max(obs1.l2m, obs2.l2m)
    return SSphereObs(s2, l2m, (l2, m2) -> obs1.get_comp(l2, m2) + obs2.get_comp(l2, m2))
end
function Base.:-(obs1 :: SSphereObs, obs2 :: SSphereObs)
    return obs1 + (-1) * obs2
end


"""
    adjoint(obs :: SSphereObs) :: SSphereObs
    
enables the Hermitian conjugate of a spherical observable.
```math
\\begin{aligned}
    𝒪^†(Ω)&=∑_{lm}(𝒪_{lm})^†\\bar{Y}^{(s)}_{lm}(Ω)=∑_{lm}(𝒪_{lm})^†(-1)^{s+m}Y^{(-s)}_{l,-m}(Ω)\\\\
    (𝒪^†)_{lm}&=(-1)^{s-m}(𝒪_{l,-m})^†
\\end{aligned}
```
"""
function Base.adjoint(obs :: SSphereObs)
    s2 = obs.s2
    l2m = obs.l2m
    obs1 = SSphereObs(-s2, l2m, (l2, m2) -> obs.get_comp(l2, -m2)' * (iseven((s2 - m2) ÷ 2) ? 1 : -1))
    return obs1
end


"""
    *(obs1 :: SSphereObs, obs2 :: SSphereObs) :: SSphereObs
    
enables the multiplication of two observable by making use of the composition of two monopole harmonics into one. 
"""
function Base.:*(obs1 :: SSphereObs, obs2 :: SSphereObs)
    s21 = obs1.s2 
    s22 = obs2.s2
    s2 = s21 + s22
    l2m1 = obs1.l2m
    l2m2 = obs2.l2m
    l2m = l2m1 + l2m2
    gc = ((l2, m2) -> sum(STerms[sum(STerms[sum(STerms[
            (iseven((s2 + m2) ÷ 2) ? 1 : -1) *
            sqrt((l21 + 1) * (l22 + 1) * (l2 + 1) / (4 * π)) *
            wigner3j(l21/2, l22/2, l2/2, -s21/2, -s22/2, s2/2) *
            wigner3j(l21/2, l22/2, l2/2, m21/2, (m2-m21)/2, -m2/2) *
            obs1.get_comp(l21, m21) * obs2.get_comp(l22, m2 - m21)
        for m21 = max(-l21, -l22 + m2) : 2 : min(l21, l22 + m2)])
        for l21 = max(abs(s21), abs(l2 - l22)) : 2 : min(l2m1, l2 + l22)])
        for l22 = abs(s22) : 2 : l2m2]))
    return SSphereObs(s2, l2m, gc)
end


"""
    Laplacian(obs :: SSphereObs[ ; norm_r2 :: Float64]) :: SSphereObs
    
Takes the Laplacian of an observable
```math
    (∇^2𝒪)_{lm}=-l(l+1)𝒪_{lm}
```

# Arguments

* `norm_r2 :: Float64` is the radius squared ``R^2`` used for normalisation. Facultative, `ObsNormRadSq` by default. If ``R≠1``, an extra factor ``1/R^2`` is included. 
"""
function Laplacian(obs :: SSphereObs ; norm_r2 :: Float64 = ObsNormRadSq)
    return SSphereObs(obs.s2, obs.l2m, (l2, m2) -> -l2 / 2 * (l2 / 2 + 1) * obs.get_comp(l2, m2) / norm_r2)
end


"""
    GetIntegral(obs :: SSphereObs[ ; norm_r2 :: Float64]) :: STerms

returns the uniform spatial integral ``∫\\mathrm{d}^2𝐫\\,𝒪(𝐫)`` in the format of a list of terms.

# Arguments

* `norm_r2 :: Float64` is the radius squared ``R^2`` used for normalisation. Facultative, `ObsNormRadSq` by default. If ``R≠1``, an extra factor ``R^2`` is included. 
"""
function GetIntegral(obs :: SSphereObs ; norm_r2 :: Float64 = ObsNormRadSq)
    return obs.get_comp(0, 0) * norm_r2 * √(4π)
end


"""
    GetComponent(obs :: SSphereObs, l :: Number, m :: Number) :: STerms

returns an angular component ``𝒪_{lm}`` of an observable in the format of a list of terms.
"""
function GetComponent(obs :: SSphereObs, l :: Number, m :: Number)
    return obs.get_comp(Int64(2 * l), Int64(2 * m))
end


"""
    FilterComponent(obs :: SSphereObs, flt) :: AngModes 

returns an observable object with certain modes filtered out.

# Arguments 

* `obs :: SSphereObs` is the original observable
* `flt` is the filter function whose input is the pair ``(l,m)`` and output is a logical that indicates whether this mode is chosen. E.g., if one wants to filter out the components with ``L^z=0``, one should put `(l, m) -> m == 0`.
"""
function FilterComponent(obs :: SSphereObs, flt) 
    l2m = obs.l2m
    s2 = obs.s2
    return SSphereObs(s2, l2m, (l2, m2) -> flt(l2 / 2, m2 / 2) ? obs.get_comp(l2, m2) : STerm[])
end


"""
    GetPointValue(obs :: SSphereObs, θ :: Float64, ϕ :: Float64) :: STerms

evaluates an observable at one point ``𝒪(θ,ϕ)`` in the format of a list of terms.
"""
function GetPointValue(obs :: SSphereObs, θ :: Float64, ϕ :: Float64)
    if (obs.s2 ≠ 0) 
        println("S ≠ 0 not supported")
        return
    end
    lm = obs.l2m ÷ 2
    Ylm = computeYlm(θ, ϕ, lmax = lm)
    tms = STerm[]
    for l = 0 : lm 
        for m = -l : l 
            tms += obs.get_comp(2 * l, 2 * m) * Ylm[(l,m)]
        end
    end
    return tms
end


"""
    GetFermionSObs(nm :: Int64, nf :: Int64, f :: Int64[ ; norm_r2 :: Float64]) :: SSphereObs

returns the fermion annihilation operator ``ψ_f``.

# Arguments

* `nf :: Int64` is the number of flavours.
* `nm :: Int64` is the number of orbitals.
* `f :: Int64` is the index of the flavour to be taken.
* `norm_r2 :: Float64` is the radius squared ``R^2`` used for normalisation. Facultative, `ObsNormRadSq` by default. If ``R≠1``, an extra factor ``1/R`` is included. 
"""
function GetFermionSObs(nm :: Int64, nf :: Int64, f :: Int64 ; norm_r2 :: Float64 = ObsNormRadSq)
    gc = (l2, m2) -> (l2 == nm - 1) ? STerms(1 / √norm_r2, [0, f + nf * ((m2 + nm - 1) ÷ 2)]) : STerm[]
    return SSphereObs(nm - 1, nm - 1, gc)
end


"""
    GetBosonSObs(nm :: Int64, nf :: Int64, f :: Int64[ ; norm_r2 :: Float64]) :: SSphereObs

returns the boson annihilation operator ``ϕ_f``.

# Arguments

* `nf :: Int64` is the number of flavours.
* `nm :: Int64` is the number of orbitals.
* `f :: Int64` is the index of the flavour to be taken.
* `norm_r2 :: Float64` is the radius squared ``R^2`` used for normalisation. Facultative, `ObsNormRadSq` by default. If ``R≠1``, an extra factor ``1/R`` is included. 
"""
function GetBosonSObs(nm :: Int64, nf :: Int64, f :: Int64 ; norm_r2 :: Float64 = ObsNormRadSq)
    gc = (l2, m2) -> (l2 == nm - 1) ? STerms(1 / √norm_r2, [0, -(f + nf * ((m2 + nm - 1) ÷ 2))]) : STerm[]
    return SSphereObs(nm - 1, nm - 1, gc)
end


"""
    GetFerDensitySObs(nm :: Int64, nf :: Int64[, mat :: Matrix{<:Number}][ ; norm_r2 :: Float64]) :: SSphereObs

returns the fermion density operator ``n_c=∑_{ff'}ψ^†_{f}M_{ff'}ψ_{f'}``

# Arguments

* `nf :: Int64` is the number of flavours.
* `nm :: Int64` is the number of orbitals.
* `mat :: Int64` is the matrix ``M_{ff'}``. Facultative, identity matrix ``\\mathbb{I}`` by default.
* `norm_r2 :: Float64` is the radius squared ``R^2`` used for normalisation. Facultative, `ObsNormRadSq` by default. If ``R≠1``, an extra factor ``1/R^2`` is included. 
"""
function GetFerDensitySObs(nm :: Int64, nf :: Int64, mat :: Matrix{<:Number} = Matrix{Float64}(I, nf, nf) ; norm_r2 :: Float64 = ObsNormRadSq)
    el = [ StoreComps(GetFermionSObs(nm, nf, f ; norm_r2)) for f = 1 : nf ]
    obs = SSphereObs(0, 0, Dict{Tuple{Int64, Int64}, STerms}())
    for f1 = 1 : nf, f2 = 1 : nf
        if abs(mat[f1, f2]) < 1E-13 continue end 
        obs += mat[f1, f2] * el'[f1] * el[f2]
    end
    return obs
end


"""
    GetBosDensitySObs(nm :: Int64, nf :: Int64[, mat :: Matrix{<:Number}][ ; norm_r2 :: Float64]) :: SSphereObs

returns the boson density operator ``n_c=∑_{ff'}ϕ^†_{f}M_{ff'}ϕ_{f'}``

# Arguments

* `nf :: Int64` is the number of flavours.
* `nm :: Int64` is the number of orbitals.
* `mat :: Int64` is the matrix ``M_{ff'}``. Facultative, identity matrix ``\\mathbb{I}`` by default.
* `norm_r2 :: Float64` is the radius squared ``R^2`` used for normalisation. Facultative, `ObsNormRadSq` by default. If ``R≠1``, an extra factor ``1/R^2`` is included. 
"""
function GetBosDensitySObs(nm :: Int64, nf :: Int64, mat :: Matrix{<:Number} = Matrix{Float64}(I, nf, nf) ; norm_r2 :: Float64 = ObsNormRadSq)
    el = [ StoreComps(GetBosonSObs(nm, nf, f ; norm_r2)) for f = 1 : nf ]
    obs = SSphereObs(0, 0, Dict{Tuple{Int64, Int64}, STerms}())
    for f1 = 1 : nf, f2 = 1 : nf
        if abs(mat[f1, f2]) < 1E-13 continue end 
        obs += mat[f1, f2] * el'[f1] * el[f2]
    end
    return obs
end
