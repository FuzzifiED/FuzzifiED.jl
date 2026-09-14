export CoupleDecomp, CoupleDecomps
export RecoupleAngMom, ConvPsPot, ContactCouple, SingleSegCouple, InsertSegment
import FuzzifiED: AngModes


"""
    CoupleDecomp

The mutable type `CoupleDecomp` records an angular-moemntum channel of a coupling — a direct product of spherical-symmetric actions on each part, coupled to a definite total angular momentum. A full operator is represented as a `CoupleDecomps`, _i. e._ a sum of such channels. _E. g._ for bipartite and tri-partite systems, a channel may take the form
```math
\\begin{aligned}
    [𝒪]_{(l_1l_2)l}&=\\Big[[𝒪_1]_{l_1}⊗[𝒪_2]_{l_2}\\Big]_l&𝒪_{lm}&=[𝒪_1]_{l_1l_1}[𝒪_2]_{l_2m_2}⟨l_1m_1,l_2m_2|lm⟩\\\\
    [𝒪]_{((l_1l_2)l_{12}l_3)l}&=\\Big[\\big[[𝒪_1]_{l_1}[𝒪_2]_{l_2}\\big]_{l_{12}}[𝒪_3]_{l_3}\\Big]_l&𝒪_{lm}&=[𝒪_1]_{l_1l_1}[𝒪_2]_{l_2m_2}[𝒪_3]_{l_3m_3}⟨l_1m_1,l_2m_2|l_{12}m_{12}⟩⟨l_{12}m_{12},l_3m_3|lm⟩
\\end{aligned}
```

# Fields

* `amd :: Vector{Union{AngModes, SAngModes, Symbol}}` records, for each segment ``p``, its angular modes ``[Φ_p]_{lm}`` — an `AngModes` on a fermionic segment, an `SAngModes` on a bosonic one, or `:Identity` on an inert segment.
* `ch :: Matrix{Int64}` is the coupling channel in the form of a ``2×N_p`` matrix.
* `coeff :: ComplexF64` is the coefficient of the channel.
* `sec :: Matrix{Int64}` records the change of quantum numbers that the term induces : `sec[iqn, p]` is the shift of the `iqn`-th diagonal quantum number on part ``p``.
"""
mutable struct CoupleDecomp
    amd :: Vector{<:Union{AngModes, SAngModes, Symbol}}
    ch :: Matrix{Int64}
    coeff :: ComplexF64
    sec :: Matrix{Int64} # sec[iqn, p]
end
function CoupleDecomp(amd :: Vector{<:Union{AngModes, SAngModes, Symbol}}, ch :: Matrix{Int64}, sec :: Matrix{Int64})
    return CoupleDecomp(amd, ch, 1.0 + 0.0im, sec)
end


"""
    CoupleDecomps = Vector{CoupleDecomp}

Alias for `Vector{CoupleDecomp}`, representing a full operator as a list of single-channel [CoupleDecomp](@ref CoupleDecomp)s.
"""
const CoupleDecomps = Vector{CoupleDecomp}



"""
    CoupleDecomps(amd :: Vector, ch :: Vector{Matrix{Int64}}, coeff :: Vector{<:Number}, sec :: Matrix{Int64}) :: CoupleDecomps

constructs a `CoupleDecomps` — one single-channel `CoupleDecomp` per coupling channel `ch[i]` with coefficient `coeff[i]` — all sharing the operators `amd` (each entry an `AngModes`, `SAngModes` or `:Identity`) and the quantum number shift `sec`. It consumes the channels and coefficients returned by [ConvPsPot](@ref ConvPsPot).
"""
function CoupleDecomps(amd :: Vector, ch :: Vector{Matrix{Int64}}, coeff :: Vector{<:Number}, sec :: Matrix{Int64} ; eltype = FuzzifiED.ElementType)
    amd1 = Union{AngModes, SAngModes, Symbol}[]
    ph = 1.0 + 0.0im
    for amdi in amd
        if (amdi === :Identity)
            push!(amd1, amdi)
            continue
        end
        amdi1 = StoreComps(deepcopy(amdi))
        if (eltype == Float64)
            coeff1 = collect(amdi1.comps)[1][2][1].coeff
            if (abs(coeff1.re / coeff1.im) < 1E-4)
                amdi1 *= 1.0im 
                ph *= -1.0im
            end
        end
        push!(amd1, amdi1)
    end
    return CoupleDecomp[ CoupleDecomp(amd1, ch[i], ph * coeff[i], sec) for i in eachindex(ch) ]
end


Base.one(CoupleDecomps, np :: Int64, nsec :: Int64) = CoupleDecomps([ :Identity for _ = 1 : np], [zeros(Int64, 2, np)], [ComplexF64(1)], zeros(Int64, nsec, np))

"""
    cpd1 + cpd2 :: CoupleDecomps
    cpd1 - cpd2 :: CoupleDecomps
    -cpd :: CoupleDecomps
    fac * cpd :: CoupleDecomps

enable the linear combination of coupling decompositions. 
"""
function Base.:+(cpd1 :: CoupleDecomps, cpd2 :: CoupleDecomps)
    return [ cpd1 ; cpd2 ]
end
function Base.:*(fac :: Number, cpd :: CoupleDecomp)
    return CoupleDecomp(cpd.amd, cpd.ch, fac * cpd.coeff, cpd.sec)
end
function Base.:*(fac :: Number, cpd :: CoupleDecomps)
    return fac .* cpd
end
function Base.:*(cpd :: CoupleDecomps, fac :: Number)
    return fac .* cpd
end
function Base.:/(cpd :: CoupleDecomps, fac :: Number)
    return (1/fac) .* cpd
end
function Base.:-(cpd :: CoupleDecomps)
    return (-1) * cpd
end
function Base.:-(cpd1 :: CoupleDecomps, cpd2 :: CoupleDecomps)
    return cpd1 + (-1) * cpd2
end
function Base.:+(cpd1 :: CoupleDecomps, cpd2 :: Vararg{CoupleDecomps})
    return cpd1 + +(cpd2...)
end


"""
    RecoupleAngMom(s1 :: Number, s2 :: Number, s3 :: Number, s4 :: Number, ps_pot0 :: Dict) :: Dict

re-couples the pseudo-potentials of a four-fermion term ``c^†_1c^†_2c_3c_4`` from the pairing channel ``(12)(34)`` into the density channel ``(14)(23)``  using Wigner's ``6j``-symbol. 

# Arguments

* `s1, s2, s3, s4 :: Number` are the spins of the four particles.
* `ps_pot0 :: Dict` is a dictionary mapping a density rank ``l`` to its coefficient ``V_l``

# Output

`ps_pot1 :: Dict` mapping the angular momentum ``l`` to the pseudoo-potential ``W_l``
"""
function RecoupleAngMom(s1 :: Number, s2 :: Number, s3 :: Number, s4 :: Number, ps_pot0 :: Dict)
    ps_pot1 = Dict{Float64, ComplexF64}()
    for j in max(abs(s1 - s4), abs(s2 - s3)) : min(s1 + s4, s2 + s3)
        W = 0.0
        for (l, V) in ps_pot0
            (abs(s1 - s2) ≤ l ≤ s1 + s2) || continue
            (abs(s3 - s4) ≤ l ≤ s3 + s4) || continue
            sixj = Float64(wigner6j(s1, s2, l, s3, s4, j))
            sixj == 0.0 && continue
            sixj *= (1.0im) ^ Int(-2l)
            W += (2l + 1) * sixj * V
        end
        ps_pot1[j] = W
    end
    return ps_pot1
end
"""
    RecoupleAngMom(s :: Number, ps_pot :: Vector{<:Number}) :: Dict

re-couples the pseudo-potentials of a four-fermion term ``c^†_1c^†_2c_3c_4`` from the pairing channel ``(12)(34)`` into the density channel ``(14)(23)``  using Wigner's ``6j``-symbol. 

# Arguments

* `s` is the single-particle spin, and 
* `ps_pot` is a vector of pseudo-potentials from `l=2s` downwards. 

# Output

`ps_pot1 :: Dict` mapping the angular momentum ``l`` to the pseudoo-potential ``W_l``
"""
function RecoupleAngMom(s :: Number, ps_pot :: Vector{<:Number})
    ps_pot0 = Dict([ 2s + 1 - i => ps_pot[i] for i ∈ eachindex(ps_pot) if abs(ps_pot[i]) > 1E-8])
    return RecoupleAngMom(s, s, s, s, ps_pot0)
end
"""
    RecoupleAngMom(s1 :: Number, s2 :: Number, s3 :: Number, s4 :: Number, ch0 :: Vector{Matrix{Int64}}, coeff0 :: Vector{<:Number}) :: Tuple{Vector{Matrix{Int64}}, Vector{<:Number}}

re-couples the angular momentum composition of a four-fermion term ``c^†_1c^†_2c_3c_4`` from the pairing channel ``(12)(34)`` into the density channel ``(14)(23)`` using Wigner's ``9j``-symbol. This is the most general method and allows non-zero total angular momentum. 

# Arguments

* `s1, s2, s3, s4 :: Number` are the spins of the four particles.
* `ch0 :: Vector{Matrix{Int64}}` are the angular-momentum channels.
* `coeff0 :: Vector` are the coefficient of each channel.

# Output

A tuple of re-coupled channels `ch0 :: Vector{Matrix{Int64}}` and co-efficients `coeff0 :: Vector{<:Number}`
"""
function RecoupleAngMom(s1 :: Number, s2 :: Number, s3 :: Number, s4 :: Number, ch0 :: Vector{Matrix{Int64}}, coeff0 :: Vector{<:Number})
    s12 = Int64(2s1)
    s22 = Int64(2s2)
    s32 = Int64(2s3)
    s42 = Int64(2s4)
    ps_pot1 = Dict{NTuple{3, Int64}, ComplexF64}()
    for i in eachindex(ch0)
        l1 = ch0[i][1, 1]
        l2 = ch0[i][1, 2]
        l = ch0[i][2, 2]
        for j1 = abs(s12 - s42) : 2 : s12 + s42, j2 = abs(s22 - s32) : 2 : s22 + s32
            (abs(j1 - j2) ≤ l ≤ j1 + j2 && iseven(j1 + j2 + l)) || continue
            ninej = Float64(float(d9j(s12, s22, l1, s42, s32, l2, j1, j2, l)))
            ninej == 0.0 && continue
            W = √((l1 + 1) * (l2 + 1) * (j1 + 1) * (j2 + 1)) * ninej * coeff0[i]
            W *= (1.0im) ^ Int(-l2 + 2s32 + s42 - s22)
            ps_pot1[(j1, j2, l)] = get(ps_pot1, (j1, j2, l), 0.0im) + W
        end
    end
    ch1 = Matrix{Int64}[]
    coeff1 = ComplexF64[]
    for ((j1, j2, l), W) in ps_pot1
        abs(W) < 1E-13 && continue
        push!(ch1, [j1  j2 ; j1  l])
        push!(coeff1, W)
    end
    return ch1, coeff1
end


"""
    ConvPsPot(ps_pot0 :: Dict) :: Tuple{Vector{Matrix{Int64}}, Vector{ComplexF64}}

converts a pseudo-potential into the coupling channels and coefficients used by a [CoupleDecomps](@ref). The co-efficient is connected to the pseudo-potential by ``Ṽ_l = (-1)^l\\sqrt{2l+1}V_l``.

# Arguments

* `ps_pot0 :: Dict` is a dictionary mapping a density rank ``l`` to its coefficient ``V_l``

# Output

* `ch :: Vector{Matrix{Int64}}` is the list of coupling channels.
* `coeff :: Vector{ComplexF64}` is the list of channel coefficients.
"""
function ConvPsPot(ps_pot0 :: Dict)
    ch = Matrix{Int64}[]
    coeff = ComplexF64[]
    for (J, V) in ps_pot0
        J2 = Int64(2 * J)
        push!(ch, [J2  J2 ; J2  0])
        push!(coeff, V * √(J2 + 1) * (1.0im) ^ (J2))
    end
    return ch, coeff
end
"""
    ConvPsPot(s :: Number, ps_pot :: Vector{<:Number}) :: Tuple{Vector{Matrix{Int64}}, Vector{ComplexF64}}

converts a pseudo-potential into the coupling channels and coefficients used by a [CoupleDecomps](@ref). The co-efficient is connected to the pseudo-potential by ``Ṽ_l = (-1)^l\\sqrt{2l+1}V_l``.

# Arguments

* `s` is the spin of a single electron, and 
* `ps_pot` is a vector of pseudo-potentials from `l=2s` downwards. 

# Output

* `ch :: Vector{Matrix{Int64}}` is the list of coupling channels.
* `coeff :: Vector{ComplexF64}` is the list of channel coefficients.
"""
function ConvPsPot(s :: Number, ps_pot :: Vector{<:Number})
    ps_pot0 = Dict([ 2s + 1 - i => ps_pot[i] for i ∈ eachindex(ps_pot) if abs(ps_pot[i]) > 1E-8])
    return ConvPsPot(ps_pot0)
end


function FuzzifiED.AngModes(obs :: SphereObs)
    return AngModes(obs.l2m, obs.get_comp)
end
function Fuzzifino.SAngModes(obs :: SSphereObs)
    return SAngModes(obs.l2m, obs.get_comp)
end
_ContactAngModes(obs :: SphereObs) = AngModes(obs)
_ContactAngModes(obs :: SSphereObs) = SAngModes(obs)


"""
    ContactCouple(obs :: Vector{<:Union{SphereObs, SSphereObs}}, sec :: Matrix{Int64}, ltot :: Int64) :: CoupleDecomps

constructs a [CoupleDecomps](@ref) for a contact term, _i. e._ the product of one spherical observable per part evaluated at the same point on the sphere
```math 
    ∫\\mathrm{d}^2𝐫\\,\\sqrt{4π}Ȳ_{lm}(𝐫)\\,Φ_1(𝐫)Φ_2(𝐫)⋯Φ_{N_p}(𝐫)
```

# Arguments

* `obs :: Vector` is the list of spherical observables acting on each segment — each a fermionic `SphereObs` or a bosonic `SSphereObs`, and the two may be mixed.
* `sec :: Matrix{Int64}` records the change of quantum numbers, `sec[iqn, p]` for the `iqn`-th quantum number on part ``p``.
* `ltot :: Int64` is twice the total angular momentum ``2l_{\\text{tot}}`` of the term. Facultative, ``0`` (a scalar) by default.

# Output

* `cpd :: CoupleDecomps` is the resulting coupling decompositions, one [CoupleDecomps](@ref) per channel.
"""
function ContactCouple(obs :: Vector, sec :: Matrix{Int64}, ltot :: Int64 = 0)
    amd = _ContactAngModes.(obs)
    np = length(obs)
    s2 = [ obsi.s2 for obsi in obs ]
    s2_ptsum = cumsum(s2)
    l_rng = [ abs(obsi.s2) : 2 : obsi.l2m for obsi in obs ]
    ch1 = Matrix{Int64}[]
    ch = Matrix{Int64}[]
    coeff = ComplexF64[]
    for li in Iterators.product(l_rng...)
        chi = FindCouplingChannels(np, collect(li), ltot)
        append!(ch1, chi)
    end
    for chi in ch1
        coeffi = 1 
        flag = true
        for p = 2 : np 
            if (chi[2, p] < abs(s2_ptsum[p]))
                flag = false 
                break 
            end
            coeffi *= clebschgordan(chi[2, p - 1]/2, -s2_ptsum[p - 1]/2, chi[1, p]/2, -s2[p]/2, chi[2, p]/2, -s2_ptsum[p]/2)
        end
        flag || continue
        coeffi *= √(prod(chi[1, :] .+ 1) / (ltot + 1)) * FuzzifiED.ObsNormRadSq / (4π) ^ (np/2 - 1)
        push!(ch, chi)
        push!(coeff, coeffi)
    end
    return CoupleDecomps(amd, ch, coeff, sec)
end


"""
    SingleSegCouple([np :: Int64, p :: Int64, ]amdp :: Union{AngModes, SAngModes}, l :: Int64, secp :: Vector{Int64}) :: CoupleDecomps

constructs a [CoupleDecomps](@ref) for an angular mode that acts only on a single part ``p`` and as the identity on all the other parts.

# Arguments

* `np :: Int64` is the number of parts. Facultative, 1 by default.
* `p :: Int64` is the index of the part on which the operator acts. Facultative, 1 by default.
* `amdp :: AngModes` or `amdp :: SAngModes` is the angular modes acting on part ``p``.
* `l :: Int64` is twice the rank ``2l`` of the operator on part ``p``.
* `secp :: Vector{Int64}` is the change of quantum numbers on part ``p``.

# Output

* `cpd :: CoupleDecomps` is the resulting coupling decomposition (a single channel).
"""
function SingleSegCouple(np :: Int64, p :: Int64, amdp :: Union{AngModes, SAngModes}, l :: Int64, secp :: Vector{Int64})
    amd = Union{AngModes, SAngModes, Symbol}[:Identity for _ = 1 : np]
    amd[p] = amdp
    ch = zeros(Int64, 2, np)
    ch[1, p] = l
    ch[2, p : end] .= l 
    sec = zeros(Int64, length(secp), np)
    sec[:, p] = secp
    return CoupleDecomps(amd, [ch], [1], sec)
end
"""
    SingleSegCouple([np :: Int64, p :: Int64, ]tms :: Union{Terms, STerms}, sec :: Vector{Int64}) :: CoupleDecomps

constructs a ``\\mathrm{SO}(3)``-spin-``0`` [CoupleDecomps](@ref) for terms that acts only on a single part ``p`` and as the identity on all the other parts.

# Arguments

* `np :: Int64` is the number of parts. Facultative, 1 by default.
* `p :: Int64` is the index of the part on which the operator acts. Facultative, 1 by default.
* `tms :: Terms` or `tms :: Terms` is the ``\\mathrm{SO}(3)``-spin-``0`` terms acting on part ``p``.
* `l :: Int64` is twice the rank ``2l`` of the operator on part ``p``.
* `secp :: Vector{Int64}` is the change of quantum numbers on part ``p``.

# Output

* `cpd :: CoupleDecomps` is the resulting coupling decomposition (a single channel).
"""
function SingleSegCouple(np :: Int64, p :: Int64, tms :: Terms, sec :: Vector{Int64})
    amdp = AngModes(0, Dict((0, 0) => tms))
    return SingleSegCouple(np, p, amdp, 0, sec)
end
function SingleSegCouple(np :: Int64, p :: Int64, tms :: STerms, sec :: Vector{Int64})
    amdp = SAngModes(0, Dict((0, 0) => tms))
    return SingleSegCouple(np, p, amdp, 0, sec)
end
SingleSegCouple(amdp :: AngModes, l :: Int64, secp :: Vector{Int64}) = SingleSegCouple(1, 1, amdp, l, secp)
SingleSegCouple(tms :: Terms, secp :: Vector{Int64}) = SingleSegCouple(1, 1, tms, secp)

"""
    InsertSegment(np :: Int64, p_rng :: Vector{Int64}, ch :: Matrix{Int64}) :: Matrix{Int64}
    InsertSegment(np :: Int64, p_rng :: Vector{Int64}, cpd :: CoupleDecomps) :: CoupleDecomps

Insert segments where the coupling acts as identity to channel(s) or CoupleDecomp(s)
```math
    H_1⊗H_2⊗⋯↦𝕀⊗⋯⊗H_1⊗𝕀⊗⋯⊗H_2⊗𝕀⊗⋯
```

# Arguments 
* `np :: Int64` is the number of segments after the insertion.
* `p_rng :: Int64` is the positions of the non-trivial segments after the insertion.
* `ch :: Matrix{Int64}`, `ch :: Vector{Matrix{Int64}}`, `cpd :: CoupleDecomp`, `cpd :: CoupleDecomps` is the channel(s) or CoupleDecomp(s) before the insertion. 
"""
function InsertSegment(np :: Int64, p_rng :: Vector{Int64}, ch :: Matrix{Int64})
    ch1 = zeros(Int64, 2, np)
    p_rng1 = [p_rng ; np + 1]
    for p = 1 : length(p_rng)
        ch1[:, p_rng[p]] = ch[:, p]
        ch1[2, p_rng[p] + 1 : p_rng1[p + 1] - 1] .= ch[2, p]
    end
    return ch1
end
InsertSegment(np :: Int64, p_rng :: Vector{Int64}, chs :: Vector{Matrix{Int64}}) = InsertSegment.(Ref(np), Ref(p_rng), chs)
function InsertSegment(np :: Int64, p_rng :: Vector{Int64}, cpd :: CoupleDecomp)
    amd1 = Union{AngModes, SAngModes, Symbol}[:Identity for _ = 1 : np]
    amd1[p_rng] = cpd.amd
    ch1 = InsertSegment(np, p_rng, cpd.ch)
    sec1 = zeros(Int64, size(cpd.sec, 1), np)
    sec1[:, p_rng] = cpd.sec
    return CoupleDecomp(amd1, ch1, cpd.coeff, sec1)
end
InsertSegment(np :: Int64, p_rng :: Vector{Int64}, cpds :: CoupleDecomps) = InsertSegment.(Ref(np), Ref(p_rng), cpds)
