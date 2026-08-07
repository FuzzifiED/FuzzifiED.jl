export JackToolkit 

module JackToolkit 

using FuzzifiED
using FuzzifiED.Fuzzifino
using LinearAlgebra

export GetJackRoots, GetSqueezedParts, GetJackCoefficients, GetJackState, GetJackStates, OrganizeJackStates

ConfToPart(nocc :: Vector{Int64}) = reduce(vcat, [ fill(m - 1, nocc[m]) for m = length(nocc) : -1 : 1 ] ; init = Int64[])
function PartToConf(part :: Vector{Int64}, nm :: Int64)
    nocc = zeros(Int64, nm)
    for j in part ; nocc[j + 1] += 1 end
    return nocc
end
function SortPart!(part :: Vector{Int64}, fermion :: Bool)
    sgn = 1
    for i = 2 : length(part)
        v = part[i]
        j = i - 1
        while j ≥ 1 && part[j] < v
            part[j + 1] = part[j]
            j -= 1
            sgn = -sgn
        end
        part[j + 1] = v
    end
    fermion || return 1
    for i = 1 : length(part) - 1
        (part[i] == part[i + 1]) && return 0
    end
    return sgn
end
function LogBinomials(n :: Int64)
    lnfct = [0.0 ; cumsum(log.(1.0 : n))]
    return [ lnfct[n + 1] - lnfct[j + 1] - lnfct[n - j + 1] for j = 0 : n ]
end


"""
    GetJackRoots(nm :: Int64, ne :: Int64, k :: Int64, r :: Int64, lz2 :: Int64 = 0 ; fermion :: Bool = true) :: Vector{Vector{Int64}}

enumerates the ``(k,r)``-admissible roots with `nm` orbitals ``N_m``, `ne` particles ``N_e`` and total ``2L^z`` `lz2`, _i.e._, the occupation vectors obeying the clustering condition
```math
    n_j+n_{j+1}+⋯+n_{j+r-1}≤k
```
together with ``n_j≤1`` for fermions and ``n_j≤k`` for bosons. Each root generates one Jack state. Returns the occupation vectors `nocc[m], m = 1 : nm` in decreasing lexicographic order.
"""
function GetJackRoots(nm :: Int64, ne :: Int64, k :: Int64, r :: Int64, lz2 :: Int64 = 0 ; fermion :: Bool = true) :: Vector{Vector{Int64}}
    roots = Vector{Int64}[]
    (ne == 0) && return iszero(lz2) ? [zeros(Int64, nm)] : roots
    (nm ≤ 0 || ne < 0) && return roots
    isodd(lz2 - ne * (nm - 1)) && return roots
    nocp_cap = fermion ? min(k, 1) : k

    # lz2_hi[m, n] / lz2_lo[m, n] : maximal / minimal 2Lz carried by n particles
    # put in the orbitals m : nm, ignoring the (k,r) rule ; -1 if impossible.
    lz2_hi = fill(0, nm + 1, ne + 1)
    lz2_lo = fill(0, nm + 1, ne + 1)
    psb = fill(false, nm + 1, ne + 1)
    for m = nm + 1 : -1 : 1, n = 0 : ne
        if (n == 0)
            psb[m, n + 1] = true
            continue
        end
        (m > nm) && continue
        nocp = min(nocp_cap, n)
        psb[m, n + 1] = psb[m + 1, n - nocp + 1]
        psb[m, n + 1] || continue
        # fill from the top for the maximum, from the bottom for the minimum
        n_hi = n ; lz2_hi[m, n + 1] = 0
        for m1 = nm : -1 : m
            nocp1 = min(nocp_cap, n_hi)
            lz2_hi[m, n + 1] += nocp1 * (2 * m1 - nm - 1)
            n_hi -= nocp1
            (n_hi == 0) && break
        end
        n_lo = n ; lz2_lo[m, n + 1] = 0
        for m1 = m : nm
            nocp1 = min(nocp_cap, n_lo)
            lz2_lo[m, n + 1] += nocp1 * (2 * m1 - nm - 1)
            n_lo -= nocp1
            (n_lo == 0) && break
        end
    end

    nocc = zeros(Int64, nm)
    function AddOrb!(m :: Int64, n_rest :: Int64, lz2_rest :: Int64)
        if (m > nm)
            (n_rest == 0 && lz2_rest == 0) && push!(roots, copy(nocc))
            return
        end
        psb[m, n_rest + 1] || return
        (lz2_rest > lz2_hi[m, n_rest + 1] || lz2_rest < lz2_lo[m, n_rest + 1]) && return
        # the (k,r) exclusion rule : the window of r orbitals ending at m
        n_win = sum(nocc[max(1, m - r + 1) : m - 1] ; init = 0)
        for nocp = min(nocp_cap, k - n_win, n_rest) : -1 : 0
            nocc[m] = nocp
            AddOrb!(m + 1, n_rest - nocp, lz2_rest - nocp * (2 * m - nm - 1))
        end
        nocc[m] = 0
    end
    AddOrb!(1, ne, lz2)
    return roots
end


"""
    GetSqueezedParts(part_rt :: Vector{Int64} ; fermion :: Bool = true) :: Vector{Vector{Int64}}

returns the partitions squeezed from the root `part_rt`, _i.e._, dominated by it, ``∑_{i≤t}μ_i≤∑_{i≤t}λ_i`` for every `t`, in decreasing lexicographic order — a refinement of the dominance order, so that `part_rt` comes first and every partition follows those that dominate it. For fermions the partitions are strictly decreasing, obtained from the bosonic ones through the dominance-preserving bijection ``μ↦μ+δ``, ``δ=(N_e-1,…,1,0)``.
"""
function GetSqueezedParts(part_rt :: Vector{Int64} ; fermion :: Bool = true)
    ne = length(part_rt)
    (ne == 0) && return [Int64[]]
    δ = fermion ? collect(ne - 1 : -1 : 0) : zeros(Int64, ne)
    λ = part_rt .- δ
    any(<(0), λ) && error("the root is not a valid partition")
    lλ = cumsum(λ)
    parts = Vector{Int64}[]
    μ = zeros(Int64, ne)
    function AddPart!(i :: Int64, prev :: Int64, wt :: Int64)
        if (i > ne)
            push!(parts, μ .+ δ)
            return
        end
        wt_rest = lλ[ne] - wt
        # the ne - i + 1 remaining parts are all ≤ μ[i] and sum up to wt_rest
        μ_lo = cld(wt_rest, ne - i + 1)
        μ_hi = min(prev, lλ[i] - wt, wt_rest)
        for v = μ_hi : -1 : μ_lo
            μ[i] = v
            AddPart!(i + 1, v, wt + v)
        end
        μ[i] = 0
    end
    AddPart!(1, λ[1], 0)
    return parts
end


"""
    GetJackCoefficients(part_rt :: Vector{Int64}, k :: Int64, r :: Int64 ; fermion :: Bool = true) :: Tuple{Vector{Vector{Int64}}, Vector{Float64}}

returns the expansion of the ``(k,r)`` model state — Laughlin ``(1,q)``, Moore-Read ``(2,2q)``, Read-Rezayi ``(k,r)`` — with root partition `part_rt` in the monomials (bosons) or the Slater determinants (fermions) of the squeezed partitions.

For bosons the state is the Jack polynomial at ``α=-(k+1)/(r-1)``, ``J^α_λ=∑_{μ≤λ}c_{λμ}m_μ`` with ``c_{λλ}=1``, whose coefficients follow from the recursion of the Laplace–Beltrami operator
```math
    c_{λμ}=\\frac{1}{E_λ-E_μ}∑_{ν}(ν_i-ν_j)c_{λν},\\qquad E_μ=∑_i\\left[\\frac{α}{2}μ_i(μ_i-1)-(i-1)μ_i\\right]
```
``ν`` running over the partitions obtained by un-squeezing a pair of ``μ``, ``(μ_i,μ_j)↦(μ_i+s,μ_j-s)`` with ``i<j`` and ``s≥1``, reordered decreasingly. For ``r=1``, ``α=∞`` and the Jack degenerates to the single monomial ``m_λ``.

For fermions the state is _not_ the Jack of the fermionic partition, but the bosonic Jack of the ``(k,r-k)`` root ``λ-δ`` times the Vandermonde determinant ``Δ=∏_{i<j}(z_i-z_j)``,
```math
    Ψ^{(k,r)}_λ=Δ(z)\\,J^{α}_{λ-δ}(z),\\qquad δ=(N_e-1,…,1,0),\\qquad α=-\\frac{k+1}{r-k-1}
```
which is the construction that gives the model states as the zero modes of the clustering Hamiltonians. The product is taken monomial by monomial into the Slater determinants ``\\mathrm{sl}_ν=\\det z_i^{ν_j}``,
```math
    Δ\\,m_μ=\\frac{1}{|\\mathrm{Aut}\\,μ|}∑_{w∈S_{N_e}}\\mathrm{sgn}(w)\\,\\det z_i^{μ_j+w_j}
```
summing over the assignments of the ``N_e`` distinct staircase powers ``w_j∈\\{0,…,N_e-1\\}`` to the parts of ``μ``, dropping those that give two equal powers ``μ_j+w_j`` and keeping one assignment per block of equal parts, which cancels ``1/|\\mathrm{Aut}\\,μ|``. _N. b._, this last step costs up to ``N_e!`` per squeezed partition and is the bottleneck for large ``N_e``.

Returns the squeezed partitions `parts`, with `parts[1] == part_rt`, and the coefficients `coeff`, with `coeff[1] == 1`.
"""
function GetJackCoefficients(part_rt :: Vector{Int64}, k :: Int64, r :: Int64 ; fermion :: Bool = true)
    ne = length(part_rt)
    (ne == 0) && return [Int64[]], [1.0]
    if (fermion)
        (r > k) || error("for fermions r must be greater than k")
    else
        (r > 0) || error("r must be positive")
    end
    # the bosonic problem : the root λ - δ at (k, r - k) for fermions, λ at (k, r) for bosons
    δ = fermion ? collect(ne - 1 : -1 : 0) : zeros(Int64, ne)
    rb = fermion ? r - k : r
    parts_b = GetSqueezedParts(part_rt .- δ ; fermion = false)
    coeff_b = zeros(Float64, length(parts_b))
    coeff_b[1] = 1.0

    # the recursion of the Laplace-Beltrami operator ; for rb = 1, α = ∞ and J = m_λ
    if (rb > 1)
        α = -(k + 1) / (rb - 1)
        id_part = Dict{Vector{Int64}, Int64}(parts_b[i] => i for i in eachindex(parts_b))
        EigenVal(μ) = sum([ .5 * α * μ[i] * (μ[i] - 1) - (i - 1) * μ[i] for i = 1 : ne ] ; init = 0.0)
        eg_rt = EigenVal(parts_b[1])
        ν = zeros(Int64, ne)
        for iμ = 2 : length(parts_b)
            μ = parts_b[iμ]
            num = 0.0
            for i = 1 : ne - 1, j = i + 1 : ne
                for s = 1 : μ[j]
                    vi = μ[i] + s
                    vj = μ[j] - s
                    (vi > parts_b[1][1]) && break # ν₁ ≥ vi cannot be dominated by λ
                    copyto!(ν, μ)
                    ν[i] = vi
                    ν[j] = vj
                    SortPart!(ν, false)
                    iν = get(id_part, ν, 0)
                    (iν == 0) && continue # ν is not dominated by λ
                    num += (vi - vj) * coeff_b[iν]
                end
            end
            den = eg_rt - EigenVal(μ)
            (abs(den) < 1E-10) && error("singular Jack : E_λ = E_μ for λ = $(parts_b[1]), μ = $μ ; the root is probably not (k,r)-admissible, or k+1 and r-1 are not coprime")
            coeff_b[iμ] = num / den
        end
    end
    fermion || return parts_b, coeff_b

    # the multiplication by the Vandermonde determinant, monomial by monomial
    coeff_f = Dict{Vector{Int64}, Float64}()
    used_w = falses(ne)
    used_pw = falses(parts_b[1][1] + ne)
    pw = zeros(Int64, ne)
    function AddStaircase!(μ :: Vector{Int64}, cf_μ :: Float64, p :: Int64, w_prev :: Int64, sgn_w :: Int64)
        if (p > ne)
            pw1 = copy(pw)
            sgn_pw = SortPart!(pw1, true)
            coeff_f[pw1] = get(coeff_f, pw1, 0.0) + sgn_w * sgn_pw * cf_μ
            return
        end
        # within a block of equal parts of μ the staircase powers are increasing
        for w = (p > 1 && μ[p] == μ[p - 1] ? w_prev + 1 : 0) : ne - 1
            used_w[w + 1] && continue
            pw[p] = μ[p] + w
            used_pw[pw[p] + 1] && continue # two equal powers : the Slater vanishes
            # each already assigned w₁ < w counts as one inversion of the sequence
            sgn_w1 = sgn_w * (-1) ^ count(used_w[1 : w])
            used_w[w + 1] = true ; used_pw[pw[p] + 1] = true
            AddStaircase!(μ, cf_μ, p + 1, w, sgn_w1)
            used_w[w + 1] = false ; used_pw[pw[p] + 1] = false
        end
    end
    for iμ in eachindex(parts_b)
        (abs(coeff_b[iμ]) < 1E-13) && continue
        AddStaircase!(parts_b[iμ], coeff_b[iμ], 1, -1, 1)
    end
    # μ ↦ μ + δ is the bijection of [GetSqueezedParts](@ref) onto the fermionic partitions
    parts = [ μ .+ δ for μ in parts_b ]
    return parts, [ get(coeff_f, p, 0.0) for p in parts ]
end



"""
    GetJackState(bs :: Basis, nm :: Int64, ne :: Int64, k :: Int64, r :: Int64, root :: Vector{Int64}) :: Vector{Float64}
    GetJackState(bs :: SBasis, nm :: Int64, ne :: Int64, k :: Int64, r :: Int64, root :: Vector{Int64}) :: Vector{Float64}

returns the normalised ``(k,r)`` model state generated by the root configuration `root` — an occupation vector of length `nm`, _e. g._, an element of [GetJackRoots](@ref) — as a vector of length `bs.dim`, _cf._ [GetJackCoefficients](@ref). The particles are fermions for a `Basis` and bosons for an `SBasis`. If `bs` carries QNOffds, the projection of the state onto the sector is returned.
"""
function GetJackState(bs :: Basis, nm :: Int64, ne :: Int64, k :: Int64, r :: Int64, root :: Vector{Int64} )
    (length(root) == nm) || error("the root must be an occupation vector of length nm")
    (sum(root) == ne) || error("the root must contain ne particles")
    parts, coeff = GetJackCoefficients(ConfToPart(root), k, r ; fermion = true)
    lnbn = LogBinomials(nm - 1)
    st = zeros(Float64, bs.dim)
    for iμ in eachindex(parts)
        (abs(coeff[iμ]) < 1E-13) && continue
        cf = 0
        lnnm = 0.0
        for j in parts[iμ]
            cf |= 1 << j
            lnnm += .5 * lnbn[j + 1]
        end
        amp = coeff[iμ] * exp(-lnnm)
        id = GetConfId(bs.cfs, cf)
        (id ≤ 0 || id > bs.cfs.ncf) && error("the configuration $(parts[iμ]) is not in the basis ; check ne and lz2")
        igr = bs.cfgr[id]
        (igr ≤ 0) && continue
        st[igr] += real(conj(bs.cffac[id])) * amp
    end
    return normalize(st)
end
function GetJackState(bs :: SBasis, nm :: Int64, ne :: Int64, k :: Int64, r :: Int64, root :: Vector{Int64})
    (length(root) == nm) || error("the root must be an occupation vector of length nm")
    (sum(root) == ne) || error("the root must contain ne particles")
    parts, coeff = GetJackCoefficients(ConfToPart(root), k, r ; fermion = false)
    lnbn = LogBinomials(nm - 1)
    lnfct = [0.0 ; cumsum(log.(1.0 : ne))]
    st = zeros(Float64, bs.dim)
    for iμ in eachindex(parts)
        (abs(coeff[iμ]) < 1E-13) && continue
        nocc = PartToConf(parts[iμ], nm)
        lnnm = 0.0
        for m = 1 : nm
            (nocc[m] == 0) && continue
            # ∏_j C_j^{n_j / 2} √(n_j !) relating m_μ to the normalised Fock state
            lnnm += .5 * (nocc[m] * lnbn[m] + lnfct[nocc[m] + 1])
        end
        amp = coeff[iμ] * exp(-lnnm)
        cfb = EncodeConfb(nocc, bs.cfs.nebm)
        id = GetSConfId(bs.cfs, 0, nocc)
        (id ≤ 0 || id > bs.cfs.ncf || bs.cfs.conff[id] ≠ 0 || bs.cfs.confb[id] ≠ cfb) &&
            error("the configuration $(parts[iμ]) is not in the basis ; check ne, lz2 and nebm")
        igr = bs.cfgr[id]
        (igr ≤ 0) && continue
        st[igr] += real(conj(bs.cffac[id])) * amp
    end
    return normalize(st)
end


"""
    GetJackStates(bs :: Basis, nm :: Int64, ne :: Int64, k :: Int64, r :: Int64, lz2 :: Int64 = 0) :: Matrix{Float64}
    GetJackStates(bs :: SBasis, nm :: Int64, ne :: Int64, k :: Int64, r :: Int64, lz2 :: Int64 = 0) :: Matrix{Float64}

returns the normalised Jack states generated by all the ``(k,r)``-admissible roots in the sector `(ne, lz2)`, _cf._ [GetJackRoots](@ref) and [GetJackState](@ref), as the columns of a `bs.dim`×`nroot` matrix. _N. b._, they are neither orthogonal nor eigenstates of the total angular momentum.
"""
function GetJackStates(bs :: Union{SBasis, Basis}, nm :: Int64, ne :: Int64, k :: Int64, r :: Int64, lz2 :: Int64 = 0)
    roots = GetJackRoots(nm, ne, k, r, lz2 ; fermion = (typeof(bs) == Basis))
    sts = zeros(Float64, bs.dim, length(roots))
    Threads.@threads for i in eachindex(roots)
        sts[:, i] = GetJackState(bs, nm, ne, k, r, roots[i])
    end
    return sts
end


"""
    OrganizeJackStates(sts :: Matrix{T}, l2_mat :: OpMat{T}) :: Tuple{Vector{T}, Matrix{T}}

takes in the Jack states generated from `GetJackStates`, orthonomalizes them and diagonalize with respect to the total angular momentum specified by `l2_mat`, and returns the vector of total angular momentum and the eigen-states. 
"""
function OrganizeJackStates(sts :: Matrix{T}, l2_mat :: OpMat{T}) where T <: Union{Float64, ComplexF64}
    sts1 = Matrix(qr(sts).Q) 
    l2_st = sts1' * l2_mat * sts1
    l2_val, Λ = eigen(Hermitian(l2_st))
    sts2 = sts1 * Λ 
    return l2_val, sts2
end
"""
    OrganizeJackStates(sts :: Matrix{T}, bs :: Basis, nm :: Int64) :: Tuple{Vector{T}, Matrix{T}}
    OrganizeJackStates(sts :: Matrix{T}, bs :: SBasis, nm :: Int64) :: Tuple{Vector{T}, Matrix{T}}

takes in the Jack states generated from `GetJackStates`, orthonomalizes them and diagonalize with respect to the total angular momentum, and returns the vector of total angular momentum and the eigen-states. The total angular momentum is generated from the basis `bs :: Basis` or `bs :: Basis` and the number of orbitals `nm`.
"""
function OrganizeJackStates(sts :: Matrix{T}, bs :: Basis, nm :: Int64) where T <: Union{Float64, ComplexF64}
    tms_l2 = GetL2Terms(nm, 1)
    l2_mat = OpMat(Operator(bs, tms_l2))
    return OrganizeJackStates(sts, l2_mat)
end 
function OrganizeJackStates(sts :: Matrix{T}, bs :: SBasis, nm :: Int64) where T <: Union{Float64, ComplexF64}
    tms_l2 = GetBosonL2STerms(nm, 1)
    l2_mat = OpMat(SOperator(bs, tms_l2))
    return OrganizeJackStates(sts, l2_mat)
end 


end