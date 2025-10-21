export QNDiag, QNOffd


"""
    QNDiag

The mutable type `QNDiag` records the information of a diagonal ``\\mathrm{U}(1)`` or ``ℤ_p`` quantum number in the form of a symmetry charge
```math
Q=∑_{o=1}^{N_o}q_on_o
```
or
```math
Q=∑_{o=1}^{N_o}q_on_o\\ \\mathrm{mod}\\ p
```
where ``i=1,…,N_U`` is the index of quantum number, ``o`` is the index of site, ``n_o=c^†_oc_o``, and ``q_o`` is a set of coefficients that must be integer valued.

# Fields 

* `name :: String` is the name of the diagonal quantum number.
* `charge :: Vector{Int64}` is the symmetry charge ``q_o`` of each site.
* `modul ::Int64` is the modulus ``p``, set to 1 for ``\\mathrm{U}(1)`` QNDiags. 

# Initialisation 

It can be initialised by the following method
```julia
QNDiag([name :: String, ]charge :: Vector{Int64}[, modul :: Int64])
```
The arguments `name` and `modul` are facultative. By default `name` is set to `\"QN\"` and `modul` is set to 1. 
"""
mutable struct QNDiag
    name :: String
    charge :: Vector{Int64}
    modul :: Int64
    QNDiag(name :: String, charge :: Vector{Int64}, modul :: Int64 = 1) = new(name, charge, modul)
    QNDiag(charge :: Vector{Int64}, modul :: Int64 = 1) = new("QN", charge, modul)
end


"""
    *(fac :: Int64, qnd :: QNDiag) :: QNDiag 
    *(qnd :: QNDiag, fac :: Int64) :: QNDiag 
    ÷(qnd :: QNDiag, fac :: Int64) :: QNDiag 
    -(qnd :: QNDiag) :: QNDiag 

returns the QNDiag multiplied or divided by an integer factor, where the charge is multiplied or integer-divided by the factor. For ``ℤ_p`` quantum numbers, their modulus will be multiplied or integer-divided by the absolute value. If `qnd.modul ÷ abs(fac) ≤ 1`, a trivial QNDiag will be returned.  
"""
function Base.:*(fac :: Int64, qnd :: QNDiag)
    return QNDiag(qnd.name, qnd.charge .* fac, qnd.modul == 1 ? 1 : (qnd.modul * abs(fac)))
end
function Base.:*(qnd :: QNDiag, fac :: Int64)
    return fac * qnd
end
function Base.:÷(qnd :: QNDiag, fac :: Int64)
    if (qnd.modul > 1 && qnd.modul ÷ fac ≤ 1)
        return QNDiag(qnd.name, qnd.charge .* 0, 1)
    else
        return QNDiag(qnd.name, qnd.charge .÷ fac, qnd.modul == 1 ? 1 : (qnd.modul ÷ abs(fac)))
    end
end
function Base.:-(qnd :: QNDiag)
    return (-1) * qnd
end


"""
    +(qnd1 :: QNDiag, qnd2 :: QNDiag) :: QNDiag 
    -(qnd1 :: QNDiag, qnd2 :: QNDiag) :: QNDiag 

returns the sum or substraction of two QNDiags, whose name is the samea as `qnd1`, charge is the same as `qnd1 ± qnd2`, and modulus is the GCD of `qnd1` and `qnd2`. If `qnd1` and `qnd2` are both ``ℤ_p`` quantum numbers and their modulus are coprime, a trivial QNDiag will be returned. 
"""
function Base.:+(qnd1 :: QNDiag, qnd2 :: QNDiag)
    if (qnd1.modul == 1)
        modul = qnd2.modul 
    elseif (qnd2.modul == 1)
        modul = qnd1.modul
    else
        modul = gcd(qnd1.modul, qnd2.modul)
        if (modul == 1) return QNDiag(qnd1.name, qnd1.charge .* 0, 1) end
    end
    return QNDiag(qnd1.name, qnd1.charge .+ qnd2.charge, modul)
end
function Base.:-(qnd1 :: QNDiag, qnd2 :: QNDiag)
    return qnd1 + (-1) * qnd2
end


"""
    QNOffd 

The mutable type `QNOffd` records the information of an off-diagonal ``ℤ_p`` quantum number in the form of a discrete transformation
```math
𝒵:\\ c_o↦ α_o^* c^{(p_o)}_{π_o},  c_o^†↦ α_o c^{(1-p_o)}_{π_o}
```
where we use a notation ``c^{(1)}=c^†`` and ``c^{0}=c`` for convenience, ``π_o`` is a permutation of ``1,…,N_o``, ``α_o`` is a coefficient, and ``p_o`` specified whether or not particle-hole transformation is performed for the site. Note that one must guarentee that all these transformations commute with each other and also commute with the diagonal QNs. 

# Arguments 

* `perm :: Vector{Int64}` is a length-``N_o`` vector that records the permutation ``π_o``.
* `ph :: Vector{Int64}` is a length-``N_o`` vector that records ``p_o`` to determine whether or not to perform a particle-hole transformation.
* `fac :: Vector{ComplexF64}` is a length-``N_o`` vector that records the factor ``α_o`` in the transformation.
* `cyc :: Int64` is the cycle ``p``. 

# Initialisation 

It can be initialised by the following method
```julia
QNOffd(perm :: Vector{Int64}[, ph :: Vector{Int64}][, fac :: Vector{ComplexF64}][, cyc :: Int64])
QNOffd(perm :: Vector{Int64}, ph_q :: Bool[, fac :: Vector{ComplexF64}])
```
The arguments `ph`, `fac` and `cyc` are facultative. By default `ph` is set all 0, `fac` is set to all 1 and `cyc` is set to 2. If `ph_q` is a bool and true, then `ph` is set to all 1. 

"""
mutable struct QNOffd
    perm :: Vector{Int64}
    ph :: Vector{Int64}
    fac :: Vector{ComplexF64}
    cyc :: Int64
    QNOffd(perm :: Vector{Int64}, ph :: Vector{Int64} = [0 for o in perm], fac :: Vector{ComplexF64} = [ComplexF64(1) for o in perm], cyc :: Int64 = 2) = new(perm, ph, fac, cyc)
    QNOffd(perm :: Vector{Int64}, fac :: Vector{ComplexF64}, cyc :: Int64 = 2) = new(perm, [0 for o in perm], fac, cyc)
    QNOffd(perm :: Vector{Int64}, cyc :: Int64) = new(perm, [0 for o in perm], [ComplexF64(1) for o in perm], cyc)
    QNOffd(perm :: Vector{Int64}, ph_q :: Bool, fac :: Vector{ComplexF64} = [ComplexF64(1) for o in perm]) = new(perm, [ph_q ? 1 : 0  for o in perm], fac, 2)
end


"""
    *(qnf1 :: QNOffd, qnf2 :: QNOffd) :: QNOffd 

returns the composition of two QNOffd transformations. The cycle is set to be the LCM of two QNOffds.
"""
function Base.:*(qnf1 :: QNOffd, qnf2 :: QNOffd)
    perm1 = [ qnf1.perm[qnf2.perm[o]] for o in eachindex(qnf1.perm)]
    ph1 = [ qnf1.ph[qnf2.perm[o]] ⊻ qnf2.ph[o] for o in eachindex(qnf1.perm) ]
    fac1 = [ (qnf2.ph[o] == 0 ? qnf1.fac[qnf2.perm[o]] : conj(qnf1.fac[qnf2.perm[o]])) * qnf2.fac[o] for o in eachindex(qnf1.perm) ]
    cyc1 = lcm(qnf1.cyc, qnf2.cyc)
    return QNOffd(perm1, ph1, fac1, cyc1)
end