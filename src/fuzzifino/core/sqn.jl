export SQNDiag, SQNOffd
export PadSQNDiag, PadSQNOffd


"""
    SQNDiag

The mutable type `SQNDiag` records the information of a diagonal ``\\mathrm{U}(1)`` or ``ℤ_p`` quantum number in the form of a symmetry charge
```math
Q=∑_{o=1}^{N_{of}}q_{f,o}n_{f,o}+∑_{o=1}^{N_{ob}}q_{b,o}n_{b,o}
```
or
```math
Q=∑_{o=1}^{N_{of}}q_{f,o}n_{f,o}+∑_{o=1}^{N_{ob}}q_{b,o}n_{b,o}\\ \\mathrm{mod}\\ p
```
where ``i=1,…,N_U`` is the index of quantum number, ``o`` is the index of site, ``N_{of}`` and ``N_{ob}`` are the number of fermionic and bosonic sites, ``n_{f,o}=f^†_of_o``, ``n_{b,o}=b^†_ob_o``, and ``q_{f,o},q_{b,o}`` are a set of symmetry charges that must be integer valued.

# Fields 

* `name :: String` is the name of the diagonal quantum number 
* `chargef :: Vector{Int64}` is the symmetry charge ``q_{f,o}`` of each site
* `chargeb :: Vector{Int64}` is the symmetry charge ``q_{b,o}`` of each site
* `modul :: Vector{Int64}` is the modulus ``p``, set to 1 for ``\\mathrm{U}(1)`` SQNDiags. 

# Initialisation 

It can be initialised by the following method
```julia
SQNDiag([name :: String, ]chargef :: Vector{Int64}, chargeb :: Vector{Int64}[, modul :: Int64]) :: SQNDiag
```
The arguments `name` and `modul` are facultative. By default `name` is set to `\"QN\"` and `modul` is set to 1. 
"""
mutable struct SQNDiag
    name :: String
    chargef :: Vector{Int64}
    chargeb :: Vector{Int64}
    modul :: Int64
    SQNDiag(name :: String, chargef :: Vector{Int64}, chargeb :: Vector{Int64}, modul :: Int64 = 1) = new(name, chargef, chargeb, modul)
    SQNDiag(chargef :: Vector{Int64}, chargeb :: Vector{Int64}, modul :: Int64 = 1) = new("QN", chargef, chargeb, modul)
end


"""
    SQNDiag(qnd :: QNDiag, nob :: Int64) :: SQNDiag

converts a pure fermionic QNDiag to a SQNDiag with the boson charges set to empty, implemented as 

    SQNDiag(qnd.name, qnd.charge, fill(0, nob), qnd.modul)
"""
SQNDiag(qnd :: QNDiag, nob :: Int64) = SQNDiag(qnd.name, qnd.charge, fill(0, nob), qnd.modul)


"""
    *(fac :: Int64, qnd :: SQNDiag) :: SQNDiag 
    *(qnd :: SQNDiag, fac :: Int64) :: SQNDiag 
    ÷(qnd :: SQNDiag, fac :: Int64) :: SQNDiag 
    -(qnd :: SQNDiag) :: SQNDiag

returns the SQNDiag multiplied or divided by an integer factor, where the charge is multiplied or integer-divided by the factor. For ``ℤ_p`` quantum numbers, their modulus will be multiplied or integer-divided by the absolute value. If `qnd.modul ÷ abs(fac) ≤ 1`, a trivial SQNDiag will be returned.  
"""
function Base.:*(fac :: Int64, qnd :: SQNDiag)
    return SQNDiag(qnd.name, qnd.chargef .* fac, qnd.chargeb .* fac, qnd.modul == 1 ? 1 : (qnd.modul * abs(fac)))
end
function Base.:*(qnd :: SQNDiag, fac :: Int64)
    return fac * qnd
end
function Base.:÷(qnd :: SQNDiag, fac :: Int64)
    if (qnd.modul > 1 && qnd.modul ÷ fac ≤ 1)
        return SQNDiag(qnd.name, qnd.chargef .* 0, qnd.chargeb .* 0, 1)
    else
        return SQNDiag(qnd.name, qnd.chargef .÷ fac, qnd.chargeb .÷ fac, qnd.modul == 1 ? 1 : (qnd.modul ÷ abs(fac)))
    end
end
function Base.:-(qnd :: SQNDiag)
    return (-1) * qnd
end

"""
    +(qnd1 :: SQNDiag, qnd2 :: SQNDiag) :: SQNDiag 
    -(qnd1 :: SQNDiag, qnd2 :: SQNDiag) :: SQNDiag 

returns the sum or substraction of two SQNDiags, whose name is the samea as `qnd1`, charge is the same as `qnd1 ± qnd2`, and modulus is the GCD of `qnd1` and `qnd2`. If `qnd1` and `qnd2` are both ``ℤ_p`` quantum numbers and their modulus are coprime, a trivial SQNDiag will be returned. 
"""
function Base.:+(qnd1 :: SQNDiag, qnd2 :: SQNDiag)
    if (qnd1.modul == 1)
        modul = qnd2.modul 
    elseif (qnd2.modul == 1)
        modul = qnd1.modul
    else
        modul = gcd(qnd1.modul, qnd2.modul)
        if (modul == 1) return SQNDiag(qnd1.name, qnd1.chargef .* 0, qnd1.chargeb .* 0, 1) end
    end
    return SQNDiag(qnd1.name, qnd1.chargef .+ qnd2.chargef, qnd1.chargeb .+ qnd2.chargeb, modul)
end
function Base.:-(qnd1 :: SQNDiag, qnd2 :: SQNDiag)
    return qnd1 + (-1) * qnd2
end


"""
    PadSQNDiag(qnd :: SQNDiag, nofl :: Int64, nobl :: Int64, nofr :: Int64, nobr :: Int64)

adds `nol` empty orbitals to the left and `nor` empty orbitals to the right, implemented as 

    SQNDiag(qnd.name, 
        [ fill(0, nofl) ; qnd.chargef ; fill(0, nofr) ], 
        [ fill(0, nobl) ; qnd.chargeb ; fill(0, nobr) ], 
        qnd.modul)
"""
PadSQNDiag(qnd :: SQNDiag, nofl :: Int64, nobl :: Int64, nofr :: Int64, nobr :: Int64) = SQNDiag(qnd.name, 
    [ fill(0, nofl) ; qnd.chargef ; fill(0, nofr) ], 
    [ fill(0, nobl) ; qnd.chargeb ; fill(0, nobr) ], 
    qnd.modul)


"""
    SQNOffd

The mutable type `SQNOffd` records the information of an off-diagonal ``ℤ_p`` quantum number in the form of a discrete transformation
```math
𝒵:\\ f_o↦ α_{f,o}^* f^{(p_{f,o})}_{π_{f,o}},  f_o^†↦α_{f,o} c^{(1-p_{f,o})}_{π_{f,o}},  b_o^†↦α_{b,o} b^†_{π_{b,o}}
```
where we use a notation ``c^{(1)}=c^†`` and ``c^{0}=c`` for convenience, ``π_{f,o},π_{b,o}`` are permutations of ``1,…,N_{of}`` or ``N_{ob}``, ``α_{f,o},α_{b,o}`` are coefficients, and ``p_{f,o}`` specified whether or not particle-hole transformation is performed for the fermionic site. Note that one must guarentee that all these transformations commute with each other and also commute with the diagonal QNs. 

# Arguments 

* `permf :: Vector{Int64}` is a length-``N_{of}`` vector that records the fermion permutation ``π_{f,o}``.
* `permb :: Vector{Int64}` is a length-``N_{ob}`` vector that records the boson permutation ``π_{b,o}``.
* `phf :: Vector{Int64}` is a length-``N_{of}`` vector that records ``p_{f,o}`` to determine whether or not to perform a particle-hole transformation
* `facf :: Vector{ComplexF64}` is a length-``N_{of}`` vector that records the factor ``α_{f,o}`` in the transformation.
* `facb :: Vector{ComplexF64}` is a length-``N_{ob}`` vector that records the factor ``α_{b,o}`` in the transformation.
* `cyc :: Int64` is the cycle ``p``. 

# Initialisation 

It can be initialised by the following method
```julia
SQNOffd(permf :: Vector{Int64}, permb :: Vector{Int64}[, phf :: Vector{Int64}][, facf :: Vector{ComplexF64}, facb :: Vector{ComplexF64}][, cyc :: Int64]) :: SQNOffd
SQNOffd(permf :: Vector{Int64}, permb :: Vector{Int64}, phf_q :: Bool[, fac :: Vector{ComplexF64}, facb :: Vector{ComplexF64}]) :: SQNOffd
```
The arguments `phf`, `facf`, `facb` and `cyc` are facultative. By default `ph` is set all 0, `facf`, `facb` is set to all 1 and `cyc` is set to 2. If `phf_q` is a bool and true, then `ph` is set to all 1. 

"""
mutable struct SQNOffd
    permf :: Vector{Int64}
    permb :: Vector{Int64}
    phf :: Vector{Int64}
    facf :: Vector{ComplexF64}
    facb :: Vector{ComplexF64}
    cyc :: Int64
    SQNOffd(permf :: Vector{Int64}, permb :: Vector{Int64}, phf :: Vector{Int64} = [0 for o in permf], facf :: Vector{ComplexF64} = [ComplexF64(1) for o in permf], facb :: Vector{ComplexF64} = [ComplexF64(1) for o in permb], cyc :: Int64 = 2) = new(permf, permb, phf, facf, facb, cyc)
    SQNOffd(permf :: Vector{Int64}, permb :: Vector{Int64}, facf :: Vector{ComplexF64}, facb :: Vector{ComplexF64}, cyc :: Int64 = 2) = new(permf, permb, [0 for o in permf], facf, facb, cyc)
    SQNOffd(permf :: Vector{Int64}, permb :: Vector{Int64}, cyc :: Int64) = new(permf, permb, [0 for o in permf], [ComplexF64(1) for o in permf], [ComplexF64(1) for o in permb], cyc)
    SQNOffd(permf :: Vector{Int64}, permb :: Vector{Int64}, phf_q :: Bool, facf :: Vector{ComplexF64} = [ComplexF64(1) for o in permf], facb :: Vector{ComplexF64} = [ComplexF64(1) for o in permb]) = new(permf, permb, [phf_q ? 1 : 0 for o in permf], facf, facb, 2)
end


"""
    SQNOffd(qnf :: QNOffd, nob :: Int64) :: SQNOffd

converts a pure fermionic QNOffd to a SQNOffd with the boson transformations set to identity, implemented as 

    SQNOffd(qnf.perm, collect(1 : nob), qnf.ph, qnf.fac, fill(ComplexF64(1), nob), qnf.cyc)
"""
SQNOffd(qnf :: QNOffd, nob :: Int64) = SQNOffd(qnf.perm, collect(1 : nob), qnf.ph, qnf.fac, fill(ComplexF64(1), nob), qnf.cyc)


"""
    *(qnf1 :: SQNOffd, qnf2 :: SQNOffd) :: SQNOffd 

returns the composition of two SQNOffd transformations. The cycle is set to be the LCM of two QNOffds.
"""
function Base.:*(qnf1 :: SQNOffd, qnf2 :: SQNOffd)
    permf1 = [ qnf1.permf[qnf2.permf[o]] for o in eachindex(qnf1.permf)]
    permb1 = [ qnf1.permb[qnf2.permb[o]] for o in eachindex(qnf1.permb)]
    phf1 = [ qnf1.phf[qnf2.permf[o]] ⊻ qnf2.phf[o] for o in eachindex(qnf1.permf) ]
    facf1 = [ (qnf2.phf[o] == 0 ? qnf1.facf[qnf2.permf[o]] : conj(qnf1.facf[qnf2.permf[o]])) * qnf2.facf[o] for o in eachindex(qnf1.permf) ]
    facb1 = [ qnf1.facb[qnf2.permb[o]] * qnf2.facb[o] for o in eachindex(qnf1.permb) ]
    cyc1 = lcm(qnf1.cyc, qnf2.cyc)
    return SQNOffd(permf1, permb1, phf1, facf1, facb1, cyc1)
end


"""
    PadSQNOffd(qnf :: SQNOffd, nofl :: Int64, nobl :: Int64, nofr :: Int64, nobr :: Int64)

adds `nofl` fermionic and `nobl` bosonic empty orbitals to the left and `nofr` fermionic and `nobr` bosonic empty orbitals to the right, implemented as 

    SQNOffd(
        [ collect(1 : nofl) ; qnf.permf .+ nofl ; collect(1 : nofr) .+ (nofl + length(qnf.permf)) ], 
        [ collect(1 : nobl) ; qnf.permb .+ nobl ; collect(1 : nobr) .+ (nobl + length(qnf.permb)) ], 
        [ fill(0, nofl) ; qnf.phf ; fill(0, nofr) ], 
        [ fill(ComplexF64(1), nofl) ; qnf.facf ; fill(ComplexF64(1), nofr) ], 
        [ fill(ComplexF64(1), nobl) ; qnf.facb ; fill(ComplexF64(1), nobr) ], 
        qnf.cyc)
"""
PadSQNOffd(qnf :: SQNOffd, nofl :: Int64, nobl :: Int64, nofr :: Int64, nobr :: Int64) = SQNOffd(
    [ collect(1 : nofl) ; qnf.permf .+ nofl ; collect(1 : nofr) .+ (nofl + length(qnf.permf)) ], 
    [ collect(1 : nobl) ; qnf.permb .+ nobl ; collect(1 : nobr) .+ (nobl + length(qnf.permb)) ], 
    [ fill(0, nofl) ; qnf.phf ; fill(0, nofr) ], 
    [ fill(ComplexF64(1), nofl) ; qnf.facf ; fill(ComplexF64(1), nofr) ], 
    [ fill(ComplexF64(1), nobl) ; qnf.facb ; fill(ComplexF64(1), nobr) ], 
    qnf.cyc)
