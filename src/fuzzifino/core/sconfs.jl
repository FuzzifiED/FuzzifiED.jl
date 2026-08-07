export SConfs
export EncodeConfb, DecodeConfb, GetSConfId


"""
    SConfs
    
The mutable type `SConfs` stores all the configurations that respects the diagonal quantum numbers (SQNDiag) and also a table to inversely look up the index from the configuration. 

# Fields

* `nof :: Int64` is the number of fermionic sites.
* `nof :: Int64` is the number of bosonic sites.
* `nebm :: Int64` is the maximal number of boson occupation.
* `ncf :: Int64` is the number of configurations.
* `conff :: Vector{Int64}` is an array of length `ncf` containing all the fermion configurations. Each configuration is expressed in a binary number. If the `o-1`-th bit of `conf[i]` is 1, then the `o`-th site in the `i`-th configuration is occupied ; if the bit is 0, then the site is empty. 
* `confb :: Vector{Int64}` is an array of length `ncf` containing all the boson configurations. Each configuration is expressed in a binary number that has ``N_{b,o}`` 1's and ``N_{eb}`` 0's and the number of 0's following each 1 records the number of bosons in that site. 
* `norf :: Int64`, `nobf :: Int64`, `lid :: Vector{Int64}` and `rid :: Vector{Int64}` contain the information of Lin table that is used to inversely look up the index `i` from the configuration. 
"""
mutable struct SConfs
    nof :: Int64 
    nob :: Int64
    norf :: Int64 
    norb :: Int64
    nebm :: Int64
    ncf :: Int64
    conff :: Vector{Int64} 
    confb :: Vector{Int64} 
    lid :: Vector{Int64}
    rid :: Vector{Int64}
end 


"""
    SConfs(nof :: Int64, nob :: Int64, nebm :: Int64, secd :: Vector{Int64}, qnd :: Vector{SQNDiag} ; , num_th :: Int64, disp_std :: Bool) :: Confs

generates the configurations from the list of QNDiags. 

# Arguments

* `nof :: Int64` is the number of fermionic sites ``N_{of}``.
* `nob :: Int64` is the number of bosonic sites ``N_{ob}``.
* `nebm :: Int64` is the maximal number of total bosons.
* `secd :: Vector{Int64}` is the set of ``Q_i`` for the selected configurations in the sector.
* `qnd :: Vector{SQNDiag}` is the set of [SQNDiags](@ref SQNDiag).
* `norf :: Int64` and `norb :: Int64` are the number of less significant bits used to generate the Lin table. Facultative, ``N_{of}/2`` and ``N_{ob}/2`` by default.
* `num_th :: Int64`, the number of threads. Facultative, `NumThreads` by default.
* `disp_std :: Bool`, whether or not the log shall be displayed. Facultative, `!SilentStd` by default. 

# Output

* `cfs :: SConfs` is a [SConfs](@ref SConfs) object.

# Note 

If your `qnd` has negative entries, QNDiags must contain the total number of particles (_i. e._, bosons plus fermions).
"""
function SConfs(nof :: Int64, nob :: Int64, nebm :: Int64, secd :: Vector{Int64}, qnd :: Vector{SQNDiag} ; norf :: Int64 = nof ÷ 2, norb :: Int64 = nob ÷ 2, num_th :: Int64 = NumThreads, disp_std :: Bool = !SilentStd)
    if nof == 0
        qnd1 = PadSQNDiag.(qnd, Ref(1), Ref(0), Ref(0), Ref(0))
        push!(qnd1, SQNDiag("Nef", [1], fill(0, nob)))
        secd1 = [secd ; 0]
        return SConfs(1, nob, nebm, secd1, qnd1 ; norb, num_th, disp_std)
    end

    nqnd = length(secd)
    binom = [ binomial(i + j, i) for i = 0 : nebm, j = 0 : nob]
    lid = Vector{Int64}(undef, 2 ^ (nof - norf) * (binom[nebm + 1, nob - norb + 1] + 1) + 1)
    ref_ncf = Ref{Int64}(0)
    secd1 = Int64[]
    qndf1 = Vector{Int64}[]
    qndb1 = Vector{Int64}[]
    modul = [ qnd_i.modul for qnd_i in qnd]

    # Check positivity 
    secp = 0
    qndfp = zeros(Int64, nof)
    qndbp = zeros(Int64, nob)
    negf = zeros(Int64, nof)
    negb = zeros(Int64, nob)
    for i = 1 : nqnd 
        (modul[i] > 1) && continue
        if (all(>=(0), [qnd[i].chargef ; qnd[i].chargeb]))
            secp += secd[i]
            qndfp += qnd[i].chargef
            qndbp += qnd[i].chargeb
        else
            for o = 1 : nof
                if (qnd[i].chargef[o] < 0) negf[o] = 1 end
            end
            for o = 1 : nob
                if (qnd[i].chargeb[o] < 0) negb[o] = 1 end
            end
        end
    end
    posq = true 
    for o = 1 : nof
        if (negf[o] == 1 && qndfp[o] == 0) posq = false end
    end
    for o = 1 : nob
        if (negb[o] == 1 && qndbp[o] == 0) posq = false end
    end

    posq || @error "Error"

    # Turn positive
    for i = 1 : nqnd 
        if (modul[i] > 1 || all(>=(0), [qnd[i].chargef ; qnd[i].chargeb])) 
            push!(secd1, secd[i]) 
            push!(qndf1, qnd[i].chargef)
            push!(qndb1, qnd[i].chargeb)
            continue
        end
        qm = minimum([ [ fld(qnd[i].chargef[o], qndfp[o]) for o = 1 : nof if qndfp[o] ≠ 0 ] ;
            [ fld(qnd[i].chargeb[o], qndbp[o]) for o = 1 : nob if qndbp[o] ≠ 0 ]])
        push!(secd1, secd[i] - qm * secp)
        push!(qndf1, qnd[i].chargef .- qm .* qndfp)
        push!(qndb1, qnd[i].chargeb .- qm .* qndbp)
    end
    qndf1_mat = hcat(qndf1...)
    qndb1_mat = hcat(qndb1...)

    @ccall Libpathino.__scfs_MOD_count_scfs(
        nof :: Ref{Int64}, nob :: Ref{Int64}, norf :: Ref{Int64}, norb :: Ref{Int64}, nebm :: Ref{Int64},
        nqnd :: Ref{Int64}, secd1 :: Ref{Int64}, qndf1_mat :: Ref{Int64}, qndb1_mat :: Ref{Int64}, modul :: Ref{Int64}, 
        ref_ncf :: Ref{Int64}, lid :: Ref{Int64}, binom :: Ref{Int64},
        num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}
    ) :: Nothing
    ncf = ref_ncf[]
    rid = Vector{Int64}(undef, 2 ^ norf * (binom[nebm + 1, norb + 1] + 1) + 1)
    conff = Vector{Int64}(undef, ncf)
    confb = Vector{Int64}(undef, ncf)
    @ccall Libpathino.__scfs_MOD_generate_scfs(
        nof :: Ref{Int64}, nob :: Ref{Int64}, norf :: Ref{Int64}, norb :: Ref{Int64}, nebm :: Ref{Int64},
        nqnd :: Ref{Int64}, secd1 :: Ref{Int64}, qndf1_mat :: Ref{Int64}, qndb1_mat :: Ref{Int64}, modul :: Ref{Int64}, 
        ncf :: Ref{Int64}, lid :: Ref{Int64}, rid :: Ref{Int64}, 
        conff :: Ref{Int64}, confb :: Ref{Int64}, binom :: Ref{Int64},
        num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}
    ) :: Nothing
        
    return SConfs(nof, nob, norf, norb, nebm, ncf, conff, confb, lid, rid)
end


"""
    EncodeConfb(nocc_b :: AbstractVector{Int64}, nebm :: Int64) :: Int64

encodes the occupation `nocc_b[o], o = 1 : nob` of the bosonic sites into the number `confb` used by `SConfs`. It reproduces the combinatorial index `encode_nb` of `scfs.f90`,
```math
    \\mathtt{confb}=∑_a C^{n^{(a)}_{eb}+o_a-1}_{n^{(a)}_{eb}}
```
where the bosons ``a=1,…,N_{eb}`` are taken from the last site to the first, ``o_a`` is the site of the ``a``-th boson and ``n^{(a)}_{eb}=\\mathtt{nebm}-a+1``.

_N.b._, the doc-string of [SConfs](@ref) still describes the older bit-string encoding ; the binary of `FuzzifiED_jll` in use implements the combinatorial index above.
"""
function EncodeConfb(nocc_b :: AbstractVector{Int64}, nebm :: Int64)
    (sum(nocc_b) ≤ nebm) || error("the number of bosons exceeds nebm")
    cfb = 0
    pos = length(nocc_b)
    neb1 = nebm
    nb1 = collect(nocc_b)
    while pos > 0
        if (nb1[pos] == 0)
            pos -= 1
            continue
        end
        cfb += binomial(neb1 + pos - 1, neb1)
        nb1[pos] -= 1
        neb1 -= 1
    end
    return cfb
end


"""
    DecodeConfb(cfb :: Int64, nob :: Int64, nebm :: Int64) :: Vector{Int64}

decodes the number `confb` used by `SConfs` back into the occupation `nocc_b[o], o = 1 : nob` of the bosonic sites, inverting [EncodeConfb](@ref). It reproduces `decode_nb` of `scfs.f90` : the largest ``C^{n_{eb}+o-1}_{n_{eb}}`` that fits into the remaining index is subtracted, which places one boson at the site ``o``.
"""
function DecodeConfb(cfb :: Int64, nob :: Int64, nebm :: Int64)
    nocc_b = zeros(Int64, nob)
    pos = nob
    neb1 = nebm
    cfb1 = cfb
    while cfb1 > 0
        (pos == 0 || neb1 == 0) && error("$cfb is not a valid confb for nob = $nob and nebm = $nebm")
        if (cfb1 < binomial(neb1 + pos - 1, neb1))
            pos -= 1
            continue
        end
        cfb1 -= binomial(neb1 + pos - 1, neb1)
        nocc_b[pos] += 1
        neb1 -= 1
    end
    return nocc_b
end


"""
    GetSConfId(cfs :: SConfs, cff :: Int64, nocc_b :: AbstractVector{Int64}) :: Int64
    GetSConfId(cfs :: SConfs, cff :: Int64, cfb :: Int64) :: Int64

inversely look up the index from the configuration, in analogy with [GetConfId](@ref) for `Confs`. It reproduces the Lin-table search `search_scf` of `scfs.f90` : the configuration is split into the `norf` least significant fermionic sites together with the bosonic sites `1 : norb` on the one hand, and the remaining fermionic sites and the bosonic sites `norb + 1 : nob` on the other,
```math
    \\mathtt{id}=\\mathtt{lid}[l]+\\mathtt{rid}[r]
```
```math
    l=2^{N_{of}-\\mathtt{norf}}\\,\\mathtt{cfb}_{>\\mathtt{norb}}+⌊\\mathtt{cff}/2^{\\mathtt{norf}}⌋+1,\\qquad
    r=2^{\\mathtt{norf}}\\,\\mathtt{cfb}_{≤\\mathtt{norb}}+(\\mathtt{cff}\\bmod 2^{\\mathtt{norf}})+1
```
where ``\\mathtt{cfb}_{≤\\mathtt{norb}}`` and ``\\mathtt{cfb}_{>\\mathtt{norb}}`` are the indices of [EncodeConfb](@ref) restricted to the two halves of the bosonic sites, both taken with the same `nebm`.

# Arguments

* `cfs :: SConfs` stores the configurations.
* `cff :: Int64` is the fermionic configuration to be looked-up expressed in a binary number, as in [GetConfId](@ref).
* `nocc_b :: AbstractVector{Int64}` is the occupation `nocc_b[o], o = 1 : cfs.nob` of the bosonic sites ; alternatively `cfb :: Int64`, the bosonic configuration already encoded by [EncodeConfb](@ref), which is then decoded by [DecodeConfb](@ref).

# Output

* `id :: Int64` is the id of the configuration such that `cfs.conff[id] == cff` and `cfs.confb[id] == cfb`. As for [GetConfId](@ref), the Lin table does not know whether the configuration belongs to the sector, so a configuration that is absent from `cfs` returns a meaningless index and it is up to the caller to check `cfs.conff[id]` and `cfs.confb[id]`.
"""
function GetSConfId(cfs :: SConfs, cff :: Int64, nocc_b :: AbstractVector{Int64})
    (length(nocc_b) == cfs.nob) || error("nocc_b must be an occupation vector of length cfs.nob")
    id_l = (EncodeConfb(view(nocc_b, cfs.norb + 1 : cfs.nob), cfs.nebm) << (cfs.nof - cfs.norf)) + (cff >> cfs.norf) + 1
    id_r = (EncodeConfb(view(nocc_b, 1 : cfs.norb), cfs.nebm) << cfs.norf) + (cff & (1 << cfs.norf - 1)) + 1
    return cfs.lid[id_l] + cfs.rid[id_r]
end
GetSConfId(cfs :: SConfs, cff :: Int64, cfb :: Int64) = GetSConfId(cfs, cff, DecodeConfb(cfb, cfs.nob, cfs.nebm))