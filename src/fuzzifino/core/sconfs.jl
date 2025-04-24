export SConfs


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

If your `qnd` has negative entries, QNDiags must contain the total number of particles (_i.e._, bosons plus fermions).
"""
function SConfs(nof :: Int64, nob :: Int64, nebm :: Int64, secd :: Vector{Int64}, qnd :: Vector{SQNDiag} ; norf :: Int64 = nof ÷ 2, norb :: Int64 = nob ÷ 2, num_th :: Int64 = NumThreads, disp_std :: Bool = !SilentStd)
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
    rid = Vector{Int64}(undef, 2 ^ norf * (binom[nebm + 1, nob - norb + 1] + 1))
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
