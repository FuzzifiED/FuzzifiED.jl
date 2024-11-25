"""
    mutable struct SConfs

This type stores all the configurations that respects the diagonal quantum numbers (SQNDiag) and also a table to inversely look up the index from the configuration. 

# Fields

* `nof :: Int64` is the number of fermionic orbital\\*flavour.
* `nof :: Int64` is the number of bosonic orbital\\*flavour.
* `ncf :: Int64` is the number of configurations.
* `conff :: Vector{Int64}` is an array of length `ncf` containing all the fermion configurations. Each configuration is expressed in a binary number. If the `o-1`-th bit of `conf[i]` is 1, then the `o`-th orbital in the `i`-th configuration is occupied ; if the bit is 0, then the orbital is empty. 
* `confb :: Vector{Int64}` is an array of length `ncf` containing all the boson configurations. Each configuration is expressed in a binary number that has ``N_{b,o}`` 1's and ``N_{eb}`` 0's and the number of 0's following each 1 records the number of bosons in that orbital. 
* `lid :: Vector{Int64}` and `rid :: Vector{Int64}` contain the information of Lin table that is used to inversely look up the index `i` from the configuration. 
"""
mutable struct SConfs
    nof :: Int64 
    nob :: Int64
    ncf :: Int64
    conff :: Vector{Int64} 
    confb :: Vector{Int64} 
    lid :: Vector{Int64}
    rid :: Vector{Int64}
end 


"""
    function SConfs(nof :: Int64, nob :: Int64, nebm :: Int64, secd :: Vector{Int64}, qnd :: Vector{SQNDiag} ; , num_th :: Int64, disp_std :: Bool) :: Confs

generates the configurations from the list of QNDiags. 

# Arguments

* `nof :: Int64` is the number of fermionic orbital\\*flavour ``N_{of}``.
* `nob :: Int64` is the number of bosonic orbital\\*flavour ``N_{on}``.
* `nebm :: Int64` is the maximal number of total bosons.
* `secd :: Vector{Int64}` is the set of ``Q_i`` for the selected configurations in the sector.
* `qnd :: Vector{SQNDiag}` is the set of [SQNDiags](@ref SQNDiag). 
* `num_th :: Int64`, the number of threads. Facultative, `NumThreads` by default. 
* `disp_std :: Bool`, whether or not the log shall be displayed. Facultative, `!SilentStd` by default. 

# Output

* `cfs :: SConfs` is a [SConfs](@ref SConfs) object.

# Note 

If your `qnd` has negative entries, QNDiags must contain the total number of particles (_i.e._, bosons plus fermions).
"""
function SConfs(nof :: Int64, nob :: Int64, nebm :: Int64, secd :: Vector{Int64}, qnd :: Vector{SQNDiag} ; num_th :: Int64 = NumThreads, disp_std :: Bool = !SilentStd)
    nqnd = length(secd)
    lid = Vector{Int64}(undef, 2 ^ nof + 1)
    ref_ncf = Ref{Int64}(0)
    secd1 = Int64[]
    qndf1 = Vector{Int64}[]
    qndb1 = Vector{Int64}[]
    modul = [ qnd_i.modul for qnd_i in qnd]
    id_ne = 0
    for i = 1 : nqnd
        if (maximum([qnd[i].chargef ; qnd[i].chargeb]) == 1 && minimum([qnd[i].chargef ; qnd[i].chargeb]) == 1)
            id_ne = i
            break
        end 
    end
    for i = 1 : nqnd
        if (modul[i] > 1 || (minimum(qnd[i].chargef) ≥ 0 && minimum(qnd[i].chargeb) ≥ 0)) 
            push!(secd1, secd[i]) 
            push!(qndf1, qnd[i].chargef)
            push!(qndb1, qnd[i].chargeb)
            continue
        end
        qm = minimum([qnd[i].chargef ; qnd[i].chargeb])
        push!(secd1, secd[i] .- secd[id_ne] * qm)
        push!(qndf1, qnd[i].chargef .- qm)
        push!(qndb1, qnd[i].chargeb .- qm)
    end
    qndf1_mat = Matrix{Int64}(reduce(hcat, qndf1))
    qndb1_mat = Matrix{Int64}(reduce(hcat, qndb1))

    @ccall Libpathino.__scfs_MOD_count_scfs(
        nof :: Ref{Int64}, nob :: Ref{Int64}, 
        nqnd :: Ref{Int64}, secd1 :: Ref{Int64}, 
        qndf1_mat :: Ref{Int64}, qndb1_mat :: Ref{Int64}, modul :: Ref{Int64}, 
        ref_ncf :: Ref{Int64}, lid :: Ref{Int64}, 
        num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}
    ) :: Nothing
    ncf = ref_ncf[]
    rid = Vector{Int64}(undef, 2 ^ (nob + nebm))
    conff = Vector{Int64}(undef, ncf)
    confb = Vector{Int64}(undef, ncf)
    @ccall Libpathino.__scfs_MOD_generate_scfs(
        nof :: Ref{Int64}, nob :: Ref{Int64}, nebm :: Ref{Int64},
        nqnd :: Ref{Int64}, secd1 :: Ref{Int64}, 
        qndf1_mat :: Ref{Int64}, qndb1_mat :: Ref{Int64}, modul :: Ref{Int64}, 
        ncf :: Ref{Int64}, lid :: Ref{Int64}, rid :: Ref{Int64}, 
        conff :: Ref{Int64}, confb :: Ref{Int64}, 
        num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}
    ) :: Nothing
    return SConfs(nof, nob, ncf, conff, confb, lid, rid)
end