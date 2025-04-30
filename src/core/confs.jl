export Confs, GetConfId 


"""
    Confs

The mutable type `Confs` stores all the configurations that respects the diagonal quantum numbers (QNDiag) and also a table to inversely look up the index from the configuration. 

# Fields

* `no :: Int64` is the number of sites.
* `ncf :: Int64` is the number of configurations.
* `conf :: Vector{Int64}` is an array of length `ncf` containing all the configurations. Each configuration is expressed in a binary number. If the `o-1`-th bit of `conf[i]` is 1, then the `o`-th site in the `i`-th configuration is occupied ; if the bit is 0, then the site is empty. 
* `nor :: Int64`, `lid :: Vector{Int64}` and `rid :: Vector{Int64}` contain the information of Lin table that is used to inversely look up the index `i` from the configuration. 
"""
mutable struct Confs
    no :: Int64 
    nor :: Int64
    ncf :: Int64
    conf :: Vector{Int64} 
    lid :: Vector{Int64}
    rid :: Vector{Int64}
end 


"""
    Confs(no :: Int64, secd :: Vector{Int64}, qnd :: Vector{QNDiag} ; nor :: Int64 = div(no, 2), num_th :: Int64, disp_std :: Bool)

generates the configurations from the list of QNDiags. 

# Arguments

* `no :: Int64` is the number of sites ``N_o``.
* `secd :: Vector{Int64}` is the set of ``Q_i`` for the selected configurations in the sector.
* `qnd :: Vector{QNDiag}` is the set of [QNDiags](@ref QNDiag).
* `nor :: Int64` is the number of less significant bits used to generate the Lin table. Facultative, ``N_o/2`` by default.
* `num_th :: Int64`, the number of threads. Facultative, `NumThreads` by default.
* `disp_std :: Bool`, whether or not the log shall be displayed. Facultative, `!SilentStd` by default. 

# Output

* `cfs :: Confs` is a [Confs](@ref Confs) object.

# Note 

If your `qnd` has negative entries, QNDiags must contain the number of electrons.
"""
function Confs(no :: Int64, secd :: Vector{Int64}, qnd :: Vector{QNDiag} ; nor :: Int64 = div(no, 2), num_th :: Int64 = NumThreads, disp_std :: Bool = !SilentStd)
    nqnd = length(secd)
    lid = Vector{Int64}(undef, 2 ^ (no - nor) + 1)
    ref_ncf = Ref{Int64}(0)
    secd1 = Int64[]
    qnd1 = Vector{Int64}[]
    modul = [ qnd_i.modul for qnd_i in qnd]
    
    # Check positivity 
    secp = 0
    qndp = zeros(Int64, no)
    neg = zeros(Int64, no)
    for i = 1 : nqnd 
        (modul[i] > 1) && continue
        if (all(>=(0), qnd[i].charge))
            secp += secd[i]
            qndp += qnd[i].charge
        else
            for o = 1 : no 
                if (qnd[i].charge[o] < 0) neg[o] = 1 end
            end
        end
    end
    posq = true 
    for o = 1 : no 
        if (neg[o] == 1 && qndp[o] == 0) 
            posq = false
        end
    end


    if (posq)
        # Turn positive
        for i = 1 : nqnd 
            if (modul[i] > 1 || all(>=(0), qnd[i].charge)) 
                push!(secd1, secd[i]) 
                push!(qnd1, qnd[i].charge)
                continue
            end
            qm = minimum([ fld(qnd[i].charge[o], qndp[o]) for o = 1 : no if qndp[o] ≠ 0 ])
            push!(secd1, secd[i] - qm * secp)
            push!(qnd1, qnd[i].charge .- qm .* qndp)
        end
        qnd1_mat = hcat(qnd1...)

        @ccall Libpath.__cfs_MOD_count_cfs(
            no :: Ref{Int64}, nor :: Ref{Int64}, 
            nqnd :: Ref{Int64}, secd1 :: Ref{Int64}, 
            qnd1_mat :: Ref{Int64}, modul :: Ref{Int64}, 
            ref_ncf :: Ref{Int64}, lid :: Ref{Int64}, 
            num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}
        ) :: Nothing
        ncf = ref_ncf[]
        rid = Vector{Int64}(undef, 2 ^ nor)
        conf = Vector{Int64}(undef, ncf)
        @ccall Libpath.__cfs_MOD_generate_cfs(
            no :: Ref{Int64}, nor :: Ref{Int64}, 
            nqnd :: Ref{Int64}, secd1 :: Ref{Int64}, 
            qnd1_mat :: Ref{Int64}, modul :: Ref{Int64}, 
            ncf :: Ref{Int64}, lid :: Ref{Int64}, rid :: Ref{Int64}, conf :: Ref{Int64}, 
            num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}
        ) :: Nothing
    else 
        # Not turn positive
        secd1 = secd
        qnd1_mat = [qnd[i].charge[o] for o = 1 : no, i = 1 : nqnd]

        @ccall Libpath.__cfs_neg_MOD_count_cfs_neg(
            no :: Ref{Int64}, nor :: Ref{Int64}, 
            nqnd :: Ref{Int64}, secd1 :: Ref{Int64}, 
            qnd1_mat :: Ref{Int64}, modul :: Ref{Int64}, 
            ref_ncf :: Ref{Int64}, lid :: Ref{Int64}, 
            num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}
        ) :: Nothing
        ncf = ref_ncf[]
        rid = Vector{Int64}(undef, 2 ^ nor)
        conf = Vector{Int64}(undef, ncf)
        @ccall Libpath.__cfs_neg_MOD_generate_cfs_neg(
            no :: Ref{Int64}, nor :: Ref{Int64}, 
            nqnd :: Ref{Int64}, secd1 :: Ref{Int64}, 
            qnd1_mat :: Ref{Int64}, modul :: Ref{Int64}, 
            ncf :: Ref{Int64}, lid :: Ref{Int64}, rid :: Ref{Int64}, conf :: Ref{Int64}, 
            num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}
        ) :: Nothing
    end

    return Confs(no, nor, ncf, conf, lid, rid)
end


"""
    GetConfId(cfs :: Confs, cf :: Int64) :: Int64

inversely look up the index from the configuration

# Arguments

* `cfs :: Confs` stores the configurations.
* `cf :: Int64` stores the configuration to be looked-up expressed in a binary number. If the `o-1`-th bit of `conf[i]` is 1, then the `o`-th site in the `i`-th configuration is occupied ; if the bit is 0, then the site is empty. 

# Output
* `id :: Int64` is the id of the configuration such that `cfs.conf[id] == cf`.
"""
function GetConfId(cfs :: Confs, cf :: Int64)
    cf_l = cf >> cfs.nor
    cf_r = cf & (1 << cfs.nor - 1)
    return cfs.lid[cf_l + 1] + cfs.rid[cf_r + 1]
end 
