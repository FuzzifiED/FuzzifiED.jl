import FuzzifiED: GetSites, GetQNDiags, TruncateQNDiag, Confs, Terms


"""
    GetQNDiags(sites :: Vector{<:Index})

Converts a `Sites` object in the `ITensors` package to a set of QNDiags. 

# Arguments 

* `sites :: Vector{<:Index}` is a `Sites` object. Only `Fermion` site type is supported, and the quantum numbers of the `0` state must be all zero. Note that this will subject to the limitation in ITensors that the number of conserved quantities must be less than 4.
"""
function GetQNDiags(sites :: Vector{<:Index})
    names = []
    moduls = Int64[]
    for site in sites 
        for qn in qn(site, 2)
            if abs(qn.modulus) == 0 continue end
            if qn.name in names continue end
            push!(names, qn.name)
            push!(moduls, abs(qn.modulus))
        end
    end
    qnd = QNDiag[]
    for i in eachindex(names)
        charge = Int64[]
        for site in sites 
            push!(charge, val(qn(site, 2), names[i]))
        end
        push!(qnd, QNDiag(String(names[i]), charge, moduls[i]))
    end
    return qnd
end 


"""
    Confs(sites :: Vector{<:Index}, sec_qn :: QN)
    Confs(sites :: Vector{<:Index}, cf_ref :: Vector{Int64})

Converts a `Sites` object in the `ITensors` package to the `Confs` object

# Arguments 

* `sites :: Vector{<:Index}` is a `Sites` object. Only `Fermion` site type is supported, and the quantum numbers of the `0` state must be all zero. Note that this will subject to the limitation in ITensors that the number of conserved quantities must be less than 4. 
* `sec_qn :: QN` is a `QN` object that specifies the the quantum number of the selected configuration. Alternatively, `cf_ref :: Vector{Int64}` is a reference configuration composed of `0` and `1`.
"""
function Confs(sites :: Vector{<:Index}, sec_qn :: QN)
    no = length(sites)
    qnd = GetQNDiags(sites)
    secd = [ val(sec_qn, qndi.name) for qndi in qnd ]
    return Confs(no, secd, qnd)
end 
function Confs(sites :: Vector{<:Index}, cf_ref :: Vector{Int64})
    no = length(sites)
    qnd = GetQNDiags(sites)
    secd = [ cf_ref' * qndi.charge for qndi in qnd ]
    return Confs(no, secd, qnd)
end 


"""
    TruncateQNDiag(qnd :: Vector{QNDiag} ; trunc_lth :: Int64, trunc_wt :: Vector{Int64}) :: Vector{QNDiag}

truncates the list of ``N_U`` QNDiags from to a number ``N'_U`` acceptable by ITensors. The new quantum numbers are 
```math
\\begin{aligned}
    &Q'_1=Q_1,\\ Q'_2=Q_2,\\ …,\\ Q'_{N'_U-1}=Q_{N'_U-1}\\\\
    &Q'_{N'_U}=λ_{N'_U}Q_{N'_U}+λ_{N'_U+1}Q_{N'_U+1}+…+λ_{N_U}Q_{N_U}
\\end{aligned}
```

# Arguments

* `qnd :: Vector{QNDiag}` stores the set of QNDiags. 
* `trunc_lth :: Int64` stores the truncated numbers of QNDiags. Facultative, 3 by default. 
* `trunc_wt :: Vecotr{Int64}` stores the ``N_U-N'_U+1`` coefficients ``λ``. Facultative, ``1,10,100,1000,…`` by default. 
"""
function TruncateQNDiag(qnd :: Vector{QNDiag} ; trunc_lth :: Int64 = 3, trunc_wt :: Vector{Int64} = [ 10 ^ (i - trunc_lth) for i = trunc_lth : length(qnd)]) 
    return [ qnd[1 : trunc_lth - 1]..., sum(qnd[trunc_lth : end]) .* trunc_wt ]
end


"""
    GetSites(qnd :: Vector{QNDiag}) :: Vector{<:Index}

returns the ITensors Sites of type "FuzzyFermion" from a set of QNDiags.
"""
function GetSites(qnd :: Vector{QNDiag})
    no = length(qnd[1].charge)
    return [ siteind("FuzzyFermion" ; o, qnd) for o = 1 : no ]
end


"""
    Terms(opsum :: OpSum)

Converts a `OpSum` object in `ITensors` to a series of terms. Note that the only operators supported are `"C"`, `"Cdag"` `"N"` and `"I"`.
"""
function Terms(opsum :: OpSum)
    tms = Term[]
    for op_prod in opsum
        cstr = []
        for op_i in ITensors.argument(op_prod)
            optype = op_i.which_op
            opsite = op_i.sites[1]
            if optype == "Cdag"
                append!(cstr, [1, opsite])
            elseif optype == "C"
                append!(cstr, [0, opsite])
            elseif optype == "N"
                append!(cstr, [1, opsite, 0, opsite])
            elseif optype == "I"
                append!(cstr, [-1, -1])
            else 
                print("The only operator supported are C, Cdag and N")
            end
        end 
        if length(ITensors.argument(op_prod)) == 0
            append!(cstr, [-1, -1])
        end
        fac = coefficient(op_prod)
        push!(tms, Term(fac, cstr))
    end
    return tms
end


"""
    OpSum(tms :: Terms)

Converts a series of terms to `OpSum` object in `ITensors`.
"""
function ITensors.Ops.OpSum(tms :: Terms)
    opsum = OpSum()
    for tm in tms
        optm = Any[tm.coeff]
        for i = 1 : 2 : length(tm.cstr)
            push!(optm, tm.cstr[i] == 0 ? "C" : "Cdag")
            push!(optm, tm.cstr[i + 1])
        end
        opsum += Tuple(optm)
    end
    return opsum
end
