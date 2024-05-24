using ITensors

function ConfsFromSite(sites :: Vector{Index{Vector{Pair{QN, Int64}}}}, qnu_s_it :: QN)
    no = length(sites)
    qn_names = []
    for site in sites 
        for qn in qn(site, 2)
            if abs(qn.modulus) != 1 continue end
            if qn.name in qn_names continue end
            push!(qn_names, qn.name)
        end
    end
    qnu_o = []
    for qn_name in qn_names 
        qnu_oi = []
        for site in sites 
            push!(qnu_oi, val(qn(site, 2), qn_name))
        end
        push!(qnu_o, qnu_oi)
    end
    qnu_s = [ val(qnu_s_it, qn_name) for qn_name in qn_names ]
    return Confs(no, qnu_s, qnu_o)
end 

function OperatorFromOpSum(bsd :: Basis, bsf :: Basis, red_q :: Int64, sym_q :: Int64, opsum :: Sum{Scaled{ComplexF64, Prod{Op}}})
    cstr_vec = []
    fac = Vector{ComplexF64}(undef, 0)
    for i = 1 : length(opsum)
        cstri = []
        for j = 1 : length(opsum[i])
            optype = opsum[i][j].which_op
            opsite = opsum[i][j].sites[1]
            if optype == "Cdag"
                append!(cstri, [1, opsite])
            elseif optype == "C"
                append!(cstri, [0, opsite])
            elseif optype == "N"
                append!(cstri, [1, opsite, 0, opsite])
            else 
                print("The only operator supported are C, Cdag and N")
            end
        end 
        if length(opsum[i]) == 0
            append!(cstri, [-1, -1])
        end
        push!(fac, coefficient(opsum[i]))
        push!(cstr_vec, cstri)
    end
    return Operator(bsd, bsf, red_q, sym_q, cstr_vec, fac)
end