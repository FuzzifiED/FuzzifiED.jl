export Term, Terms
export ParticleHole, NormalOrder, SimplifyTerms


"""
    Term

The mutable type `Term` records a term that looks like ``Uc^{(p_1)}_{o_1}c^{(p_2)}_{o_2}… c^{(p_l)}_{o_l}`` in an operator

# Fields

* `coeff :: ComplexF64` records the coefficient ``U``.
* `cstr :: Vector{Int64}` is a length-``2l`` vector ``(p_1,o_1,p_2,o_2,… p_l,o_l)`` recording the operator string.

# Initialisation

It can be generated by the function

    Term(coeff :: Number, cstr :: Vector{Int64})
"""
mutable struct Term 
    coeff :: ComplexF64 
    cstr :: Vector{Int64}
end 


""" 
    Terms

`Terms` is an alias for `Vector{Term}` for convenience

# Initialisation

    Terms(coeff :: Number, cstr :: Vector{Int64})

Gives a `Terms` with a single `Term`.

# Special elements

The zero and identity terms are defined.

    zero(Terms) = Term[]
    one(Terms) = Terms(1, [-1, -1])
"""
const Terms = Vector{Term}
Vector{T}(coeff :: Number, cstr :: Vector{Int64}) where T <: Term = [Term(coeff, cstr)]
Base.zero( :: Type{Terms}) = Term[]
Base.one( :: Type{Terms}) = Terms(1, [-1, -1])


"""
    *(fac :: Number, tms :: Terms) :: Terms
    -(tms :: Terms) :: Terms
    *(tms :: Terms, fac :: Number) :: Terms
    /(tms :: Terms, fac :: Number) :: Terms

Return the product of a collection of terms with a number. 
"""
function Base.:*(fac :: Number, tms :: Terms)
    return [ Term(fac * tm.coeff, tm.cstr) for tm in tms ]
end
function Base.:-(tms :: Terms)
    return (-1) * tms
end
function Base.:*(tms :: Terms, fac :: Number)
    return fac * tms
end
function Base.:/(tms :: Terms, fac :: Number)
    return (1 / fac) * tms
end


"""
    +(tms1 :: Terms, tms2 :: Terms) :: Terms
    -(tms1 :: Terms, tms2 :: Terms) :: Terms

Return the naive sum of two series of terms by taking their union. 
"""
function Base.:+(tms1 :: Terms, tms2 :: Terms)
    return [ tms1 ; tms2 ]
end
function Base.:-(tms1 :: Terms, tms2 :: Terms)
    return tms1 + (-tms2)
end
function Base.:+(tms1 :: Terms, tms2 :: Vararg{Terms})
    return tms1 + +(tms2...)
end


"""
    *(tms1 :: Terms, tms2 :: Terms) :: Terms
    ^(tms :: Terms, pow :: Int64) :: Terms

Return the naive product of two series of terms or the power of one terms. The number of terms equals the product of the number of terms in `tms1` and `tms2`. For each term in `tms1` ``Uc^{(p_1)}_{o_1}…`` and `tms2` ``U'c^{(p'_1)}_{o'_1}…``, a new term is formed by taking ``UU'c^{(p_1)}_{o_1}… c^{(p'_1)}_{o'_1}…``.
"""
function Base.:*(tms1 :: Terms, tms2 :: Terms)
    return Terms(vcat([ Term(tm1.coeff * tm2.coeff, [tm1.cstr ; tm2.cstr])
        for tm1 in tms1, tm2 in tms2 ]...))
end
function Base.:*(tms1 :: Terms, tms2 :: Vararg{Terms})
    return tms1 * *(tms2...)
end
function Base.:^(tms :: Terms, pow :: Int64)
    if pow == 1
        return tms
    else
        return tms * tms ^ (pow - 1)
    end
end

function Base.adjoint(tm :: Term)
    nc = length(tm.cstr)
    cstr1 = [ isodd(i) ? 1 - tm.cstr[nc - i] : tm.cstr[nc + 2 - i] for i = 1 : nc]
    return Term(conj(tm.coeff), cstr1)
end
"""
    adjoint(tm :: Term) :: Term
    adjoint(tms :: Terms) :: Terms

Return the Hermitian conjugate of a series of terms. For each term ``Uc^{(p_1)}_{o_1}c^{(p_2)}_{o_2}… c^{(p_l)}_{o_l}``, the adjoint is ``\\bar{U}c^{(1-p_l)}_{o_l}… c^{(1-p_2)}_{o_2}c^{(1-p_1)}_{o_1}``.
"""
function Base.adjoint(tms :: Terms)
    return adjoint.(tms)
end

function ParticleHole(tm :: Term)
    nc = length(tm.cstr)
    cstr1 = [ isodd(i) ? 1 - tm.cstr[i] : tm.cstr[i] for i = 1 : nc]
    return Term(tm.coeff, cstr1)
end
"""
    ParticleHole(tm :: Term) :: Term
    ParticleHole(tms :: Terms) :: Terms

Return the particle-hole transformation of a series of terms. For each term ``Uc^{(p_1)}_{o_1}c^{(p_2)}_{o_2}… c^{(p_l)}_{o_l}``, the transformation results in ``Uc^{(1-p_1)}_{o_1}c^{(1-p_2)}_{o_2}…c^{(1-p_l)}_{o_l}``.
"""
function ParticleHole(tms :: Terms)
    return ParticleHole.(tms)
end


"""
    NormalOrder(tm :: Term) :: Terms

rearrange a term such that 
* the creation operators must be commuted in front of the annihilation operator 
* the site index of the creation operators are in ascending order and the annihilation operators in descending order. 
return a list of terms whose result is equal to the original term. 
"""
function NormalOrder(tm :: Term)

    coeff0 = tm.coeff
    cstr0 = tm.cstr
    for i = 1 : 2 : length(cstr0) - 3
        if (cstr0[i] == -1) 
            cstr1 = deepcopy(cstr0)
            deleteat!(cstr1, i : i + 1)
            return(NormalOrder(Term(coeff0, cstr1)))
        end
        if (cstr0[i] == 0 && cstr0[i + 2] == 1)
            if (cstr0[i + 1] == cstr0[i + 3])
                cstr_nrm = deepcopy(cstr0)
                cstr_com = deepcopy(cstr0)
                cstr_nrm[i : i + 1], cstr_nrm[i + 2 : i + 3] = cstr_nrm[i + 2 : i + 3], cstr_nrm[i : i + 1]
                deleteat!(cstr_com, i : i + 3)
                return([ NormalOrder(Term(-coeff0, cstr_nrm)) ; 
                    NormalOrder(Term(coeff0, cstr_com))])
            else
                cstr_nrm = deepcopy(cstr0)
                cstr_nrm[i : i + 1], cstr_nrm[i + 2 : i + 3] = cstr_nrm[i + 2 : i + 3], cstr_nrm[i : i + 1]
                return(NormalOrder(Term(-coeff0, cstr_nrm)))
            end
        elseif (cstr0[i] == cstr0[i + 2])
            if (cstr0[i + 1] == cstr0[i + 3]) return Term[] end 
            if ((cstr0[i] == 1) == (cstr0[i + 1] > cstr0[i + 3]))
                cstr_nrm = deepcopy(cstr0)
                cstr_nrm[i : i + 1], cstr_nrm[i + 2 : i + 3] = cstr_nrm[i + 2 : i + 3], cstr_nrm[i : i + 1]
                return(NormalOrder(Term(-coeff0, cstr_nrm)))
            end
        end
    end
    if length(cstr0) == 0 return Terms(coeff0, [-1, -1]) end
    if (cstr0[end - 1] == -1 && length(cstr0) > 2) return Terms(coeff0, cstr0[1 : end - 2]) end
    return Term[tm]
end


"""
    SimplifyTerms(tms :: Terms ; cutoff :: Float64 = eps(Float64)) :: Terms

simplifies the sum of terms such that 
* each term is normal ordered,
* like terms are combined, and terms with zero coefficients are removed.

# Argument 

* `cutoff :: Float64` is the cutoff such that terms with smaller absolute value of coefficients will be neglected. Facultative, `eps(Float64)` by default. 
"""
function SimplifyTerms(tms :: Terms ; cutoff :: Float64 = eps(Float64)) :: Terms
    dictlock = [ ReentrantLock() for i = 1 : 64 ]
    dict_tms = [ Dict{Vector{Int64}, ComplexF64}() for i = 1 : 64 ]
    
    Threads.@threads for tm in tms 
        tm1 = NormalOrder(tm)
        for tmi in tm1 
            id = 1 + tmi.cstr[2] & 63
            Threads.lock(dictlock[id]) 
            try
                if haskey(dict_tms[id], tmi.cstr) 
                    dict_tms[id][tmi.cstr] += tmi.coeff
                else
                    dict_tms[id][tmi.cstr] = tmi.coeff
                end
            finally
                Threads.unlock(dictlock[id])
            end
        end
    end
    tms_f = Term[ 
        Term(coeff_i, cstr_i)
        for dict_tms_i in dict_tms 
        for (cstr_i, coeff_i) in dict_tms_i if abs(coeff_i) > cutoff
    ]
    return tms_f
end
