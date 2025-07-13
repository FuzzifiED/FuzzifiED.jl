export GetIntMatrix, GetDenIntTerms, GetPairIntTerms, GetPolTerms, GetL2Terms, GetC2Terms


"""
    GetIntMatrix(nm :: Int64, ps_pot :: Vector{<:Number}) :: Array{ComplexF64, 3}

Gives the interaction matrix ``U_{m_1,m_2,m_3,m_4}`` from the pseudopotentials.
    
# Argument

* `nm :: Int64` is the number of orbitals.
* `ps_pot :: Vector{<:Number}` is the vector of non-zero pseudopotentials.

# Output
* A `nm`×`nm`×`nm` array giving the interaction matrix ``U_{m_1,m_2,m_3,-m_1-m_2-m_3}``.
"""
function GetIntMatrix(nm :: Int64, ps_pot :: Vector{<:Number})
    int_el = zeros(ComplexF64, nm, nm, nm)
    s = .5 * (nm - 1)
    for m1 in 1 : nm 
        m1r = m1 - s - 1
        for m2 in 1 : nm 
            m2r = m2 - s - 1
            for m3 = 1 : nm
                m3r = m3 - s - 1
                m4 = m1 + m2 - m3 
                if (m4 <= 0 || m4 > nm)
                    continue
                end
                m4r = m4 - s - 1
                for l in eachindex(ps_pot)
                    if (abs(m1r + m2r) > nm - l || abs(m3r + m4r) > nm - l)
                        break
                    end 
                    int_el[m1, m2, m3] += ps_pot[l] * (2 * nm - 2 * l + 1) * wigner3j(s, s, nm - l, m1r, m2r, -m1r - m2r) * wigner3j(s, s, nm - l, m4r, m3r, -m3r - m4r)
                end 
            end 
        end 
    end
    return int_el
end


"""
    GetDenIntTerms(nm :: Int64, nf :: Int64[, ps_pot :: Vector{<:Number}][, mat_a :: Matrix{<:Number}[, mat_b :: Matrix{<:Number}]][ ; m_kept :: Vector{Int64}]) :: Terms

Return the normal-ordered density-density term in the Hamiltonian 
```math 
∑_{\\{m_i,f_i\\}}U_{m_1m_2m_3m_4}M^A_{f_1f_4}M^B_{f_2f_3}c^{†}_{m_1f_1}c^{†}_{m_2f_2}c_{m_3f_3}c_{m_4f_4}.
```

# Arguments 

* `nm :: Int64` is the number of orbitals.
* `nf :: Int64` is the number of flavours.
* `ps_pot :: Vector{<:Number}` is a list of numbers specifying the pseudopotentials for the interacting matrix ``U_{m_1m_2m_3m_4}``. Facultative, `[1.0]` by default. 
* `mat_a :: Matrix{<:Number}` is a `nf`×`nf` matrix specifying ``M^A_{ff'}``. Facultative, ``I_{N_f}`` by default. 
* `mat_b :: Matrix{<:Number}` is a `nf`×`nf` matrix specifying ``M^B_{ff'}``. Facultative, the Hermitian conjugate of `mat_a` by default. 
* `m_kept :: Vector{Int64}` is a list of orbitals that range from 1 to `nm`. Facultative, if specified, only terms for which all ``m_i`` are in the list are kept. 
"""
function GetDenIntTerms(nm :: Int64, nf :: Int64, ps_pot :: Vector{<:Number}, mat_a :: Matrix{<:Number} = Matrix{Float64}(I, nf, nf), mat_b :: Matrix{<:Number} = Matrix(mat_a') ; m_kept :: Vector{Int64} = collect(1 : nm))
    no = nm * nf
    int_el = GetIntMatrix(nm, ps_pot)
    tms = Term[]
    for o1 = 1 : no
        m1 = div(o1 - 1, nf) + 1
        f1 = mod(o1 - 1, nf) + 1
        if (m1 ∉ m_kept) continue end
        for o2 = 1 : no 
            m2 = div(o2 - 1, nf) + 1
            f2 = mod(o2 - 1, nf) + 1
            if (o1 == o2) continue end
            if (m2 ∉ m_kept) continue end
            for o3 = 1 : no
                m3 = div(o3 - 1, nf) + 1
                f3 = mod(o3 - 1, nf) + 1
                if (m3 ∉ m_kept) continue end
                if (abs(mat_b[f2, f3]) < 1E-13) continue end 
                m4 = m1 + m2 - m3 
                if (m4 <= 0 || m4 > nm) continue end
                if (m4 ∉ m_kept) continue end
                for f4 = 1 : nf 
                    if (abs(mat_a[f1, f4]) < 1E-13) continue end
                    o4 = (m4 - 1) * nf + f4
                    if (o3 == o4) continue end
                    val = mat_a[f1, f4] * mat_b[f2, f3] * int_el[m1, m2, m3]
                    if (abs(val) < 1E-15) continue end 
                    push!(tms, Term(val, [1, o1, 1, o2, 0, o3, 0, o4]))
                end
            end
        end
    end
    return SimplifyTerms(tms)
end
GetDenIntTerms(nm :: Int64, nf :: Int64, mat_a :: Matrix{<:Number} = Matrix{Float64}(I, nf, nf), mat_b :: Matrix{<:Number} = Matrix(mat_a') ; m_kept :: Vector{Int64} = collect(1 : nm)) = GetDenIntTerms(nm, nf, [1.0], mat_a, mat_b ; m_kept)


"""
    GetDenIntTerms(nm :: Int64, nf :: Int64, ps_pot :: Vector{<:Number}, mat_a :: Vector{<:AbstractMatrix{<:Number}}[, mat_b :: Vector{<:AbstractMatrix{<:Number}}][ ; m_kept :: Vector{Int64}]) :: Terms

Return the sum of a series of normal-ordered density-density term in the Hamiltonian 
```math 
∑_{\\{m_i,f_i,α\\}}U_{m_1m_2m_3m_4}(M^A_{α})_{f_1f_4}(M^B_{α})_{f_2f_3}c^{†}_{m_1f_1}c^{†}_{m_2f_2}c_{m_3f_3}c_{m_4f_4}.
```

# Arguments 

* `nm :: Int64` is the number of orbitals.
* `nf :: Int64` is the number of flavours.
* `ps_pot :: Vector{<:Number}` is a list of numbers specifying the pseudopotentials for the interacting matrix ``U_{m_1m_2m_3m_4}``. Facultative, `[1.0]` by default.
* `mat_a :: Vector{<:AbstractMatrix{<:Number}}` is a vector of `nf`×`nf` matrix specifying ``(M^A_{α})_{ff'}``. Facultative, ``I_{N_f}`` by default. 
* `mat_b :: Vector{<:AbstractMatrix{<:Number}}` is a vector of `nf`×`nf` matrix specifying ``(M^B_{α})_{ff'}``. Facultative, the Hermitian conjugate of `mat_a` by default. 
* `m_kept :: Vector{Int64}` is a list of orbitals that range from 1 to `nm`. Facultative, if specified, only terms for which all ``m_i`` are in the list are kept. 
"""
function GetDenIntTerms(nm :: Int64, nf :: Int64, ps_pot :: Vector{<:Number}, mats_a :: Vector{<:AbstractMatrix{<:Number}}, mats_b :: Vector{<:AbstractMatrix{<:Number}} = [Matrix(mat_a') for mat_a in mats_a] ; m_kept :: Vector{Int64} = collect(1 : nm))
    return sum([GetDenIntTerms(nm, nf, ps_pot, mats_a[i], mats_b[i] ; m_kept) for i in eachindex(mats_a)])
end
GetDenIntTerms(nm :: Int64, nf :: Int64, mats_a :: Vector{<:AbstractMatrix{<:Number}}, mats_b :: Vector{<:AbstractMatrix{<:Number}} = [Matrix(mat_a') for mat_a in mats_a] ; m_kept :: Vector{Int64} = collect(1 : nm)) = GetDenIntTerms(nm, nf, [1.0], mats_a, mats_b ; m_kept)


"""
    GetPairIntTerms(nm :: Int64, nf :: Int64, ps_pot :: Vector{<:Number}, mat_a :: Matrix{<:Number}[, mat_b :: Matrix{<:Number}][ ; m_kept :: Vector{Int64}]) :: Terms

Return the normal-ordered pair-pair interaction term in the Hamiltonian 
```math 
∑_{\\{m_i,f_i\\}}U_{m_1m_2m_3m_4}M^A_{f_1f_2}M^B_{f_3f_4}c^{†}_{m_1f_1}c^{†}_{m_2f_2}c_{m_3f_3}c_{m_4f_4}.
```

# Arguments 

* `nm :: Int64` is the number of orbitals.
* `nf :: Int64` is the number of flavours.
* `ps_pot :: Vector{<:Number}` is a list of numbers specifying the pseudopotentials for the interacting matrix ``U_{m_1m_2m_3m_4}``. 
* `mat_a :: Matrix{<:Number}` is a `nf`×`nf` matrix specifying ``M^A_{ff'}``. Facultative, ``I_{N_f}`` by default. 
* `mat_b :: Matrix{<:Number}` is a `nf`×`nf` matrix specifying ``M^B_{ff'}``. Facultative, the Hermitian conjugate of `mat_a` by default. 
* `m_kept :: Vector{Int64}` is a list of orbitals that range from 1 to `nm`. Facultative, if specified, only terms for which all ``m_i`` are in the list are kept. 
"""
function GetPairIntTerms(nm :: Int64, nf :: Int64, ps_pot :: Vector{<:Number}, mat_a :: Matrix{<:Number}, mat_b :: Matrix{<:Number} = Matrix(mat_a') ; m_kept :: Vector{Int64} = collect(1 : nm))
    no = nm * nf
    int_el = GetIntMatrix(nm, ps_pot)
    tms = Term[]
    for o1 = 1 : no
        m1 = div(o1 - 1, nf) + 1
        f1 = mod(o1 - 1, nf) + 1
        if (m1 ∉ m_kept) continue end
        for o2 = 1 : no 
            m2 = div(o2 - 1, nf) + 1
            f2 = mod(o2 - 1, nf) + 1
            if (m2 ∉ m_kept) continue end
            if (o1 == o2) continue end
            if (abs(mat_a[f1, f2]) < 1E-13) continue end 
            for o3 = 1 : no
                m3 = div(o3 - 1, nf) + 1
                f3 = mod(o3 - 1, nf) + 1
                if (m3 ∉ m_kept) continue end
                m4 = m1 + m2 - m3 
                if (m4 <= 0 || m4 > nm) continue end
                if (m4 ∉ m_kept) continue end
                for f4 = 1 : nf 
                    if (abs(mat_b[f3, f4]) < 1E-13) continue end
                    o4 = (m4 - 1) * nf + f4
                    if (o3 == o4) continue end
                    val = mat_a[f1, f2] * mat_b[f3, f4] * int_el[m1, m2, m3]
                    if (abs(val) < 1E-15) continue end 
                    push!(tms, Term(val, [1, o1, 1, o2, 0, o3, 0, o4]))
                end
            end
        end
    end
    return SimplifyTerms(tms)
end
GetPairIntTerms(nm :: Int64, nf :: Int64, mat_a :: Matrix{<:Number}, mat_b :: Matrix{<:Number} = Matrix(mat_a') ; m_kept :: Vector{Int64} = collect(1 : nm)) = GetPairIntTerms(nm, nf, [1.0], mat_a, mat_b ; m_kept)


"""
    GetPolTerms(nm :: Int64, nf :: Int64[, mat :: Matrix{<:Number}][ ; fld_m :: Vector{<:Number}]) :: Terms

Return the polarisation term in the Hamiltonian 
```math 
∑_{mff'}c^{†}_{mf}M_{ff'}c_{mf'}.
```

# Arguments 

* `nm :: Int64` is the number of orbitals.
* `nf :: Int64` is the number of flavours. 
* `mat :: Matrix{<:Number}` is a `nf`×`nf` matrix specifying ``M_{ff'}``. Facultative, ``I_{N_f}`` by default. 
* `fld_m :: Vector{<:Number}` gives an orbital dependent polarisation
```math 
∑_{mff'}h_mc^{†}_{mf}M_{ff'}c_{mf'}
```
Facultative. 
"""
function GetPolTerms(nm :: Int64, nf :: Int64, mat :: Matrix{<:Number} = Matrix{Float64}(I, nf, nf) ; fld_m :: Vector{<:Number} = fill(1, nm))
    no = nm * nf
    tms = Term[]
    for o1 = 1 : no
        m1 = div(o1 - 1, nf) + 1
        f1 = mod(o1 - 1, nf) + 1
        for f2 = 1 : nf 
            if abs(mat[f1, f2]) < 1E-13 continue end
            o2 = (m1 - 1) * nf + f2
            push!(tms, Term(mat[f1, f2] * fld_m[m1], [1, o1, 0, o2]))
        end
    end
    return SimplifyTerms(tms)
end


"""
    GetL2Terms(nm :: Int64, nf :: Int64) :: Terms

Return the terms for the total angular momentum.

# Arguments
* `nm :: Int64` is the number of orbitals.
* `nf :: Int64` is the number of flavours.
"""
function GetL2Terms(nm :: Int64, nf :: Int64)
    s = (nm - 1) / 2.0
    no = nm * nf
    tms_lz = 
        [ begin m = div(o - 1, nf)
            Term(m - s, [1, o, 0, o])
        end for o = 1 : no ]
    tms_lp = 
        [ begin m = div(o - 1, nf)
            Term(sqrt(m * (nm - m)), [1, o, 0, o - nf])
        end for o = nf + 1 : no ]
    tms_lm = tms_lp' 
    tms_l2 = tms_lz * tms_lz - tms_lz + tms_lp * tms_lm
    return SimplifyTerms(tms_l2)
end


"""
    GetC2Terms(nm :: Int64, nf :: Int64, mat_gen :: Vector{Matrix{<:Number}}[, mat_tr :: Vector{Matrix{<:Number}}]) :: Terms

Return the terms for the quadratic Casimir of the flavour symmetry.
```math
    C_2=∑_{imm'}\\frac{(c^†_{mf_1}G_{i,f_1f_2}c_{mf_2})(c^†_{m'f_3}G^†_{i,f_3f_4}c_{m'f_4})}{\\operatorname{tr}G_i^†G_i}-∑_{imm'}\\frac{(c^†_{mf_1}T_{i,f_1f_2}c_{mf_2})(c^†_{m'f_3}T^†_{i,f_3f_4}c_{m'f_4})}{\\operatorname{tr}T_i^†T_i}
```
where ``G_i`` are the generator matrices, and ``T_i`` are the trace matrices. 

# Arguments
* `nm :: Int64` is the number of orbitals.
* `nf :: Int64` is the number of flavours.
* `mat_gen :: Vector{Matrix{Number}})` is a list of the matrices that gives the generators. It will automatically be normalised such that its square traces to \$1/2\$. 
* `mat_tr :: Vector{Matrix{Number}})` is a list of trace matrices that will be normalised automatically and substracted. Facultative.
"""
function GetC2Terms(nm :: Int64, nf :: Int64, mat_gen :: Vector{<:AbstractMatrix{<:Number}}, mat_tr :: Vector{<:AbstractMatrix{<:Number}} = Matrix{Float64}[])
    return SimplifyTerms(
        sum([GetPolTerms(nm, nf, Matrix(mati')) * GetPolTerms(nm, nf, mati) / tr(mati' * mati) * 0.5 for mati in mat_gen])
         - (isempty(mat_tr) ? Term[] : sum([GetPolTerms(nm, nf, Matrix(mati')) * GetPolTerms(nm, nf, mati) / tr(mati' * mati) * 0.5 for mati in mat_tr])))
end 