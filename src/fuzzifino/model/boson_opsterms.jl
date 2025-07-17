export GetBosonDenIntSTerms, GetBosonPairIntSTerms, GetBosonPolSTerms, GetBosonL2STerms, GetBosonC2STerms, GetL2STerms


"""
    GetBosonDenIntSTerms(nm :: Int64, nf :: Int64[, ps_pot :: Vector{<:Number}][, mat_a :: Matrix{<:Number}[, mat_b :: Matrix{<:Number}]]) :: Terms

Return the normal-ordered density-density term in the Hamiltonian for bosons
```math 
∑_{\\{m_i,f_i\\}}U_{m_1m_2m_3m_4}M^A_{f_1f_4}M^B_{f_2f_3}b^{†}_{m_1f_1}b^{†}_{m_2f_2}b_{m_3f_3}b_{m_4f_4}.
```

# Arguments 

* `nm :: Int64` is the number of orbitals.
* `nf :: Int64` is the number of flavours.
* `ps_pot :: Vector{<:Number}` is a list of numbers specifying the pseudopotentials for the interacting matrix ``U_{m_1m_2m_3m_4}``. Facultative, `[1.0]` by default. 
* `mat_a :: Matrix{<:Number}` is a `nf`×`nf` matrix specifying ``M^A_{ff'}``. Facultative, ``I_{N_f}`` by default. 
* `mat_b :: Matrix{<:Number}` is a `nf`×`nf` matrix specifying ``M^B_{ff'}``. Facultative, the Hermitian conjugate of `mat_a` by default. 
"""
function GetBosonDenIntSTerms(nm :: Int64, nf :: Int64, ps_pot :: Vector{<:Number}, mat_a :: Matrix{<:Number} = Matrix{Float64}(I, nf, nf), mat_b :: Matrix{<:Number} = Matrix(mat_a'))
    no = nm * nf
    int_el = GetIntMatrix(nm, ps_pot)
    tms = STerm[]
    for o1 = 1 : no
        m1 = div(o1 - 1, nf) + 1
        f1 = mod(o1 - 1, nf) + 1
        for o2 = 1 : no 
            m2 = div(o2 - 1, nf) + 1
            f2 = mod(o2 - 1, nf) + 1
            for o3 = 1 : no
                m3 = div(o3 - 1, nf) + 1
                f3 = mod(o3 - 1, nf) + 1
                if (abs(mat_b[f2, f3]) < 1E-13) continue end 
                m4 = m1 + m2 - m3 
                if (m4 <= 0 || m4 > nm) continue end
                for f4 = 1 : nf 
                    if (abs(mat_a[f1, f4]) < 1E-13) continue end
                    o4 = (m4 - 1) * nf + f4
                    val = mat_a[f1, f4] * mat_b[f2, f3] * int_el[m1, m2, m3]
                    if (abs(val) < 1E-15) continue end 
                    push!(tms, STerm(val, [1, -o1, 1, -o2, 0, -o3, 0, -o4]))
                end
            end
        end
    end
    return SimplifyTerms(tms)
end
GetBosonDenIntSTerms(nm :: Int64, nf :: Int64, mat_a :: Matrix{<:Number} = Matrix{Float64}(I, nf, nf), mat_b :: Matrix{<:Number} = Matrix(mat_a')) = GetBosonDenIntSTerms(nm, nf, [1.0], mat_a, mat_b)


"""
    GetBosonDenIntSTerms(nm :: Int64, nf :: Int64, ps_pot :: Vector{<:Number}, mat_a :: Vector{<:AbstractMatrix{<:Number}}[, mat_b :: Vector{<:AbstractMatrix{<:Number}}]) :: Terms

Return the sum of a series of normal-ordered density-density term in the Hamiltonian for bosons
```math 
∑_{\\{m_i,f_i,α\\}}U_{m_1m_2m_3m_4}(M^A_{α})_{f_1f_4}(M^B_{α})_{f_2f_3}b^{†}_{m_1f_1}b^{†}_{m_2f_2}b_{m_3f_3}b_{m_4f_4}.
```

# Arguments 

* `nm :: Int64` is the number of orbitals.
* `nf :: Int64` is the number of flavours.
* `ps_pot :: Vector{<:Number}` is a list of numbers specifying the pseudopotentials for the interacting matrix ``U_{m_1m_2m_3m_4}``. Facultative, `[1.0]` by default.
* `mat_a :: Vector{<:AbstractMatrix{<:Number}}` is a vector of `nf`×`nf` matrix specifying ``(M^A_{α})_{ff'}``. Facultative, ``I_{N_f}`` by default. 
* `mat_b :: Vector{<:AbstractMatrix{<:Number}}` is a vector of `nf`×`nf` matrix specifying ``(M^B_{α})_{ff'}``. Facultative, the Hermitian conjugate of `mat_a` by default. 
"""
function GetBosonDenIntSTerms(nm :: Int64, nf :: Int64, ps_pot :: Vector{<:Number}, mats_a :: Vector{<:AbstractMatrix{<:Number}}, mats_b :: Vector{<:AbstractMatrix{<:Number}} = [Matrix(mat_a') for mat_a in mats_a])
    return sum([GetBosonDenIntSTerms(nm, nf, ps_pot, mats_a[i], mats_b[i]) for i in eachindex(mats_a)])
end
GetBosonDenIntSTerms(nm :: Int64, nf :: Int64, mats_a :: Vector{<:AbstractMatrix{<:Number}}, mats_b :: Vector{<:AbstractMatrix{<:Number}} = [Matrix(mat_a') for mat_a in mats_a]) = GetBosonDenIntSTerms(nm, nf, [1.0], mats_a, mats_b)


"""
    GetBosonPairIntSTerms(nm :: Int64, nf :: Int64, ps_pot :: Vector{<:Number}, mat_a :: Matrix{<:Number}[, mat_b :: Matrix{<:Number}]) :: Terms

Return the normal-ordered pair-pair interaction term in the Hamiltonian for bosons
```math 
∑_{\\{m_i,f_i\\}}U_{m_1m_2m_3m_4}M^A_{f_1f_2}M^B_{f_3f_4}b^{†}_{m_1f_1}b^{†}_{m_2f_2}b_{m_3f_3}b_{m_4f_4}.
```

# Arguments 

* `nm :: Int64` is the number of orbitals.
* `nf :: Int64` is the number of flavours.
* `ps_pot :: Vector{<:Number}` is a list of numbers specifying the pseudopotentials for the interacting matrix ``U_{m_1m_2m_3m_4}``. 
* `mat_a :: Matrix{<:Number}` is a `nf`×`nf` matrix specifying ``M^A_{ff'}``. Facultative, ``I_{N_f}`` by default. 
* `mat_b :: Matrix{<:Number}` is a `nf`×`nf` matrix specifying ``M^B_{ff'}``. Facultative, the Hermitian conjugate of `mat_a` by default. 
"""
function GetBosonPairIntSTerms(nm :: Int64, nf :: Int64, ps_pot :: Vector{<:Number}, mat_a :: Matrix{<:Number}, mat_b :: Matrix{<:Number} = Matrix(mat_a'))
    no = nm * nf
    int_el = GetIntMatrix(nm, ps_pot)
    tms = STerm[]
    for o1 = 1 : no
        m1 = div(o1 - 1, nf) + 1
        f1 = mod(o1 - 1, nf) + 1
        for o2 = 1 : no 
            m2 = div(o2 - 1, nf) + 1
            f2 = mod(o2 - 1, nf) + 1
            if (abs(mat_a[f1, f2]) < 1E-13) continue end
            for o3 = 1 : no
                m3 = div(o3 - 1, nf) + 1
                f3 = mod(o3 - 1, nf) + 1
                m4 = m1 + m2 - m3 
                if (m4 <= 0 || m4 > nm) continue end
                for f4 = 1 : nf 
                    if (abs(mat_b[f3, f4]) < 1E-13) continue end
                    o4 = (m4 - 1) * nf + f4
                    val = mat_a[f1, f2] * mat_b[f3, f4] * int_el[m1, m2, m3]
                    if (abs(val) < 1E-15) continue end 
                    push!(tms, STerm(val, [1, -o1, 1, -o2, 0, -o3, 0, -o4]))
                end
            end
        end
    end
    return SimplifyTerms(tms)
end
GetBosonPairIntSTerms(nm :: Int64, nf :: Int64, mat_a :: Matrix{<:Number}, mat_b :: Matrix{<:Number} = Matrix(mat_a')) = GetBosonPairIntSTerms(nm, nf, [1.0], mat_a, mat_b)


"""
    GetBosonPolSTerms(nm :: Int64, nf :: Int64[, mat :: Matrix{<:Number}][ ; fld_m :: Vector{<:Number}]) :: STerms

Return the polarisation term in the Hamiltonian for bosons
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
function GetBosonPolSTerms(nm :: Int64, nf :: Int64, mat :: Matrix{<:Number} = Matrix{Float64}(I, nf, nf) ; fld_m :: Vector{<:Number} = fill(1, nm))
    no = nm * nf
    tms = STerm[]
    for o1 = 1 : no
        m1 = div(o1 - 1, nf) + 1
        f1 = mod(o1 - 1, nf) + 1
        for f2 = 1 : nf 
            if abs(mat[f1, f2]) < 1E-13 continue end
            o2 = (m1 - 1) * nf + f2
            push!(tms, STerm(mat[f1, f2] * fld_m[m1], [1, -o1, 0, -o2]))
        end
    end
    return SimplifyTerms(tms)
end

"""
    GetL2STerms(nmf :: Int64, nff :: Int64, nmb :: Int64, nfb :: Int64) :: STerms

Return the terms for the total angular momentum for bosons and .

# Arguments
* `nmf :: Int64` is the number of fermion orbitals.
* `nff :: Int64` is the number of fermion flavours.
* `nmb :: Int64` is the number of boson orbitals.
* `nfb :: Int64` is the number of boson flavours.
"""
function GetL2STerms(nmf :: Int64, nff :: Int64, nmb :: Int64, nfb :: Int64)
    tms_lzf = [ STerm(m - (nmf + 1) / 2, [1, (m - 1) * nff + f, 0, (m - 1) * nff + f]) for m = 1 : nmf for f = 1 : nff ]
    tms_lzb = [ STerm(m - (nmb + 1) / 2, [1, -((m - 1) * nfb + f), 0, -((m - 1) * nfb + f)]) for m = 1 : nmb for f = 1 : nfb ]

    tms_lpf = [ STerm(sqrt((nmf - m) * m), [1, m * nff + f, 0, (m - 1) * nff + f]) for m = 1 : nmf - 1 for f = 1 : nff ]
    tms_lpb = [ STerm(sqrt((nmb - m) * m), [1,-(m * nfb + f), 0,-((m - 1) * nfb + f)]) for m = 1 : nmb - 1 for f = 1 : nfb ]

    tms_lz = tms_lzf + tms_lzb
    tms_lp = tms_lpf + tms_lpb
    tms_lm = tms_lp'
    
    return SimplifyTerms(tms_lz * tms_lz - tms_lz + tms_lp * tms_lm)
end ;


"""
    GetBosonL2STerms(nm :: Int64, nf :: Int64) :: STerms

Return the terms for the total angular momentum for bosons.

# Arguments
* `nm :: Int64` is the number of orbitals.
* `nf :: Int64` is the number of flavours.
"""
GetBosonL2STerms(nm :: Int64, nf :: Int64) = GetL2STerms(0, 0, nm, nf)


"""
    GetBosonC2STerms(nm :: Int64, nf :: Int64, mat_gen :: Vector{Matrix{<:Number}}[, mat_tr :: Vector{Matrix{<:Number}}]) :: STerms

Return the terms for the quadratic Casimir of the flavour symmetry for bosons.
```math
    C_2=∑_{imm'}\\frac{(b^†_{mf_1}G_{i,f_1f_2}b_{mf_2})(b^†_{m'f_3}G^†_{i,f_3f_4}b_{m'f_4})}{2\\operatorname{tr}G_i^†G_i}-∑_{imm'}\\frac{(b^†_{mf_1}T_{i,f_1f_2}b_{mf_2})(b^†_{m'f_3}T^†_{i,f_3f_4}b_{m'f_4})}{2\\operatorname{tr}T_i^†T_i}
```
where ``G_i`` are the generator matrices, and ``T_i`` are the trace matrices. 

# Arguments
* `nm :: Int64` is the number of orbitals.
* `nf :: Int64` is the number of flavours.
* `mat_gen :: Vector{Matrix{Number}}` is a list of the matrices that gives the generators. It will automatically be normalised such that its square traces to \$1/2\$. 
* `mat_tr :: Vector{Matrix{Number}}` is a list of trace matrices that will be normalised automatically and substracted. Facultative.
"""
function GetBosonC2STerms(nm :: Int64, nf :: Int64, mat_gen :: Vector{<:AbstractMatrix{<:Number}}, mat_tr :: Vector{<:AbstractMatrix{<:Number}} = Matrix{Float64}[])
    return SimplifyTerms(
        sum([GetBosonPolSTerms(nm, nf, Matrix(mati')) * GetBosonPolSTerms(nm, nf, mati) / tr(mati' * mati) * 0.5 for mati in mat_gen])
         - (isempty(mat_tr) ? STerm[] : sum([GetBosonPolSTerms(nm, nf, Matrix(mati')) * GetBosonPolSTerms(nm, nf, mati) / tr(mati' * mati) * 0.5 for mati in mat_tr])))
end 
