export GetTorusLz2QNDiag, GetTorusIntMatrix, GetTorusTranlQNOffd, GetTorusDenIntTerms, GetTorusPairIntTerms


""" 
    GetTorusLz2QNDiag(nm :: Int64, nf :: Int64) :: QNDiag 

Return the QNDiag of twice the angular momentum ``2L_z\\mod 2N_m``, with an offset such that one fully filled LLL has angular momentum 0, implemented as 
```julia
QNDiag("Lz", [2 * m - 1 + mod(nm, 2) for m = 1 : nm for f = 1 : nf], 2 * nm)
```
"""
GetTorusLz2QNDiag(nm :: Int64, nf :: Int64) = QNDiag("Lz", [2 * m - 1 + mod(nm, 2) for m = 1 : nm for f = 1 : nf], 2 * nm)


""" 
    GetTorusTranlQNOffd(nm :: Int64, nf :: Int64, Δ :: Int64 = 1) :: QNOffd

Return the QNOffd of translation 
```math 
    c_m↦c_{m+Δ\\mod N_m}
```
on the torus, implemented as 
```julia
    QNOffd([((m - 1 + Δ) % nm) * nf + f for m = 1 : nm for f = 1 : nf], nm)
```
"""
GetTorusTranlQNOffd(nm :: Int64, nf :: Int64, Δ :: Int64 = 1) = QNOffd([((m - 1 + Δ) % nm) * nf + f for m = 1 : nm for f = 1 : nf], nm)

function LaguerrePsPot(ps_pot :: Vector{<:Number}, q2 :: Number)
    lag = zeros(Float64, length(ps_pot))
    lag[1] = 1.0
    (length(ps_pot) == 1) && return sum(lag .* ps_pot)
    lag[2] = 1.0 - q2
    for k = 1 : length(ps_pot) - 2 
        lag[k + 2] = ((2 * k + 1 - q2) * lag[k + 1] - k * lag[k]) / (k + 1)
    end
    return sum(lag .* ps_pot)
end


"""
    GetTorusIntMatrix(nm :: Int64, lx :: Number, ps_pot :: Vector{<:Number}) :: Array{ComplexF64, 3}
    GetTorusIntMatrix(nm :: Int64, ps_pot :: Vector{<:Number} ; aspect_ratio :: Number = 1.0) :: Array{ComplexF64, 3}

Gives the interaction matrix ``U_{m_1m_2m_3m_4}`` from the pseudopotentials.
```math
    U_{m_1m_2m_3m_4}=δ'_{m_1+m_2,m_3+m_4}\\frac{1}{N_m}∑_{l,𝐪}δ'_{m_1-m_4,t}U_lL_l(q^2)e^{-q^2/2}e^{2πis(m_1-m_3)/N_m}
```
where ``𝐪=(2πs/L_x,2πt/L_y), s,t∈ℤ``, ``L_xL_y=2πN_m`` and the Kronecker ``δ`` is defined in a sense of mod ``N_m``

# Argument

* `nm :: Int64` is the number of orbitals.
* `lx :: Number` is the length along ``x``-direction. Facultative.
* `ps_pot :: Vector{<:Number}` is the vector of non-zero pseudopotentials.
* `aspect_ratio :: Number` is the ratio ``L_y/L_x``. Facultative, at most one of `lx` and `aspect_ratio` is given. If both are omitted, ``L_x=L_y`` is taken. 

# Output
* A `nm`×`nm`×`nm` array giving the interaction matrix ``U_{m_1m_2m_3m_4}`` where ``m_4=m_1+m_2-m_3\\mod N_m``.
"""
function GetTorusIntMatrix(nm :: Int64, lx :: Number, ps_pot :: Vector{<:Number})
    int_el = zeros(ComplexF64, nm, nm, nm)
    b = 2 * π * nm / lx
    sm = ceil(lx / (2 * π) * √(-2 * log(eps(Float64))))
    km = ceil(b / (2 * π * nm) * √(-2 * log(eps(Float64))))
    for m1 in 1 : nm, m2 in 1 : nm, m3 = 1 : nm
        m4 = mod1(m1 + m2 - m3, nm)
        val = ComplexF64(0) 
        for k = -km : km, s = -sm : sm
            qx = 2 * π / lx * s 
            qy = 2 * π / b * (m1 - m4 + k * nm)
            q2 = qx ^ 2 + qy ^ 2
            val += exp(-q2 / 2) * LaguerrePsPot(ps_pot, q2) * exp(2im * π * s * (m1 - m3) / nm)
        end
        int_el[m1, m2, m3] = val / nm
    end
    return int_el
end
GetTorusIntMatrix(nm :: Int64, ps_pot :: Vector{<:Number} ; aspect_ratio :: Number = 1.0) = GetTorusIntMatrix(nm, √(2 * π * nm / aspect_ratio), ps_pot)


"""
    GetTorusDenIntTerms(nm :: Int64, nf :: Int64, lx :: Number, ps_pot :: Vector{<:Number}[, mat_a :: Matrix{<:Number}[, mat_b :: Matrix{<:Number}]]) :: Terms
    GetTorusDenIntTerms(nm :: Int64, nf :: Int64, ps_pot :: Vector{<:Number}[, mat_a :: Matrix{<:Number}[, mat_b :: Matrix{<:Number}]] ; aspect_ratio :: Number = 1.0) :: Terms

Return the normal-ordered density-density term in the Hamiltonian 
```math 
∑_{\\{m_i,f_i\\}}U_{m_1m_2m_3m_4}M^A_{f_1f_4}M^B_{f_2f_3}c^{†}_{m_1f_1}c^{†}_{m_2f_2}c_{m_3f_3}c_{m_4f_4}.
```

# Arguments 

* `nm :: Int64` is the number of orbitals.
* `nf :: Int64` is the number of flavours.
* `lx :: Number` is the length along ``x``-direction. Facultative.
* `ps_pot :: Vector{<:Number}` is a list of numbers specifying the pseudopotentials for the interacting matrix ``U_{m_1m_2m_3m_4}``. 
* `mat_a :: Matrix{<:Number}` is a `nf`×`nf` matrix specifying ``M^A_{ff'}``. Facultative, ``I_{N_f}`` by default. 
* `mat_b :: Matrix{<:Number}` is a `nf`×`nf` matrix specifying ``M^B_{ff'}``. Facultative, the Hermitian conjugate of `mat_a` by default. 
* `aspect_ratio :: Number` is the ratio ``L_y/L_x``. Facultative, at most one of `lx` and `aspect_ratio` is given. If both are omitted, ``L_x=L_y`` is taken. 
"""
function GetTorusDenIntTerms(nm :: Int64, nf :: Int64, lx :: Number, ps_pot :: Vector{<:Number}, mat_a :: Matrix{<:Number} = Matrix{Float64}(I, nf, nf), mat_b :: Matrix{<:Number} = Matrix(mat_a'))
    no = nm * nf 
    int_el = GetTorusIntMatrix(nm, lx, ps_pot)
    tms = Term[]
    o = zeros(Int64, 4)
    for o[1] = 1 : no, o[2] = 1 : no, o[3] = 1 : no
        (o[1] == o[2]) && continue
        m = cld.(o, nf)
        f = mod1.(o, nf)
        (abs(mat_b[f[2], f[3]]) < 1E-13) && continue 
        m[4] = mod1(m[1] + m[2] - m[3], nm)
        for f[4] = 1 : nf 
            (abs(mat_a[f[1], f[4]]) < 1E-13) && continue 
            o[4] = (m[4] - 1) * nf + f[4]
            (o[3] == o[4]) && continue
            val = mat_a[f[1], f[4]] * mat_b[f[2], f[3]] * int_el[m[1], m[2], m[3]]
            (abs(val) < 1E-15) && continue
            push!(tms, Term(val, [1, o[1], 1, o[2], 0, o[3], 0, o[4]]))
        end
    end
    return SimplifyTerms(tms)
end
GetTorusDenIntTerms(nm :: Int64, nf :: Int64, ps_pot :: Vector{<:Number}, mat_a :: Matrix{<:Number} = Matrix{Float64}(I, nf, nf), mat_b :: Matrix{<:Number} = Matrix(mat_a') ; aspect_ratio :: Number = 1.0) = GetTorusDenIntTerms(nm, nf, √(2 * π * nm / aspect_ratio), ps_pot, mat_a, mat_b)


"""
    GetTorusPairIntTerms(nm :: Int64, nf :: Int64, lx :: Number, ps_pot :: Vector{<:Number}, mat_a :: Matrix{<:Number}[, mat_b :: Matrix{<:Number}]) :: Terms
    GetTorusPairIntTerms(nm :: Int64, nf :: Int64, ps_pot :: Vector{<:Number}, mat_a :: Matrix{<:Number}[, mat_b :: Matrix{<:Number}] ; aspect_ratio :: Number = 1.0) :: Terms

Return the normal-ordered pair-pair interaction term in the Hamiltonian 
```math 
∑_{\\{m_i,f_i\\}}U_{m_1m_2m_3m_4}M^A_{f_1f_2}M^B_{f_3f_4}c^{†}_{m_1f_1}c^{†}_{m_2f_2}c_{m_3f_3}c_{m_4f_4}.
```

# Arguments 

* `nm :: Int64` is the number of orbitals.
* `nf :: Int64` is the number of flavours.
* `lx :: Number` is the length along ``x``-direction. Facultative.
* `ps_pot :: Vector{<:Number}` is a list of numbers specifying the pseudopotentials for the interacting matrix ``U_{m_1m_2m_3m_4}``. 
* `mat_a :: Matrix{<:Number}` is a `nf`×`nf` matrix specifying ``M^A_{ff'}``. Facultative, ``I_{N_f}`` by default. 
* `mat_b :: Matrix{<:Number}` is a `nf`×`nf` matrix specifying ``M^B_{ff'}``. Facultative, the Hermitian conjugate of `mat_a` by default. 
* `m_kept :: Vector{Int64}` is a list of orbitals that range from 1 to `nm`. Facultative, if specified, only terms for which all ``m_i`` are in the list are kept. 
* `aspect_ratio :: Number` is the ratio ``L_y/L_x``. Facultative, at most one of `lx` and `aspect_ratio` is given. If both are omitted, ``L_x=L_y`` is taken. 
"""
function GetTorusPairIntTerms(nm :: Int64, nf :: Int64, lx :: Number, ps_pot :: Vector{<:Number}, mat_a :: Matrix{<:Number}, mat_b :: Matrix{<:Number} = Matrix(mat_a'))
    no = nm * nf 
    int_el = GetTorusIntMatrix(nm, lx, ps_pot)
    tms = Term[]
    o = zeros(Int64, 4)
    for o[1] = 1 : no, o[2] = 1 : no, o[3] = 1 : no
        (o[1] == o[2]) && continue
        m = cld.(o, nf)
        f = mod1.(o, nf)
        (abs(mat_a[f[1], f[2]]) < 1E-13) && continue 
        m[4] = mod1(m[1] + m[2] - m[3], nm)
        for f[4] = 1 : nf 
            (abs(mat_b[f[3], f[4]]) < 1E-13) && continue 
            o[4] = (m[4] - 1) * nf + f[4]
            (o[3] == o[4]) && continue
            val = mat_a[f[1], f[2]] * mat_b[f[3], f[4]] * int_el[m[1], m[2], m[3]]
            (abs(val) < 1E-15) && continue
            push!(tms, Term(val, [1, o[1], 1, o[2], 0, o[3], 0, o[4]]))
        end
    end
    return SimplifyTerms(tms)
end
GetTorusPairIntTerms(nm :: Int64, nf :: Int64, ps_pot :: Vector{<:Number}, mat_a :: Matrix{<:Number}, mat_b :: Matrix{<:Number} = Matrix(mat_a') ; aspect_ratio :: Number = 1.0) = GetTorusPairIntTerms(nm, nf, √(2 * π * nm / aspect_ratio), ps_pot, mat_a, mat_b)