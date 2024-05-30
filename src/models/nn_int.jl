"""
    function GetSnBasis(cfs :: Confs, nf :: Int64 ; qn_r :: Int64 = 0, perm :: Vector, qn_z :: Vector{<:Number}) :: Basis

Return the basis where the ``\\pi``-rotation along ``y``-axis ``\\mathscr{R}`` and certain permutationss of flavour are implemented. Quantum numbers set to zero signify that they are not conserved. 

# Arguments

- `cfs :: Confs` is the configurations generated by [`GetLzConfs`](@ref) or [`GetLzZnConfs`](@ref).
- `nf :: Int64` is the number of flavours
- `qn_r :: Int64` is the quantum number for  ``\\pi`` rotation along ``y``-axis compared with the ground state. Facultive, 0 by default.
- `perm :: Vector{Vector{Int64}}` is a list where each element specifies a permutation of flavour indices (from 1 to ``N_f``) in the format of a cycle. Facultive, empty by default.
- `qn_z :: Vector{<:Number}` is a list where each element specifies the quantum number under the flavour permutation. Facultive, empty by default. 
"""
function GetSnBasis(cfs :: Confs, nf :: Int64 ; qn_r :: Int64 = 0, perm :: Vector = [], qn_z :: Vector{<:Number} = Number[]) 
    no = cfs.no
    nm = no ÷ nf
    qn_r1 = qn_r
    if (nm % 4 >= 2) qn_r1 = -qn_r end
    cyc = Vector{Int64}(undef, 0)
    qnz_s = Vector{ComplexF64}(undef, 0)
    perm_o = []
    ph_o = []
    fac_o = []
    if qn_r != 0
        push!(perm_o, [begin 
                f1 = mod(o - 1, nf) ; m1 = div(o - 1, nf)
                1 + f1 + nf * (nm - 1 - m1)
            end for o = 1 : no])
        push!(ph_o, fill(0, no)) 
        push!(fac_o, fill(ComplexF64(1), no)) 
        push!(qnz_s, qn_r1)
        push!(cyc, 2)
    end
    for i in eachindex(perm)
        if qn_z[i] == 0 continue end
        f_perm = collect(1 : nf)
        for j in 1 : length(perm[i]) - 1
            f_perm[perm[i][j]] = perm[i][j + 1] 
        end
        f_perm[perm[i][end]] = perm[i][1]
        push!(perm_o, [begin 
                f1 = (o - 1) % nf + 1 ; m1 = (o - 1) ÷ nf
                f_perm[f1] + nf * m1
            end for o = 1 : no])
        push!(ph_o, fill(0, no))
        push!(fac_o, fill(ComplexF64(1), no)) 
        push!(qnz_s, qn_z[i]) 
        push!(cyc, length(perm[i]))
    end
    @show qnz_s 
    @show cyc 
    @show perm_o 
    # Generate the basis and print the dimension
    return Basis(cfs, qnz_s, cyc, perm_o, ph_o, fac_o)
end


"""
    GetDenIntTerms(nm :: Int64, nf :: Int64 ; ps_pot :: Vector{<:Number}, mat_a :: Matrix{<:Number}, mat_b :: Matrix{<:Number}) :: Vector{Term}

Return the normal-ordered density-density term in the Hamiltonian 
```math 
\\sum_{\\{m_i,f_i\\}}U_{m_1m_2m_3m_4}M^A_{f_1f_4}M^B_{f_2f_3}c^{\\dagger}_{m_1f_1}c^{\\dagger}_{m_2f_2}c_{m_3f_3}c_{m_4f_4}
```

# Arguments 

- `nm :: Int64` is the number of orbitals.
- `nf :: Int64` is the number of flavours.
- `ps_pot :: Vector{<:Number}` is a list of numbers specifying the pseudopotentials for the interacting matrix ``U_{m_1m_2m_3m_4}``. Facultive, `[1.0]` by default.
- `mat_a :: Matrix{<:Number}` is a `nf`\\*`nf` matrix specifying ``M^A_{ff'}``. Facultive, ``I_{N_f}`` by default. 
- `mat_b :: Matrix{<:Number}` is a `nf`\\*`nf` matrix specifying ``M^B_{ff'}``. Facultive, the Hermitian conjugate of `mat_a` by default. 

"""
function GetDenIntTerms(nm :: Int64, nf :: Int64 ; ps_pot :: Vector{<:Number} = [1.0], mat_a :: Matrix{<:Number} = Matrix{Float64}(I, nf, nf), mat_b :: Matrix{<:Number} = Matrix(mat_a'))
    no = nm * nf
    int_el = GetIntMatrix(nm, ps_pot)
    tms = Vector{Term}(undef, 0)
    # Go through all the m1-up, m2-down, m3-down, m4-up and m4 = m1 + m2 - m3
    for o1 = 1 : no
        m1 = div(o1 - 1, nf) + 1
        f1 = mod(o1 - 1, nf) + 1
        for o2 = 1 : no 
            m2 = div(o2 - 1, nf) + 1
            f2 = mod(o2 - 1, nf) + 1
            if (o1 == o2) continue end
            # if (f1 < f2) continue end # f1 >= f2
            # if (f1 == f2 && m1 <= m2) continue end 
            for o3 = 1 : no
                m3 = div(o3 - 1, nf) + 1
                f3 = mod(o3 - 1, nf) + 1
                if (abs(mat_b[f2, f3]) < 1E-13) continue end 
                m4 = m1 + m2 - m3 
                if (m4 <= 0 || m4 > nm) continue end
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


"""
    GetPairIntTerms(nm :: Int64, nf :: Int64 ; ps_pot :: Vector{<:Number}, mat_a :: Matrix{<:Number}, mat_b :: Matrix{<:Number}) :: Vector{Term}

Return the normal-ordered pair-pair interaction term in the Hamiltonian 
```math 
\\sum_{\\{m_i,f_i\\}}U_{m_1m_2m_3m_4}M^A_{f_1f_2}M^B_{f_4f_3}c^{\\dagger}_{m_1f_1}c^{\\dagger}_{m_2f_2}c_{m_3f_3}c_{m_4f_4}
```

# Arguments 

- `nm :: Int64` is the number of orbitals.
- `nf :: Int64` is the number of flavours.
- `ps_pot :: Vector{<:Number}` is a list of numbers specifying the pseudopotentials for the interacting matrix ``U_{m_1m_2m_3m_4}``. Facultive, `[1.0]` by default.
- `mat_a :: Matrix{<:Number}` is a `nf`\\*`nf` matrix specifying ``M^A_{ff'}``. Facultive, ``I_{N_f}`` by default. 
- `mat_b :: Matrix{<:Number}` is a `nf`\\*`nf` matrix specifying ``M^B_{ff'}``. Facultive, the Hermitian conjugate of `mat_a` by default. 

"""
function GetPairIntTerms(nm :: Int64, nf :: Int64 ; ps_pot :: Vector{<:Number} = [1.0], mat_a :: Matrix{<:Number}, mat_b :: Matrix{<:Number} = Matrix(mat_a'))
    no = nm * nf
    int_el = GetIntMatrix(nm, ps_pot)
    tms = Vector{Term}(undef, 0)
    # Go through all the m1-up, m2-down, m3-down, m4-up and m4 = m1 + m2 - m3
    for o1 = 1 : no
        m1 = div(o1 - 1, nf) + 1
        f1 = mod(o1 - 1, nf) + 1
        for o2 = 1 : no 
            m2 = div(o2 - 1, nf) + 1
            f2 = mod(o2 - 1, nf) + 1
            if (o1 == o2) continue end
            if (abs(mat_a[f1, f2]) < 1E-13) continue end 
            # if (f1 < f2) continue end # f1 >= f2
            # if (f1 == f2 && m1 <= m2) continue end 
            for o3 = 1 : no
                m3 = div(o3 - 1, nf) + 1
                f3 = mod(o3 - 1, nf) + 1
                m4 = m1 + m2 - m3 
                if (m4 <= 0 || m4 > nm) continue end
                for f4 = 1 : nf 
                    if (abs(mat_b[f4, f3]) < 1E-13) continue end
                    o4 = (m4 - 1) * nf + f4
                    if (o3 == o4) continue end
                    val = mat_a[f1, f2] * mat_b[f4, f3] * int_el[m1, m2, m3]
                    if (abs(val) < 1E-15) continue end 
                    push!(tms, Term(val, [1, o1, 1, o2, 0, o3, 0, o4]))
                end
            end
        end
    end
    return SimplifyTerms(tms)
end

"""
    GetPolTerms(nm :: Int64, nf :: Int64 ; mat :: Matrix{<:Number} = Matrix{Float64}(I, nf, nf)) :: Vector{Term}

Return the polarisation term in the Hamiltonian 
```math 
\\sum_{mff'}c^{\\dagger}_{mf}M_{ff'}c_{mf'}
```

# Arguments 

- `nm :: Int64` is the number of orbitals ;
- `nf :: Int64` is the number of flavours ; 
- `mat :: Matrix{<:Number}` is a `nf`\\*`nf` matrix specifying ``M_{ff'}``. Facultive, ``I_{N_f}`` by default. 

"""
function GetPolTerms(nm :: Int64, nf :: Int64 ; mat :: Matrix{<:Number} = Matrix{Float64}(I, nf, nf))
    no = nm * nf
    tms = Vector{Term}(undef, 0)
    for o1 = 1 : no
        m1 = div(o1 - 1, nf) + 1
        f1 = mod(o1 - 1, nf) + 1
        for f2 = 1 : nf 
            if abs(mat[f1, f2]) < 1E-13 continue end
            o2 = (m1 - 1) * nf + f2
            push!(tms, Term(mat[f1, f2], [1, o1, 0, o2]))
        end
    end
    return SimplifyTerms(tms)
end