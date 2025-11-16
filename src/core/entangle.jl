export StateDecompMat, GetEntSpec


"""
    StateDecompMat(st :: Vector{<:Number}, bs0 :: Basis, bsa :: Basis, bsb :: Basis, amp_oa :: Vector{ComplexF64}, amp_ob :: Vector{ComplexF64}) :: Matrix{ComplexF64}

Decompose a state ``|ψ⟩=v_I|I⟩`` into a direct-product basis of two subsystems ``|ψ⟩=M_{JI}|I_A⟩|J_B⟩``

# Arguments 

* `st :: Vector{<:Number}` is the state to be decomposed into direct-product basis of two subsystems.
* `bs0 :: Basis` is the total basis. 
* `bsa :: Basis` is the basis for the subsystem A.
* `bsb :: Basis` is the basis for the subsystem B. 
* `amp_oa :: Vector{ComplexF64}` is a complex list of length `no` that specifies the amplitute of each orbital in the subsystem A. For a non-local basis, we decompose each electron into creation operators in two subsystems ``c^†_o=a_{o,A}c^†_{o,A}+a_{o,B}c^†_{o,B}`` and this list specifies ``a_{o,A}``. This is equivalent to ``√{ℱ_{m,A}}`` in [PRB 85, 125308 (2012)](https://dx.doi.org/10.1103/PhysRevB.85.125308) with an extra phase factor. 
* `amp_ob :: Vector{ComplexF64}` is a complex list of length `no` that specifies the amplitute of each orbital in the subsystem B. 

# Output

A complex matrix of dimension `bsb.dim * bsa.dim` that corresponds to the state in the decomposed basis ``|ψ⟩=M_{JI}|I_A⟩|J_B⟩``. This is equivalent to ``R_{μν}^A/√p`` in [PRB 85, 125308 (2012)](https://dx.doi.org/10.1103/PhysRevB.85.125308). After calculating all the sectors, the reduced density matrix will be ``ρ_B=⊕𝐌𝐌^†``.
"""
function StateDecompMat(st :: Vector{<:Number}, bs0 :: Basis, bsa :: Basis, bsb :: Basis, amp_oa :: Vector{<:Number}, amp_ob :: Vector{<:Number})
    st_dcp = Matrix{ComplexF64}(undef, bsb.dim, bsa.dim) ;
    @ccall Libpath.__ent_MOD_decomp_basis(bs0.cfs.no :: Ref{Int64}, bs0.cfs.nor :: Ref{Int64}, 
        bs0.cfs.ncf :: Ref{Int64}, bs0.dim :: Ref{Int64}, bs0.cfs.conf :: Ref{Int64}, bs0.cfs.lid :: Ref{Int64}, bs0.cfs.rid :: Ref{Int64}, bs0.szz :: Ref{Int64}, bs0.cfgr :: Ref{Int64}, bs0.cffac :: Ref{ComplexF64}, bs0.grel :: Ref{Int64}, bs0.grsz :: Ref{Int64}, 
        bsa.cfs.ncf :: Ref{Int64}, bsa.dim :: Ref{Int64}, bsa.cfs.conf :: Ref{Int64}, bsa.cfs.lid :: Ref{Int64}, bsa.cfs.rid :: Ref{Int64}, bsa.szz :: Ref{Int64}, bsa.cfgr :: Ref{Int64}, bsa.cffac :: Ref{ComplexF64}, bsa.grel :: Ref{Int64}, bsa.grsz :: Ref{Int64}, 
        bsb.cfs.ncf :: Ref{Int64}, bsb.dim :: Ref{Int64}, bsb.cfs.conf :: Ref{Int64}, bsb.cfs.lid :: Ref{Int64}, bsb.cfs.rid :: Ref{Int64}, bsb.szz :: Ref{Int64}, bsb.cfgr :: Ref{Int64}, bsb.cffac :: Ref{ComplexF64}, bsb.grel :: Ref{Int64}, bsb.grsz :: Ref{Int64}, 
        ComplexF64.(amp_oa) :: Ref{ComplexF64}, ComplexF64.(amp_ob) :: Ref{ComplexF64}, ComplexF64.(st) :: Ref{ComplexF64}, st_dcp :: Ref{ComplexF64}) :: Nothing
    return st_dcp
end


"""
    GetEntSpec(st :: Vector{<:Number}, bs0 :: Basis, secd_lst :: Vector{Vector{Vector{Int64}}}, secf_lst :: Vector{Vector{Vector{<:Number}}} ; qnd_a :: Vector{QNDiag}, qnd_b :: Vector{QNDiag} = qnd_a, qnf_a :: Vector{QNOffd}, qnf_b :: Vector{QNOffd} = qnf_a, amp_oa :: Vector{<:Number}, amp_ob :: Vector{<:Number} = sqrt.(1 .- abs.(amp_oa .^ 2))) :: Dict{@NamedTuple{secd_a, secf_a, secd_b, secf_b}, Vector{Float64}}

# Arguments 

* `st :: Vector{<:Number}` is the state to be decomposed into direct-product basis of two subsystems.
* `bs0 :: Basis` is the total basis. 
* `secd_lst :: Vector{Vector{Vector{Int64}}}` gives the list of QNDiag sectors of subsystems to be calculated. Each of its elements is a two element vector ; the first specifies the sector for subsystem A, and the second specifies the sector for subsystem B. 
* `secf_lst :: Vector{Vector{Vector{ComplexF64}}}` gives the list of QNOffd sectors of subsystems to be calculated. Each of its elements is a two element vector ; the first specifies the sector for subsystem A, and the second specifies the sector for subsystem B. 
* `qnd_a :: Vector{QNDiag}, qnd_b :: Vector{QNDiag} = qnd_a, qnf_a :: Vector{QNOffd}, qnf_b :: Vector{QNOffd}` specifies the diagonal and off-diagonal quantum numbers of the subsystems A and B. `qnd_b` and `qnf_b` are facultative and the same as `qnd_a` and `qnf_a` by default. 
* `amp_oa :: Vector{ComplexF64}` and `amp_ob :: Vector{ComplexF64}` are complex lists of length `no` that specify the amplitute of each orbital in the subsystems A and B. For a non-local basis, we decompose each electron into creation operators in two subsystems ``c^†_o=a_{o,A}c^†_{o,A}+a_{o,B}c^†_{o,B}`` and this list specifies ``a_{o,A}``. This is equivalent to ``√{ℱ_{m,A}}`` in [PRB 85, 125308 (2012)](https://dx.doi.org/10.1103/PhysRevB.85.125308) with an extra phase factor. 

# Output

A dictionary whose keys are named tuples that specify the sector containing entries `secd_a`, `secf_a`, `secd_b`, `secf_b` and values are lists of eigenvalues of the density matrix in those sectors. 
"""
function GetEntSpec(st :: Vector{<:Number}, bs0 :: Basis, secd_lst :: Vector{Vector{Vector{Int64}}}, secf_lst :: Union{ Vector{Vector{Vector{Int64}}}, Vector{Vector{Vector{Float64}}}, Vector{Vector{Vector{ComplexF64}}} } ; qnd_a :: Vector{QNDiag}, qnd_b :: Vector{QNDiag} = qnd_a, qnf_a :: Vector{QNOffd}, qnf_b :: Vector{QNOffd} = qnf_a, amp_oa :: Vector{<:Number}, amp_ob :: Vector{<:Number} = sqrt.(1 .- abs.(amp_oa .^ 2)), disp_std = false)
    no = bs0.cfs.no 
    nor = bs0.cfs.nor
    dictlock = ReentrantLock()
    EntSpec = Dict{@NamedTuple{secd_a :: Vector{<:Number}, secf_a :: Vector{<:Number}, secd_b :: Vector{<:Number}, secf_b :: Vector{<:Number}}, Vector{Float64}}()
    Threads.@threads for secd in secd_lst
        cfsa = Confs(no, secd[1], qnd_a ; nor, num_th = 1, disp_std)
        if (cfsa.ncf == 0) continue end 
        cfsb = Confs(no, secd[2], qnd_b ; nor, num_th = 1, disp_std)
        if (cfsb.ncf == 0) continue end 
        for secf in secf_lst 
            bsa = Basis(cfsa, secf[1], qnf_a ; num_th = 1, disp_std)
            if (bsa.dim == 0) continue end 
            bsb = Basis(cfsb, secf[2], qnf_b ; num_th = 1, disp_std)
            if (bsb.dim == 0) continue end 
            st_dcp = StateDecompMat(st, bs0, bsa, bsb, amp_oa, amp_ob) 
            ent_spec = abs.(svdvals(st_dcp)) .^ 2
            Threads.lock(dictlock) 
            try
                EntSpec[(secd_a = secd[1], secf_a = secf[1], secd_b = secd[2], secf_b = secf[2])] = ent_spec
            finally
                Threads.unlock(dictlock)
            end
        end
    end
    return EntSpec
end