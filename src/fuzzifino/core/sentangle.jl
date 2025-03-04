import FuzzifiED: StateDecompMat, GetEntSpec


"""
    StateDecompMat(st :: Vector{<:Number}, bs0 :: SBasis, bsa :: SBasis, bsb :: SBasis, amp_ofa :: Vector{<:Number}, amp_oba :: Vector{<:Number}, amp_ofb :: Vector{<:Number}, amp_obb :: Vector{<:Number}) :: Matrix{ComplexF64}

Decompose a state ``|ψ⟩=v_I|I⟩`` into a direct-product basis of two subsystems ``|ψ⟩=M_{JI}|I_A⟩|J_B⟩``

# Arguments 

* `st :: Vector{<:Number}` is the state to be decomposed into direct-product basis of two subsystems.
* `bs0 :: SBasis` is the total basis. 
* `bsa :: SBasis` is the basis for the subsystem A.
* `bsb :: SBasis` is the basis for the subsystem B. 
* `amp_ofa :: Vector{ComplexF64}` is a complex list of length `no` that specifies the fermionic amplitute of each orbital in the subsystem A. For a non-local basis, we decompose each electron into creation operators in two subsystems ``c^†_o=a_{o,A}c^†_{o,A}+a_{o,B}c^†_{o,B}`` and this list specifies ``a_{o,A}``. This is equivalent to ``√{ℱ_{m,A}}`` in [PRB 85, 125308 (2012)](https://dx.doi.org/10.1103/PhysRevB.85.125308) with an extra phase factor. 
* `amp_oba :: Vector{ComplexF64}` is a complex list of length `no` that specifies the bosonic amplitute of each orbital in the subsystem A. 
* `amp_ofb :: Vector{ComplexF64}` is a complex list of length `no` that specifies the fermionic amplitute of each orbital in the subsystem B. 
* `amp_obb :: Vector{ComplexF64}` is a complex list of length `no` that specifies the bosonic amplitute of each orbital in the subsystem B. 

# Output

A complex matrix of dimension `bsb.dim * bsa.dim` that corresponds to the state in the decomposed basis ``|ψ⟩=M_{JI}|I_A⟩|J_B⟩``. This is equivalent to ``R_{μν}^A/√p`` in [PRB 85, 125308 (2012)](https://dx.doi.org/10.1103/PhysRevB.85.125308). After calculating all the sectors, the reduced density matrix will be ``ρ_B=⊕\\mathbf{M}\\mathbf{M}^†``.
"""
function StateDecompMat(st :: Vector{<:Number}, bs0 :: SBasis, bsa :: SBasis, bsb :: SBasis, amp_ofa :: Vector{<:Number}, amp_oba :: Vector{<:Number}, amp_ofb :: Vector{<:Number}, amp_obb :: Vector{<:Number})
    st_dcp = Matrix{ComplexF64}(undef, bsb.dim, bsa.dim) ;
    binom_0 = [ binomial(i + j, i) for i = 0 : bs0.cfs.nebm, j = 0 : bs0.cfs.nob]
    binom_a = [ binomial(i + j, i) for i = 0 : bsa.cfs.nebm, j = 0 : bsa.cfs.nob]
    binom_b = [ binomial(i + j, i) for i = 0 : bsb.cfs.nebm, j = 0 : bsb.cfs.nob]
    binom_c = [ binomial(i + j, i) for i = 0 : bsa.cfs.nebm, j = 0 : bsb.cfs.nebm]
    @ccall Libpathino.__sent_MOD_decomp_sbasis(bs0.cfs.nof :: Ref{Int64}, bs0.cfs.nob :: Ref{Int64}, bs0.cfs.norf :: Ref{Int64}, bs0.cfs.norb :: Ref{Int64}, 
        bs0.cfs.nebm :: Ref{Int64}, bs0.cfs.ncf :: Ref{Int64}, bs0.dim :: Ref{Int64}, bs0.cfs.conff :: Ref{Int64}, bs0.cfs.confb :: Ref{Int64}, bs0.cfs.lid :: Ref{Int64}, bs0.cfs.rid :: Ref{Int64}, bs0.szz :: Ref{Int64}, bs0.cfgr :: Ref{Int64}, bs0.cffac :: Ref{ComplexF64}, bs0.grel :: Ref{Int64}, bs0.grsz :: Ref{Int64}, binom_0 :: Ref{Int64}, 
        bsa.cfs.nebm :: Ref{Int64}, bsa.cfs.ncf :: Ref{Int64}, bsa.dim :: Ref{Int64}, bsa.cfs.conff :: Ref{Int64}, bsa.cfs.confb :: Ref{Int64}, bsa.cfs.lid :: Ref{Int64}, bsa.cfs.rid :: Ref{Int64}, bsa.szz :: Ref{Int64}, bsa.cfgr :: Ref{Int64}, bsa.cffac :: Ref{ComplexF64}, bsa.grel :: Ref{Int64}, bsa.grsz :: Ref{Int64}, binom_a :: Ref{Int64}, 
        bsb.cfs.nebm :: Ref{Int64}, bsb.cfs.ncf :: Ref{Int64}, bsb.dim :: Ref{Int64}, bsb.cfs.conff :: Ref{Int64}, bsb.cfs.confb :: Ref{Int64}, bsb.cfs.lid :: Ref{Int64}, bsb.cfs.rid :: Ref{Int64}, bsb.szz :: Ref{Int64}, bsb.cfgr :: Ref{Int64}, bsb.cffac :: Ref{ComplexF64}, bsb.grel :: Ref{Int64}, bsb.grsz :: Ref{Int64}, binom_b :: Ref{Int64}, 
        ComplexF64.(amp_ofa) :: Ref{ComplexF64}, ComplexF64.(amp_oba) :: Ref{ComplexF64}, ComplexF64.(amp_ofb) :: Ref{ComplexF64}, ComplexF64.(amp_obb) :: Ref{ComplexF64}, binom_c :: Ref{Int64}, ComplexF64.(st) :: Ref{ComplexF64}, st_dcp :: Ref{ComplexF64}) :: Nothing
    return st_dcp
end


"""
    GetEntSpec(st :: Vector{<:Number}, bs0 :: SBasis, secd_lst :: Vector{Vector{Vector{Int64}}}, secf_lst :: Vector{Vector{Vector{<:Number}}} ; qnd_a :: Vector{SQNDiag}, qnd_b :: Vector{SQNDiag} = qnd_a, qnf_a :: Vector{SQNOffd}, qnf_b :: Vector{SQNOffd} = qnf_a, amp_oa :: Vector{<:Number}, amp_ob :: Vector{<:Number} = sqrt.(1 .- abs.(amp_oa .^ 2))) :: Dict{@NamedTuple{secd_a, secf_a, secd_b, secf_b}, Vector{Float64}}

# Arguments 

* `st :: Vector{<:Number}` is the state to be decomposed into direct-product basis of two subsystems.
* `bs0 :: SBasis` is the total basis. 
* `secd_lst :: Vector{Vector{Vector{Int64}}}` gives the list of QNDiag sectors of subsystems to be calculated. Each of its elements is a two element vector ; the first specifies the sector for subsystem A, and the second specifies the sector for subsystem B. 
* `secf_lst :: Vector{Vector{Vector{ComplexF64}}}` gives the list of QNOffd sectors of subsystems to be calculated. Each of its elements is a two element vector ; the first specifies the sector for subsystem A, and the second specifies the sector for subsystem B. 
* `qnd_a :: Vector{SQNDiag}, qnd_b :: Vector{SQNDiag} = qnd_a, qnf_a :: Vector{QNOffd}, qnf_b :: Vector{QNOffd}` specifies the diagonal and off-diagonal quantum numbers of the subsystems A and B. `qnd_b` and `qnf_b` are facultative and the same as `qnd_a` and `qnf_a` by default. 
* `amp_oa :: Vector{ComplexF64}` and `amp_ob :: Vector{ComplexF64}` are complex lists of length `no` that specify the amplitute of each orbital in the subsystems A and B. For a non-local basis, we decompose each electron into creation operators in two subsystems ``c^†_o=a_{o,A}c^†_{o,A}+a_{o,B}c^†_{o,B}`` and this list specifies ``a_{o,A}``. This is equivalent to ``√{ℱ_{m,A}}`` in [PRB 85, 125308 (2012)](https://dx.doi.org/10.1103/PhysRevB.85.125308) with an extra phase factor. 

# Output

A dictionary whose keys are named tuples that specify the sector containing entries `secd_a`, `secf_a`, `secd_b`, `secf_b` and values are lists of eigenvalues of the density matrix in those sectors. 
"""
function GetEntSpec(st :: Vector{<:Number}, bs0 :: SBasis, secd_lst :: Vector{Vector{Vector{Int64}}}, secf_lst :: Union{Vector{Vector{Vector{Int64}}}, Vector{Vector{Vector{Float64}}}, Vector{Vector{Vector{ComplexF64}}} } ; qnd_a :: Vector{SQNDiag}, qnd_b :: Vector{SQNDiag} = qnd_a, qnf_a :: Vector{SQNOffd}, qnf_b :: Vector{SQNOffd} = qnf_a, amp_ofa :: Vector{<:Number}, amp_oba :: Vector{<:Number}, amp_ofb :: Vector{<:Number} = sqrt.(1 .- abs.(amp_ofa .^ 2)), amp_obb :: Vector{<:Number} = sqrt.(1 .- abs.(amp_oba .^ 2)), disp_std = false)
    nof = bs0.cfs.nof 
    nob = bs0.cfs.nob 
    norf = bs0.cfs.norf
    norb = bs0.cfs.norb
    nebm = bs0.cfs.nebm 
    dictlock = ReentrantLock()
    EntSpec = Dict{@NamedTuple{secd_a :: Vector{<:Number}, secf_a :: Vector{<:Number}, secd_b :: Vector{<:Number}, secf_b :: Vector{<:Number}}, Vector{Float64}}()
    Threads.@threads for secd in secd_lst 
        cfsa = SConfs(nof, nob, nebm, secd[1], qnd_a ; norf, norb, num_th = 1, disp_std)
        if (cfsa.ncf == 0) continue end 
        cfsb = SConfs(nof, nob, nebm, secd[2], qnd_b ; norf, norb, num_th = 1, disp_std)
        if (cfsb.ncf == 0) continue end 
        for secf in secf_lst 
            bsa = SBasis(cfsa, secf[1], qnf_a ; num_th = 1, disp_std)
            if (bsa.dim == 0) continue end 
            bsb = SBasis(cfsb, secf[2], qnf_b ; num_th = 1, disp_std)
            if (bsb.dim == 0) continue end 
            st_dcp = StateDecompMat(st, bs0, bsa, bsb, amp_ofa, amp_oba, amp_ofb, amp_obb) 
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