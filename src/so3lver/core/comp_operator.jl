export CompOperator, BuildCompOperator
import FuzzifiED: GetEigensystem


"""
    CompOperator{Float64}
    CompOperator{ComplexF64}

The mutable type `CompOperator` represents a composite operator — such as the Hamiltonian — acting on a [CompSpace](@ref CompSpace) of definite total angular momentum. It combines the reduced matrix elements of the per-part [SegOperators](@ref SegOperator) with the ``9j`` re-coupling coefficients that relate the coupled basis of the initial and final composite spaces. The operator is stored in a block-structured form. It can be multiplied formally to a state.

# Fields

* `cpspd :: CompSpace` and `cpspf :: CompSpace` are the initial and final composite spaces.
* `nd :: Int64` is the number of decomposition channels of the operator.
* `ltot :: Int64` is twice the total angular momentum ``2l_{\\text{tot}}`` carried by the operator (``0`` for a scalar such as the Hamiltonian).
* `coeff :: Vector{T}` is the coefficient of each decomposition channel.
* `sgop :: Matrix{SegOperator}` is the ``N_p×N_d`` matrix of segment operators produced by [BuildSegOperators](@ref BuildSegOperators).
* `colptr :: Matrix{Int64}` and `rowid :: Vector{Vector{Int64}}` store the allowed blocks of sectors for each channel. `colptr[:, d]` and `rowid[d]` bear the format of a CSC sparse matrix.
* `idel :: Vector{Matrix{Int64}}` locates the matrix element block in the SegOperators for each segment. It takes three indices `idel[d][p, e]` where `d` is the channel index, `e` is the element index, and `p` is the part index. 
* `mat9j :: Vector{Vector{Matrix{Float64}}}` stores the precomputed ``9j`` re-coupling coefficients, including the fermionic reordering sign. It takes four indices `mat9j[d][e][ich, jch]` where `ich`, `jch` is the channel indices within the sectors specified by `e`.
* `wklist :: Vector{Tuple{Int64, Int64, Int64}}` is the flattened list of the `(jsec, d, e)` matrix element blocks, ordered from the most to the least expensive. The threads of the operator application take their blocks greedily from this list.
"""
mutable struct CompOperator{T <: Union{Float64, ComplexF64}}
    cpspd :: CompSpace
    cpspf :: CompSpace
    nd :: Int64
    ltot :: Int64
    coeff :: Vector{T}
    sgop :: Matrix{SegOperator}
    colptr :: Matrix{Int64}
    rowid :: Vector{Vector{Int64}}
    idel :: Vector{Matrix{Int64}}
    mat9j :: Vector{Vector{Matrix{Float64}}}
    wklist :: Vector{Tuple{Int64, Int64, Int64}}
end


"""
    BuildCompOperator(cpspd :: CompSpace{T}[, cpspf :: CompSpace{T}], cpd :: CoupleDecomps[, sgop :: Matrix{SegOperator}][, ltot :: Int64] ; ident_seg :: Vector{Int64}, full_mat :: Bool, num_th :: Int64) :: CompOperator

constructs a [CompOperator](@ref CompOperator) from the composite spaces, the coupling decompositions `cpd` and the segment operators `sgop`. It computes and stores the ``9j`` re-coupling coefficient between every pair of initial and final coupling channels and every decomposition channel together with the sign arising from fermion parity.

# Arguments

* `cpspd :: CompSpace{T}` is the initial composite space.
* `cpspf :: CompSpace{T}` is the final composite space. Facultative, the same as `cpspd` by default.
* `cpd :: CoupleDecomps` is the list of coupling decompositions ; it must be the same one used to build `sgop`.
* `sgop :: Matrix{SegOperator}` is the matrix of segment operators. Facultative, if omitted, the segment operators will be automatically generated from [`BuildSegOperators`](@ref).
* `ltot :: Int64` is twice the total angular momentum ``2l_{\\text{tot}}`` carried by the operator. Facultative, ``0`` (a scalar) by default.
* `ident_seg :: Vector{Int64}`, `full_mat :: Bool`, and `num_th :: Int64` are forwarded to [`BuildSegOperators`](@ref) ; they are accepted only when `sgop` is omitted. 
* `disp_std :: Bool`, whether or not the log shall be displayed. Facultative, `!SilentStd` by default. 

# Output

* `cpop :: CompOperator` is the resulting composite operator.
"""
function BuildCompOperator(cpspd :: CompSpace{T}, cpspf :: CompSpace{T}, cpd :: CoupleDecomps, sgop :: Matrix{SegOperator}, ltot :: Int64 = 0 ; disp_std = !FuzzifiED.SilentStd) where T <: Union{Float64, ComplexF64}
    nd = length(cpd)
    ch = [ cpd[d].ch for d = 1 : nd ]
    coeff = [ cpd[d].coeff for d = 1 : nd ]
    np = cpspd.np

    colptr = zeros(Int64, size(cpspd.idsec, 2) + 1, nd)
    colptr[1, :] .= 1
    rowid = [ Int64[] for d = 1 : nd ]
    idel = [ Matrix{Int64}(undef, np, 0) for d = 1 : nd]
    Threads.@threads :greedy for d = 1 : nd 
        index = 0
        for jsec in axes(cpspd.idsec, 2)
            idsecj = cpspd.idsec[:, jsec]
            idel_rng = [ sgop[p, d].colptr[idsecj[p]] : sgop[p, d].colptr[idsecj[p] + 1] - 1 for p = 1 : np ]
            for ideli in Iterators.product(idel_rng...)
                idseci = [ sgop[p, d].rowid[ideli[p]] for p = 1 : np]
                isec_rng = searchsorted(eachcol(cpspf.idsec), idseci)
                isempty(isec_rng) && continue
                isec = isec_rng[1]

                index += 1 
                push!(rowid[d], isec)
                idel[d] = hcat(idel[d], collect(ideli))
            end
            colptr[jsec + 1, d] = index + 1
        end
    end

    mat9j = [ Vector{Matrix{Float64}}(undef, colptr[end, d] - 1) for d = 1 : nd ]
    wkcost = [ Vector{Int64}(undef, colptr[end, d] - 1) for d = 1 : nd ]
    Threads.@threads :greedy for (jsec, d) in collect(Iterators.product(axes(cpspd.idsec, 2), 1 : nd))
        idsecj = cpspd.idsec[:, jsec]
        dimj = [ _ChannelDims(cpspd, jsec, j) for j in eachindex(cpspd.chs[jsec]) ]
        for e = colptr[jsec, d] : colptr[jsec + 1, d] - 1
            isec = rowid[d][e]
            dimi = [ _ChannelDims(cpspf, isec, i) for i in eachindex(cpspf.chs[isec]) ]
            pfh = [ mod(cpd[d].sec[1, p], 2) for p = 1 : np ]
            pfj = [ mod(cpspd.sgsp[p].sec[1, idsecj[p]], 2) for p = 1 : np]
            pftot = sum([pfj[p] * sum(pfh[p + 1 : end]) for p = 1 : np]) % 2
            chh = ch[d]
            mat = Matrix{Float64}(undef, length(cpspf.chs[isec]), length(cpspd.chs[jsec]))
            for j in eachindex(cpspd.chs[jsec])
                chj = cpspd.chs[jsec][j]
                for i in eachindex(cpspf.chs[isec])
                    chi = cpspf.chs[isec][i]
                    fac = 1
                    for p = 2 : np 
                        fac = fac * float(d9j(chj[2,p-1], chj[1,p], chj[2,p],  chh[2,p-1], chh[1,p], chh[2,p],  chi[2,p-1], chi[1,p], chi[2,p])) * √((chj[2,p]+1) * (chh[2,p]+1) * (chi[2,p-1]+1) * (chi[1,p]+1))
                    end
                    (pftot == 1) && (fac = -fac)
                    mat[i, j] = fac
                end
            end
            mat9j[d][e] = mat
            cost = 0
            for j in eachindex(cpspd.chs[jsec]), i in eachindex(cpspf.chs[isec])
                (abs(mat[i, j]) < √eps(Float64)) && continue
                cost += _KronMulCost(dimi[i], dimj[j])
            end
            wkcost[d][e] = cost
        end
    end

    wklist = [ (jsec, d, e) for d = 1 : nd for jsec in axes(cpspd.idsec, 2) for e = colptr[jsec, d] : colptr[jsec + 1, d] - 1 ]
    sort!(wklist ; by = wk -> wkcost[wk[2]][wk[3]], rev = true)
    disp_std && @info "FINISH GENERATING COMP OPERATOR"

    return CompOperator{T}(cpspd, cpspf, nd, ltot, coeff, sgop, colptr, rowid, idel, mat9j, wklist)
end

BuildCompOperator(cpspd :: CompSpace{T}, cpd :: CoupleDecomps, sgop :: Matrix{SegOperator}, ltot :: Int64 = 0) where T <: Union{Float64, ComplexF64} = BuildCompOperator(cpspd, cpspd, cpd, sgop, ltot)

function BuildCompOperator(cpspd :: CompSpace{T}, cpspf :: CompSpace{T}, cpd :: CoupleDecomps, ltot :: Int64 = 0 ; ident_seg :: Vector{Int64} = collect(1 : cpspd.np), full_mat :: Bool = false, num_th :: Int64 = FuzzifiED.NumThreads) where T <: Union{Float64, ComplexF64}
    sgop = BuildSegOperators(cpspd.sgsp, cpspf.sgsp, cpd ; full_mat, ident_seg, num_th)
    return BuildCompOperator(cpspd, cpspf, cpd, sgop, ltot)
end

BuildCompOperator(cpspd :: CompSpace{T}, cpd :: CoupleDecomps, ltot :: Int64 = 0 ; ident_seg :: Vector{Int64} = collect(1 : cpspd.np), full_mat :: Bool = false, num_th :: Int64 = FuzzifiED.NumThreads) where T <: Union{Float64, ComplexF64} = BuildCompOperator(cpspd, cpspd, cpd, ltot ; full_mat, ident_seg, num_th)


"""
    *(cpop :: CompOperator{T}, std :: Vector{T} ; num_th :: Int64) :: Vector{T}
    *(stf :: LinearAlgebra.Adjoint{T, Vector{T}}, cpop :: CompOperator{T}, std :: Vector{T} ; num_th :: Int64) :: Vector{T}

applies the composite operator `cpop` to a state `std` of the initial composite space and returns the resulting state of the final composite space or calculates its inner product between an initial and a final state. The action is evaluated block by block : for every decomposition channel and every pair of coupling channels it takes the Kronecker product of the corresponding per-part reduced matrix element blocks, weighted by the channel coefficient and the ``9j`` re-coupling factor. The number of threads used `num_th` is by default `NumThreads`.
"""
function Base.:*(cpop :: CompOperator{T}, std :: Vector{T} ; num_th = FuzzifiED.NumThreads) where T <: Union{Float64, ComplexF64}
    th_lock = ReentrantLock()
    stf = zeros(T, cpop.cpspf.dim)
    np = cpop.cpspd.np

    nwk = length(cpop.wklist)

    maxblk = 1
    for pts in cpop.cpspf.ptr_st, i in 1 : length(pts) - 1
        maxblk = max(maxblk, pts[i + 1] - pts[i])
    end
    nth = max(1, min(num_th, nwk))

    next_wk = Threads.Atomic{Int64}(0)

    nth_blas = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    @sync for _ = 1 : nth
        Threads.@spawn begin
            stf1 = zeros(T, cpop.cpspf.dim)
            scratch = Vector{T}(undef, maxblk)
            idlj = Vector{Int64}(undef, np)
            idli = Vector{Int64}(undef, np)
            blocks = Vector{Matrix{T}}(undef, np)
            scr2 = T[]
            nwk_th = 0
            while true
                iwk = Threads.atomic_add!(next_wk, 1) + 1
                iwk > nwk && break
                nwk_th += 1
                jsec, d, e = cpop.wklist[iwk]
                idsecj = @view cpop.cpspd.idsec[:, jsec]
                coeff = cpop.coeff[d]
                idel_sg = @view cpop.idel[d][:, e]
                isec = cpop.rowid[d][e]
                idseci = @view cpop.cpspf.idsec[:, isec]
                mat9j_e = cpop.mat9j[d][e]
                for jch in eachindex(cpop.cpspd.chs[jsec])
                    jrng = cpop.cpspd.ptr_st[jsec][jch] : cpop.cpspd.ptr_st[jsec][jch + 1] - 1
                    for p = 1 : np
                        lj = cpop.cpspd.chs[jsec][jch][1, p]
                        idlj[p] = cpop.cpspd.sgsp[p].l_lookup[idsecj[p]][lj]
                    end
                    for ich in eachindex(cpop.cpspf.chs[isec])
                        fac9j = mat9j_e[ich, jch]
                        abs(fac9j) < √eps(Float64) && continue
                        irng = cpop.cpspf.ptr_st[isec][ich] : cpop.cpspf.ptr_st[isec][ich + 1] - 1
                        for p = 1 : np
                            li = cpop.cpspf.chs[isec][ich][1, p]
                            idli[p] = cpop.cpspf.sgsp[p].l_lookup[idseci[p]][li]
                        end
                        for p = 1 : np
                            blocks[p] = cpop.sgop[p, d].elmat[idel_sg[p]][idli[p], idlj[p]]
                        end
                        tmp = @view scratch[1 : length(irng)]
                        _KronMul!(tmp, blocks, (@view std[jrng]), scr2)
                        @views stf1[irng] .+= (coeff * fac9j) .* tmp
                    end
                end
            end
            if nwk_th > 0
                lock(th_lock) do
                    stf .+= stf1
                end
            end
        end
    end
    BLAS.set_num_threads(nth_blas)
    return stf
end
Base.:*(stf :: LinearAlgebra.Adjoint{T, Vector{T}}, cpop :: CompOperator{T}, std :: Vector{T}) where T <: Union{Float64, ComplexF64} = stf * (cpop * std)


"""
    Matrix(cpop :: CompOperator{T} ; disp_std :: Bool, num_th :: Int64) :: Matrix{T}

materializes the composite operator `cpop` into a dense matrix. The number of threads used `num_th` is by default `NumThreads`.
"""
function Base.Matrix(cpop :: CompOperator{T} ; disp_std = !FuzzifiED.SilentStd, num_th = FuzzifiED.NumThreads) where T <: Union{Float64, ComplexF64}
    np = cpop.cpspd.np
    mat = zeros(T, cpop.cpspf.dim, cpop.cpspd.dim)

    nwk = length(cpop.wklist)
    nth = max(1, min(num_th, nwk))

    sec_lock = [ ReentrantLock() for _ in axes(cpop.cpspd.idsec, 2) ]
    next_wk = Threads.Atomic{Int64}(0)
    nth_blas = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    @sync for _ = 1 : nth
        Threads.@spawn begin
            idlj = Vector{Int64}(undef, np)
            idli = Vector{Int64}(undef, np)
            blocks = Vector{Matrix{T}}(undef, np)
            while true
                iwk = Threads.atomic_add!(next_wk, 1) + 1
                iwk > nwk && break
                jsec, d, e = cpop.wklist[iwk]
                idsecj = @view cpop.cpspd.idsec[:, jsec]
                coeff = cpop.coeff[d]
                idel_sg = @view cpop.idel[d][:, e]
                isec = cpop.rowid[d][e]
                idseci = @view cpop.cpspf.idsec[:, isec]
                mat9j_e = cpop.mat9j[d][e]
                for jch in eachindex(cpop.cpspd.chs[jsec])
                    jrng = cpop.cpspd.ptr_st[jsec][jch] : cpop.cpspd.ptr_st[jsec][jch + 1] - 1
                    for p = 1 : np
                        lj = cpop.cpspd.chs[jsec][jch][1, p]
                        idlj[p] = cpop.cpspd.sgsp[p].l_lookup[idsecj[p]][lj]
                    end
                    for ich in eachindex(cpop.cpspf.chs[isec])
                        fac9j = mat9j_e[ich, jch]
                        abs(fac9j) < √eps(Float64) && continue
                        irng = cpop.cpspf.ptr_st[isec][ich] : cpop.cpspf.ptr_st[isec][ich + 1] - 1
                        for p = 1 : np
                            li = cpop.cpspf.chs[isec][ich][1, p]
                            idli[p] = cpop.cpspf.sgsp[p].l_lookup[idseci[p]][li]
                        end
                        for p = 1 : np
                            blocks[p] = cpop.sgop[p, d].elmat[idel_sg[p]][idli[p], idlj[p]]
                        end
                        blk = (np == 1) ? blocks[1] : kron(blocks...)
                        lock(sec_lock[jsec]) do
                            @views mat[irng, jrng] .+= (coeff * fac9j) .* blk
                        end
                    end
                end
            end
        end
    end
    BLAS.set_num_threads(nth_blas)
    disp_std && @info "FINISH GENERATING MATRIX OF DIMENSION $(size(mat, 1)) * $(size(mat, 2))"
    return mat
end


function _KronMul!(y :: AbstractVector{T}, As, x :: AbstractVector{T}, scr :: Vector{T}) where T
    N = length(As)
    if N == 1
        mul!(y, As[1], x)
    elseif N == 2
        A1, A2 = As[1], As[2]
        m1, n1 = size(A1) ; m2, n2 = size(A2)
        len = m2 * n1
        (length(scr) < len) && resize!(scr, len)
        X = reshape(x, n2, n1)
        W = reshape(view(scr, 1 : len), m2, n1)
        mul!(W, A2, X)
        mul!(reshape(y, m2, m1), W, transpose(A1))
    else
        copyto!(y, _KronVec(As, x))
    end
    return y
end

# The dimension that each part contributes to the channel `ich` of the sector `isec`
function _ChannelDims(cpsp :: CompSpace, isec :: Int64, ich :: Int64)
    idseci = @view cpsp.idsec[:, isec]
    chi = cpsp.chs[isec][ich]
    return [ begin
        il = cpsp.sgsp[p].l_lookup[idseci[p]][chi[1, p]]
        cpsp.sgsp[p].ptr_st[idseci[p]][il + 1] - cpsp.sgsp[p].ptr_st[idseci[p]][il]
    end for p = 1 : cpsp.np ]
end

# The number of multiplications performed by `_KronMul!`
function _KronMulCost(mpt :: Vector{Int64}, npt :: Vector{Int64})
    np = length(mpt)
    cost = 0
    sufm = 1
    for p = np : -1 : 1
        pren = 1
        for q = 1 : p - 1
            pren *= npt[q]
        end
        cost += pren * sufm * mpt[p] * npt[p]
        sufm *= mpt[p]
    end
    return cost
end

function _KronVec(As, x :: AbstractVector{T}) where T
    (length(As) == 1) && return As[1] * x
    A1 = As[1] ; m1, n1 = size(A1)
    rest = @view As[2 : end]
    nrest = prod(size(B, 2) for B in rest)
    mrest = prod(size(B, 1) for B in rest)
    X = reshape(x, nrest, n1)
    W = Matrix{T}(undef, mrest, n1)
    for c in 1 : n1
        @views W[:, c] = _KronVec(rest, X[:, c])
    end
    return vec(W * transpose(A1))
end


"""
    GetEigensystem(cpop :: CompOperator{T}, nst :: Int64[, alg :: Type{<:KrylovKit.KrylovAlgorithm}] ; tol :: Float64, ncv :: Int64, initvec :: Vector{T}, full_mat :: Bool, proj_sym :: Function, num_th :: Int64, kwargs...) :: Tuple{Vector{T}, Matrix{T}}

computes the lowest `nst` eigen-values and eigen-states of the composite operator `cpop` through `KrylovKit.eigsolve`. 

# Arguments

* `cpop :: CompOperator{T}` is the composite operator (usually the Hamiltonian).
* `nst :: Int64` is the number of eigen-states to compute.
* `alg :: Type{<:KrylovKit.KrylovAlgorithm}`. Facultative. It specifies an algorithm that is directly fed into `KrylovKit.eigsolve`. _E. g._ `LockedLanczos` in the [modified version of KrylovKit].
* `tol :: Float64` is the tolerance of the eigen-solver. Facultative, `1E-8` by default.
* `ncv :: Int64` is the dimension of the Krylov subspace. Facultative, `max(2 * nst, nst + 10)` by default.
* `initvec :: Vector{T}` is the initial vector. Facultative, a random vector after projection by default.
* `full_mat :: Bool`, whether the operator is first materialized into a dense matrix and this matrix is handed to `eigsolve`. Facultative, `false` by default.
* `proj_sym :: Function`. To implement a off-diagonal symmetry, feed in a function from a vector to a vector that projects it to the desired sector. Facultative, `identity` by default. Only supported when `full_mat` is `false`. _N. b._, some unphysical null states may appear. _E. g._, to target the sector that is even under the permutation of the first two segments, feed in `st -> (st .+ PermFirstSecondSegs(cpsp, st)) ./ 2`. 
* `num_th :: Int64` is the number of threads used in matrix multiplication. Facultative, `NumThreads` by default.
* `kwargs...` are further key-word arguments forwarded to `eigsolve`, _e. g._, `ishermitian = true` for complex and `issymmetric = true` for real matrix. 

# Output

* `eigval :: Vector{T}` is the vector of the `nst` lowest eigen-values.
* `eigvec :: Matrix{T}` is the matrix whose columns are the corresponding eigen-states.
"""
function FuzzifiED.GetEigensystem(cpop :: CompOperator{T}, nst :: Int64 ; tol :: Float64 = 1E-8, ncv :: Int64 = max(2 * nst, nst + 10), proj_sym :: Function = identity, initvec = proj_sym(rand(T, cpop.cpspd.dim)), full_mat :: Bool = false, num_th = FuzzifiED.NumThreads, disp_std = !FuzzifiED.SilentStd, kwargs...) where T <: Union{ComplexF64,Float64}
    verbosity = disp_std ? 3 : 0
    if full_mat
        fmul = Matrix(cpop ; disp_std)
    else
        fmul = x -> *(cpop, proj_sym(x) ; num_th)
    end
    eigval, eigvec, info = eigsolve(fmul, initvec, nst, :SR ; tol, krylovdim = ncv, verbosity, kwargs...)
    return Vector{T}(eigval), Matrix{T}(hcat(eigvec...))
end

function FuzzifiED.GetEigensystem(cpop :: CompOperator{T}, nst :: Int64, alg :: Type{<:KrylovKit.KrylovAlgorithm} ; tol :: Float64 = 1E-8, ncv :: Int64 = min(max(2 * nst, nst + 10), 100), proj_sym :: Function = identity, initvec = proj_sym(rand(T, cpop.cpspd.dim)), full_mat :: Bool = false, num_th = FuzzifiED.NumThreads, disp_std = !FuzzifiED.SilentStd, kwargs...) where T <: Union{ComplexF64,Float64}
    verbosity = disp_std ? 3 : 0
    if full_mat
        fmul = Matrix(cpop ; disp_std)
    else
        fmul = x -> *(cpop, proj_sym(x) ; num_th)
    end
    eigval, eigvec, info = eigsolve(fmul, initvec, nst, :SR, alg(; tol, krylovdim = ncv, verbosity, kwargs...))
    return Vector{T}(eigval), Matrix{T}(hcat(eigvec...))
end