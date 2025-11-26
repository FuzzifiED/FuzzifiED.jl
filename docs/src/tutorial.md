# FuzzifiED Explained in a Tutorial

To demonstrate the usage of FuzzifiED interfaces for ED and DMRG, in this section, we use a tutorial that calculates the eigenstates for the Ising model on the fuzzy sphere. Specifically, it

1. calculates the lowest eigenstates in the symmetry sector $L^z=0$ and $(\mathcal{P},\mathcal{Z},\mathcal{R})=(+,+,+)$,
2. measures their total angular momenta, and 
3. calcultes the OPE coefficient $f_{\sigma\sigma\epsilon}=\langle \sigma|n^z_{00}|\epsilon\rangle/\langle \sigma|n^z_{00}|0\rangle$.

Four versions of the tutorial code are provided : 

1. [`tutorial_ising.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/tutorial_ising.jl) — the ED code that uses the built-in models.
2. [`tutorial_ising_primitive.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/tutorial_ising_primitive.jl) — The ED code that uses only the core functions.
3. [`tutorial_ising_dmrg.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/tutorial_ising_dmrg.jl) — the DMRG code that converts the format into ITensor.
4. [`tutorial_ising_dmrg_easysweep.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/tutorial_ising_dmrg_easysweep.jl) — the DMRG code that uses the Easy-Sweep extension.

The examples can be found in the directory [`examples`](https://github.com/FuzzifiED/FuzzifiED.jl/tree/main/examples). We also append in the end [a list of given examples](@ref List-of-Examples) at the end of the page. 

## ED with FuzzifiED

In this section, we briefly describe the procedure for exact diagonalisation (ED) calculation and give an instruction for using FuzzifiED for ED. 

Practically, the ED calculation can be divided into 4 steps.

* Construct a many-body basis that respect a given set of quantum numbers. Specifically, in FuzzifiED we support quantum numbers of commuting $\mathrm{U}(1)$ or discrete $\mathbb{Z}_p$ symmetries.
* Construct the sparse matrix corresponding to the Hamiltonian in the basis above. 
* Find the lowest eigenstates and their corresponding eigenenergies of the sparse matrix.
* Making measurements on the eigenstates. This including the total angular momentum, density operators, entanglement, _etc._

### Setup

Before starting the calculation, we need to input the setup for the system, including the number of flavours $N_f$, orbitals $N_m$ and sites $N_o$ 

* A 'flavour' is labelled by $f$. The number of flavours is $N_f$.
* An 'orbital' is specified by the magnetic quantum number labelled by $m$. The number of orbitals is $N_m=2s+1$.
* A 'site' is specific by both the flavour and the orbital index $o=(f,m)$. The number of sites is $N_o=N_mN_f$. In practice, we label the sites with an integer from $1$ to $N_o$. We store the sites in an ascending order of first $m$ and then $f$ : $o=(m+s)N_f+f$.

In the example of Ising model with $s=5.5$,
```julia
nm = 12
nf = 2
no = nm * nf
```

FuzzifiED also provides several environment parameters that defines how FuzzifiED works, _viz._ [FuzzifiED.ElementType](@ref), [FuzzifiED.NumThreads](@ref) and [FuzzifiED.SilentStd](@ref).

### Constructing the Configurations

The first step for the ED calculation is to construct the basis tha respects the symmetries of the Hamiltonian. This is divided into two steps : (1) generate the 'configurations' that carry the diagonal quantum numbers, and (2) generate the 'basis' that also carry the off-diagonal quantum numbers (under discrete transformations). The '_configurations_' are the collection of states that can be written as direct product of occupied $|1\rangle$ or empty $|0\rangle$ on each site and carries certain diagonal quantum numbers (QNDiag). 

The QNDiags supported by FuzzifiED are the charges of $\mathrm{U}(1)$ or $\mathbb{Z}_p$ symmetry in the form of 
```math
\begin{aligned}
    Q&=\sum_oq_on_o&\mathrm{U}(1)&\textrm{ symmetry}\\
    Q&=\sum_oq_on_o\mod p&\mathbb{Z}_p&\textrm{ symmetry}
\end{aligned}
```
where $n_o=c^\dagger_oc_o$ is the particle number on each site, and $q_o$ is the charge that each orbital carries. FuzzifiED restricts $q_o$ to be integer-valued. In FuzzifiED, the QNDiags are recorded in the mutable type [`QNDiag`](@ref). Several useful QNDiags are [built-in](@ref Diagonal-Quantum-Numbers-on-the-Fuzzy-Sphere). 

The collection of configurations is generated from the QNDiags. It is recorded in the mutable type [`Confs`](@ref) 
and can be constructed by the method 
```julia
Confs(no :: Int64, secd :: Vector{Int64}, qnd :: Vector{QNDiag})
```
where `qnd` is the array of QNDiags, and `secd` is the array of charges $Q$ of each QNDiag. Here each configuration is stored as a binary number with $N_o$ bits. If the $o$-th site in the configuration is occupied, the $(o-1)$-th bit of the configuration is $1$; if the site is empty, then the bit is $0$. Besides the storation of the configuration, we also need a reverse look-up process that returns the index from the binary string. This is realised by a Lin table.

In the example of Ising model, there are two QNDiags, _viz._ the particle number and the angular momentum. 
```math
\begin{aligned}
Q_1&=N_e,& q_{1,m\sigma}&=1\\
Q_2&=2L_z,&q_{2,m\sigma}&=2m
\end{aligned}
```
The full code to generating the configurations in the $L_z=0$ sector is 
```julia
qnd = [ 
    QNDiag(fill(1, no)), 
    QNDiag([ 2 * m - nm - 1 for m = 1 : nm for f = 1 : nf ])
]
cfs = Confs(no, [nm, 0], qnd)
```
Alternatively, using the built-in models, 
```julia
qnd = [ 
    GetNeQNDiag(no), 
    GetLz2QNDiag(nm, nf) 
]
cfs = Confs(no, [nm, 0], qnd)
```

### Constructing the Basis

Having constructed the configurations, we now construct the basis of the Hilbert space. The `_basis_' is the collection of states that are linear combinations of the configuration carrying certain diagonal and $\mathbb{Z}_p$ off-diagonal quantum numbers (QNOffd). 

The QNOffds supported by FuzzifiED are the $\mathbb{Z}_p$ symmetry that are in the form of 
```math
    \mathcal{Z}:\ c_o\to \alpha_o^* c^{(p_o)}_{\pi_o},\quad c_o^\dagger\to \alpha_o c^{(1-p_o)}_{\pi_o}
```
where we use a notation $c^{(1)}=c^\dagger$ and $c^{(0)}=c$ for convenience, $\pi_o$ is a permutation of the sites $1,\dots N_o$, $\alpha_o$ is a coefficient, and $p_o$ specified whether or not particle-hole transformation is performed for the site. Note that one must guarentee that all these transformations commute with each other and also commute with the diagonal quantum numbers. In FuzzifiED, the QNOffds are recorded in the mutable type [`QNOffd`](@ref). Several useful QNOffds are [built-in](@ref Off-diagonal-Quantum-Numbers-on-the-Fuzzy-Sphere).

After implementing the QNOffds, a state in the new basis should look like 
```math
    |I\rangle=\lambda_{i_{I1}}|i_{I1}\rangle+\lambda_{i_{I2}}|i_{I2}\rangle+\cdots+\lambda_{i_{Im_I}}|i_{Im_I}\rangle
```
where the $|i\rangle$'s are configurations, and $|I\rangle$ is a linear combination of them. This process can be regarded as organising the configurations into groups of size $m_I$. 

In FuzzifiED, the basis $\{|I\rangle\}$ is recorded in the mutable type [`Basis`](@ref). It can be constructed by the methods 
```julia
Basis(cfs :: Confs, secf :: Vector{ComplexF64}, qnf :: Vector{QNOffd})
Basis(cfs :: Confs)
```
where `secf` records the eigenvalue of each transformation, typically in the form $e^{i2\pi q/p}$ where $p$ is the cycle and $q$ is the $\mathbb{Z}_p$ charge. 

In the example of Ising model, There are three $\mathbb{Z}_2$ symmetries, _viz._ the particle-hole transformation $\mathcal{P}$, the $\pi$-rotation along the $y$-axis $\mathcal{R}_y$, and the flavour (Ising) symmetry $\mathcal{Z}$
```math
\begin{aligned}
    \mathscr{P}:c^\dagger_{\sigma m}&\to\sigma c_{-\sigma,m}\\
    \mathscr{Z}:c^\dagger_{\sigma m}&\to c^\dagger_{-\sigma,m}\\
    \mathscr{R}_y:c^\dagger_{\sigma m}&\to c^\dagger_{\sigma,-m}\\
\end{aligned}
```
The code to generate the basis in the all-positive sector is 
```julia
qnf = [
    # Parity (Particle-hole)
    QNOffd([ (m - 1) * nf + [2, 1][f] for m = 1 : nm for f = 1 : nf ], true, 
        ComplexF64[ [-1, 1][f] for m = 1 : nm for f = 1 : nf ]),
    # Flavour symmetry
    QNOffd([ (m - 1) * nf + [2, 1][f] for m = 1 : nm for f = 1 : nf ]),
    # Y-axis pi-rotation
    QNOffd([ (nm - m) * nf + f for m = 1 : nm for f = 1 : nf], 
        ComplexF64[ iseven(m) ? 1 : -1 for m = 1 : nm for f = 1 : nf ])
]
bs = Basis(cfs, [1, 1, 1], qnf) 
```
Alternatively, using the built-in functions
```julia
qnf = [ 
    GetParityQNOffd(nm, 2, [2, 1], [-1, 1]), 
    GetFlavPermQNOffd(nm, 2, [2, 1]), 
    GetRotyQNOffd(nm, 2) 
]
bs = Basis(cfs, [1, 1, 1], qnf)
# The second argument gives the eigenvalues under the transformations, for Z_2 put 1,-1 ; for Z_n put exp(2im*pi*q/p)
```

### Recording the Many-Body Operator Terms

Having constructed the basis, we now construct the many-body operators. A general many-body operator can be written as
```math
    \mathscr{O}=\sum_{t=1}^{N_t}U_tc^{(p_{t1})}_{o_{t1}}c^{(p_{t2})}_{o_{t2}}\dots c^{(p_{tl_t})}_{o_{tl_t}}
```
where $c^{(0)}=c$ and $c^{(1)}=c^\dagger$. In FuzzifiED, this is recorded as an array of `Term`, and each `Term` records the building block $Uc^{(p_{1})}_{o_{1}}c^{(p_{2})}_{o_{2}}\dots c^{(p_{l})}_{o_{l}}$. It can be initialised by the method 
```julia
Term(coeff :: ComplexF64, cstr :: Vector{Int64})
```
The addition and multiplication of terms are supported, and the terms can be simplified by the method [`SimplifyTerms`](@ref)
```julia
SimplifyTerms(tms :: Vector{Term})
```
In FuzzifiED, several useful operator terms are [built-in](@ref Operators-on-the-Fuzzy-Sphere).

In the example of Ising model, the full code that records the Hamiltonian is 
```julia
using WignerSymbols
# Input the parameters of the Hamiltonian
ps_pot = [ 4.75, 1. ] * 2.
h = 3.16
tms_hmt = Term[]
# Go through all the m1-up, m2-down, m3-down, m4-up and m4 = m1 + m2 - m
m = zeros(Int64, 4)
for m[1] = 0 : nm - 1, m[2] = 0 : nm - 1, m[3] = 0 : nm - 1
    m[4] = m[1] + m[2] - m[3]
    (m[4] < 0 || m[4] >= nm) && continue
    f = [0, 1, 1, 0]
    o = m .* nf .+ f .+ 1
    mr = m .- s
    
    # Calculate the matrix element val from pseudopotentials
    val = ComplexF64(0)
    for l in eachindex(ps_pot)
        (abs(mr[1] + mr[2]) > nm - l || abs(mr[3] + mr[4]) > nm - l) && break 
        val += ps_pot[l] * (2 * nm - 2 * l + 1) * wigner3j(s, s, nm - l, mr[1], mr[2], -mr[1] - mr[2]) * wigner3j(s, s, nm - l, mr[4], mr[3], -mr[3] - mr[4])
    end 
    # Record the interaction term val, "Cdag", o1, "Cdag", o2, "C", o3, "C", o4
    tms_hmt += Terms(val, [1, o[1], 1, o[2], 0, o[3], 0, o[4]])
end 
for m = 0 : nm - 1
    o = m * nf .+ [1, 2]
    # Record the transverse field term
    tms_hmt += Terms(-h, [1, o[1], 0, o[2]])
    tms_hmt += Terms(-h, [1, o[2], 0, o[1]])
end
```
Alternatively, using the built-in functions
```julia
σ1 = [ 1 0 ; 0 0 ]
σ2 = [ 0 0 ; 0 1 ]
σx = [ 0 1 ; 1 0 ]
σz = [ 1 0 ; 0 -1]
ps_pot = [ 4.75, 1.0 ] * 2.0
fld_h = 3.16
tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, 2, ps_pot, σ1, σ2)
    - fld_h * GetPolTerms(nm, 2, σx) 
)
```

We also need to construct the total angular momentum. It is defined as 
```math
    L^2=L^+L^-+(L^z)^2-L^z,
```
as $c_m$ carries the $\mathrm{SO}(3)$ spin-$s$ representation, 
```math
    L^z=\sum_{mf}mc_m^\dagger c_m,\quad L^\pm=\sum_{mf}\sqrt{(s\mp m)(s\pm m+1)}c^\dagger_{m\pm 1}c_m
```
we can first construt its building blocks and use the addition and multiplication of the terms
```julia
tms_lz = [ Term(m - s - 1, [1, (m - 1) * nf + f, 0, (m - 1) * nf + f]) for m = 1 : nm for f = 1 : nf ]
tms_lp = [ Term(sqrt((nm - m) * m), [1, m * nf + f, 0, (m - 1) * nf + f]) for m = 1 : nm - 1 for f = 1 : nf ]
tms_lm = tms_lp' 
tms_l2 = SimplifyTerms(tms_lz * tms_lz - tms_lz + tms_lp * tms_lm)
```
Alternatively, using the built-in functions,
```julia
tms_l2 = GetL2Terms(nm, nf)
```

### Generating the Sparse Matrix

Having gotten the terms in the many-body operator, we now need to generate the matrix elements given the initial and final basis and find its eigenstates. 

In FuzzifiED, the mutable type [`Operator`](@ref) records the terms together with information about its symmetry and the basis of the state it acts on and the basis of the resulting state. 
be initialised with the [method](@ref Operator(bsd :: Basis, bsf :: Basis, terms :: Terms ; red_q :: Int64 = 0, sym_q :: Int64 = 0)) 
```julia
Operator(bsd :: Basis[, bsf :: Basis], terms :: Vector{Term} ; red_q :: Int64, sym_q :: Int64)
```
In FuzzifiED, the sparse matrix is stored in the mutable type [`OpMat{T}`](@ref OpMat) where `T` is the type of the elements (`ComplexF64` or `Float64`). It can be generated from the method
```julia
OpMat[{T}](op :: Operator)
```

### Finding the Eigenstates

After generating the sparse matrix, the method [`GetEigensystem`](@ref) uses the Fortran Arpack package to calculate its lowest eigenstates. 

In the example of Ising model, the full code to calculate the lowest $N_\mathrm{st}=10$ eigenstates from the basis and the terms is 
```julia
nst = 10
hmt = Operator(bs, tms_hmt)
hmt_mat = OpMat(hmt)
enrg, st = GetEigensystem(hmt_mat, nst)
```

### Write the Sparse Matrix into a HDF5 File 

It is sometimes needed to write the sparse matrix into file to avoid extra effort to generate it again in another calculation. 

```julia
using HDF5 
f = h5open("data.h5", "cw")
write(f, "hmt_mat", hmt_mat)
close(f)
```

To read from the file

```julia
f = h5open("data.h5", "r") 
hmt_mat = read(f, "hmt_mat", OpMat{ComplexF64})
close(f)
```

Apart from `OpMat`, the supported types for writing include `Confs`, `Basis`, `Terms`, `Operator`, `OpMat{ComplexF64}` and `OpMat{Float64}`. 

### Inner Product of States, Operators and Transformations

Having obtained the eigenstates, we need to make measurements on it. The simplest kind of measurements is the inner product of a many body operator with two states $\langle j|\mathcal{O}|i\rangle$. FuzzifiED supports the inner product and vector product of `Operator` and `OpMat{T}` with vectors that represent the state
```julia
(op :: Operator) * (st_d :: Vector{T}) :: Vector{T}
(mat :: OpMat{T}) * (st_d :: Vector{T}) :: Vector{T}
(st_f :: Vector{T}) * (op :: Operator) * (st_d :: Vector{T}) :: T
(st_f :: Vector{T}) * (mat :: OpMat{T}) * (st_d :: Vector{T}) :: T
```

For example, the code to measure the angular momenta of each state is
```julia
tms_l2 = GetL2Terms(nm, 2)
l2 = Operator(bs, tms_l2)
l2_mat = OpMat(l2)
l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)]
```

### Go through All the Sectors

One can repeat the calculation for all the sectors and records the results
```julia
result = []
for P in [1, -1], Z in [1, -1], R in [1, -1]
    bs = Basis(cfs, [P, Z, R], qnf)
    hmt = Operator(bs, tms_hmt)
    hmt_mat = OpMat(hmt)
    enrg, st = GetEigensystem(hmt_mat, 10)

    l2 = Operator(bs, tms_l2)
    l2_mat = OpMat(l2)
    l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)]

    for i in eachindex(enrg)
        push!(result, round.([enrg[i], l2_val[i], P, Z], digits = 6))
    end
end
```
We then sort the eigenstates, find the energy of ground state and stress tensor, and calibrate the scaling dimensions. 
```julia
sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_T = filter(st -> st[2] ≊ 6 && st[3] ≊ 1 && st[4] ≊ 1, result)[1][1]
result_dim = [ [ 3 * (st[1] - enrg_0) / (enrg_T - enrg_0) ; st] for st in result ]
display(permutedims(hcat(result_dim...)))
```

### Measuring Local Observables

Local observables are a kind of particularly useful operators on the fuzzy sphere. Their value at a point on the sphere can be decomposed into spherical components, and the multiplication of the components follows the triple integral formula of monopole spherical harmonics
```math
\begin{aligned}
    \mathcal{O}(\hat{\mathbf{n}})&=\sum_{lm}Y^{(s)}_{lm}(\hat{\mathbf{n}})\mathcal{O}_{lm}\\
    (\mathcal{O}_1\mathcal{O}_2)_{lm}&=\sum_{l_1l_2m_1m_2}(\mathcal{O}_1)_{l_1m_1}(\mathcal{O}_2)_{l_2m_2}\\
    &\qquad\qquad\times(-1)^{s+m}\sqrt{\frac{(2l_1+1)(2l_2+1)(2l_3+1)}{4\pi}}\begin{pmatrix}l_1&l_2&l\\m_1&m_2&-m\end{pmatrix}\begin{pmatrix}l_1&l_2&l\\-s_1&-s_2&s\end{pmatrix}
\end{aligned}
```
In FuzzifiED, they are stored in the type [`SphereObs`](@ref) and can be initialised from either a function or a dictionary that specifies the components
```julia
SphereObs(s2 :: Int64, l2m :: Int64, get_comp :: Function)
SphereObs(s2 :: Int64, l2m :: Int64, comps :: Dict)
```
Their adjoint, addition, multiplication and [Laplacian](@ref) are supported. The related functions are [`StoreComps`](@ref) that stores all the components, [`GetComponent`](@ref) and [`GetPointValue`](@ref) that evaluate a spherical component $\mathcal{O}_{lm}$ or value at one point $\mathcal{O}(\hat{\mathbf{n}})$. Several important types of spherical observables are built-in in FuzzifiED, _viz._, [electron](@ref GetElectronObs), [density operator](@ref GetDensityObs) and [pairing operator](@ref GetPairingObs).

In the example of Ising model, to calculate the OPE coefficient $f_{\sigma\sigma\epsilon}=\langle \sigma|n^z_{00}|\epsilon\rangle/\langle \sigma|n^z_{00}|0\rangle$, one need to first calculate the eigenstates in the $\mathbb{Z}_2$-odd sector
```julia
bs_m = Basis(cfs, [1, -1, 1], qnf)
hmt_m = Operator(bs_m, bs_m, tms_hmt ; red_q = 1, sym_q = 1) 
hmt_mat_m = OpMat(hmt_m)
enrg_m, st_m = GetEigensystem(hmt_mat_m, 10)
st0 = st[:, 1] 
ste = st[:, 2] 
sts = st_m[:, 1]
```
and then construct the density operator
```julia
obs_nz = GetDensityObs(nm, 2, σz)
tms_nz00 = SimplifyTerms(GetComponent(obs_nz, 0.0, 0.0))
nz00 = Operator(bs, bs_m, tms_nz00 ; red_q = 1) 
f_sse = abs((sts' * nz00 * ste) / (sts' * nz00 * st0))
```

Besides the spherical observable, we also provide a type [`AngModes`](@ref) that superposes under the rule of angular momentum superposition instead of spherical harmonics triple integral
```math
    (\mathscr{A}_1\mathscr{A}_2)_{lm}=\sum_{l_1m_1l_2m_2}(\mathcal{A}_1)_{l_1m_1}(\mathcal{A}_2)_{l_2m_2}\langle l_1m_1l_2m_2|lm\rangle.
```
The interfaces are similar.

### Measuring the Entanglement

A non-local quantity that bears particular significance is the entanglement. To calculate the entanglement, we divide the sphere into two parts $A$ and $B$. The reduced density matrix of part $A$ is obtained by tracing the density matrix over the part $B$
```math
    \rho_A(\Psi)=\operatorname{tr}_B|\Psi\rangle\langle\Psi|
```
The entanglement entropy is $S=-\operatorname{tr}\rho_A\log\rho_A$ and the entanglement spectrum is the collection of eigenvalues of $\rho_A$ taken negative logarithm. 

The detail of the calculation is given in [PRB 85, 125308 (2012)](https://dx.doi.org/10.1103/PhysRevB.85.125308). Here we only sketch the process. The creation operator in each orbital is divided into the creation on $A$ part and the creation on $B$ part. 
```math
    c^\dagger_o=\alpha_oc^\dagger_{o,A}+\beta_mc^\dagger_{o,B}
```
where $|\alpha_o|^2+|\beta_o|^2=1$. For the cut in orbital space $m_c$, 
```math
    \alpha_{mf}=\Theta(m_c-m)
```
where $\Theta$ is the Heaviside function ; for the cut in real space along latitude circle $\theta_c$,
```math
    \alpha_{mf}=\Beta_{\cos^2\theta_c/2}(s-m+1,s+m+1)^{1/2}
```
where $\Beta$ is the incomplete beta function. 

To calculate the reduced density matrix, we decompose the state into the direct-product basis of two subsystems
```math
    |\Psi\rangle=\sum_{K_0}v_{K_0}|K_0\rangle=\sum_{I_AJ_B}M_{I_AJ_B}|I_A\rangle|J_B\rangle
```
where the indices $K_0\in\mathcal{H},I_A\in\mathscr{H_A},J_B\in\mathscr{H_B}$ are in the overall Hilbert space and the Hilbert space of subsystem $A$ and $B$. The density matrix is then 
```math
    \rho_A=\mathbf{M}\mathbf{M}^\dagger
```
and the entanglement spectrum can be obtained from the SVD decomposition of the $\mathbf{M}$ matrix. Like the Hamiltonian, the $\mathbf{M}$ matrix is block diagonal, and each block carries different quantum numbers of the Hilbert spaces of $A$ and $B$ subsystem. 

In FuzzifiED, the decomposition of states into matrix $M_{I_AJ_B}$ is done by the funciton [StateDecompMat](@ref), and the calculation of entanglement spectrum is done by the funciton [GetEntSpec](@ref).

In the example of Ising model, to calculate the entanglement entropy cut from the equator, we first need to specify the quantum numbers of the subsystems : the conservation of $N_e$, $L_z$ and the $\mathbb{Z}_2$ symmetry.
```julia
qnd_a = [ GetNeQNDiag(no), GetLz2QNDiag(nm, nf) ]
qnf_a = [ GetFlavPermQNOffd(nm, nf, [2, 1]) ]
```
we then specify the sectors to calculate : The number of electrons in subsystem $A$ run from $0$ to $N_m$ ; the angular momenta in subsystem $A$ can take all permitted values ; for subsystem $B$, $N_{e,B}=N_m-N_{e,A}$, $L_{z,B}=-L_{z,A}$ ; the $\mathbb{Z}_2$ sectors of the two subsystems are the same. 
```julia
secd_lst = Vector{Vector{Int64}}[]
for nea = 0 : nm 
    neb = nm - nea 
    for lza = -min(nea, neb) * (nm - 1) : 2 : min(nea, neb) * (nm - 1)
        lzb = -lza 
        push!(secd_lst, [[nea, lza], [neb, lzb]])
    end
end
secf_lst = [ [[1], [1]], [[-1], [-1]] ]
```
Finally, we specify the list of amplitute $\alpha_m$.
```julia
amp_oa = [ sqrt(beta_inc(m, nm - m + 1, 0.5)) for f = 1 : 2 for m = 1 : nm]
```
To calculate the entanglement spectrum, 
```julia
ent_spec = GetEntSpec(st_g, bs, secd_lst, secf_lst ; qnd_a, qnf_a, amp_oa)
```
The entanglement entropy can be calculated by collecting all the eigenvalues of the density matrix.
```julia
eig_rho = vcat(values(ent_spec)...)
ent_entropy = -sum(eig_rho .* log.(eig_rho))
```

## DMRG with FuzzifiED

Having introduced ED, we now turn to density matrix renormalisation group (DMRG) that deals with larger systems. We briefly describe its procedure  and give an instruction for using FuzzifiED for DMRG. 

Practically, the `dmrg` function in ITensor package automatically uses DMRG to optimise a matrix product state (MPS) to be the lowest eigenstate of a Hermitian Hamiltonian represented as a matrix product operator (MPO). To generate the input of the function, one needs to 

* construct a set of sites that carries a certain set of QNDiags,
* construct a MPO representing the Hamiltonian on the sites from a set of terms (or OpSum in ITensor), and 
* construct an initial MPS on the sites in the desired symmetry sector.

In FuzzifiED, a new SiteType [`"FuzzyFermion"`](@ref ITensors.space) is defined that behaves similar to the built-in `"Fermion"` type and the set of sites can be generated by the function [`GetSites`](@ref).

In the example of Ising model, for convenience we exchange the Pauli matrices $\sigma^x$ and $\sigma^z$ so that the two flavours carry $\mathbb{Z}_2$-charge $0$ and $1$. The sites can be constructed by 
```julia
nm = 12
nf = 2
no = nm * nf
sites = GetSites([
    GetNeQNDiag(nm * nf), 
    GetLz2QNDiag(nm, nf),
    GetZnfChargeQNDiag(nm, nf)
])
```

In ITensor, the MPO is generated from an OpSum and the sites. The [OpSum](@ref) can be directly converted from the array of terms. In the example of Ising model, 
```julia
σx = [  0  1 ;  1  0 ]
σz = [  1  0 ;  0 -1 ]
ps_pot = [4.75, 1.] ./ 2
tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, 2, ps_pot) - 
    GetDenIntTerms(nm, 2, ps_pot, σx) - 
    3.16 * GetPolTerms(nm, nf, σz)
)
os_hmt = OpSum(tms_hmt)
hmt = MPO(os_hmt, sites)
```

To calculate the $\mathbb{Z}_2$-even $L^z=0$ sector, the initial state can be taken as the all the $\mathbb{Z}_2$-even sites being filled and all the $\mathbb{Z}_2$-odd sites being empty. (Note that ITensor takes the string `"1"` instead of the number `1` as occupied and `"0"` instead of `0` as filled.) 
```julia
cfi_p = [ [1, 0][f] for m = 1 : nm for f = 1 : nf ]
sti_p = MPS(sites, string.(cfi_p))
```

Having these ingrediants ready, we can call the `dmrg` function. To ensure performance, the maximal bond dimension should be increased gradually and the noise decreased gradually to 0. An example that deals with maximal bond dimension 500 is 
```julia
E0, st0 = dmrg(hmt, sti_p ; 
    nsweeps = 10, 
    maxdim = [10,20,50,100,200,500], 
    noise = [1E-4,3E-5,1E-5,3E-6,1E-6,3E-7], 
    cutoff = [1E-8])
```
To generate a $\mathbb{Z}_2$-odd initial state, we can simply flip the spin on the first orbital
```julia
cfi_m = [ m == 1 ? [0, 1][f] : [1, 0][f] for m = 1 : nm for f = 1 : nf ]
sti_m = MPS(sites, string.(cfi_m))
Es, sts = dmrg(hmt, sti_m ; 
    nsweeps = 10, 
    maxdim = [10,20,50,100,200,500], 
    noise = [1E-4,3E-5,1E-5,3E-6,1E-6,3E-7], 
    cutoff = [1E-8])
```
The first excited $\mathbb{Z}_2$-even state can be generated by adding a projector $w|0\rangle\langle0|$ to the MPO 
```julia
Ee, ste = dmrg(hmt, [st0], sti_p ; 
    nsweeps = 10, 
    maxdim = [10,20,50,100,200,500], 
    noise = [1E-4,3E-5,1E-5,3E-6,1E-6,3E-7], 
    cutoff = [1E-8], 
    weight = 100)
```

The inner product can be measured by the ITensor function `inner`. For example, to measure the angular momentum $L^2$ of the ground state,
```julia
tms_l2 = GetL2Terms(nm, 2)
l2 = MPO(OpSum(tms_l2), sites)
val_l20 = inner(st0', l2, st0)
```
To measure the OPE coefficient $f_{\sigma\sigma\epsilon}=\langle \sigma|n^x_{00}|\epsilon\rangle/\langle \sigma|n^x_{00}|0\rangle$. (Note that the indices $x$ and $z$ have already been exchanged here.)
```julia
obs_nx = GetDensityObs(nm, 2, σx)
tms_nx00 = SimplifyTerms(GetComponent(obs_nx, 0.0, 0.0))
nx00 = MPO(OpSum(tms_nx00), sites)
f_sse = abs(inner(sts', nx00, ste) / inner(sts', nx00, st0))
```

### The Easy-Sweep Extension

The extension Easy-Sweep facilitates the management of DMRG process. It automatically records the intermediate results and recover these results if a job is stopped and run again on HPC. It also manages the gradual increase of maximal bond dimensions and the determination of convergence by the criteria of energy. This extension contains the following functions : [`GetMPOSites`](@ref), [`GetMPO`](@ref), [`SweepOne`](@ref), [`EasySweep`](@ref).

To use the this extension, one need to use the packages `ITensors`, `ItensorMPS` and `HDF5`. A path need to be created a priori to store the result HDF5 files. We recommend using the package `ITensorMPOConstruction` to generate the MPO, which can be installed through 
```julia
using Pkg ; Pkg.add(url="https://github.com/ITensor/ITensorMPOConstruction.jl.git"); 
```
```julia
using FuzzifiED
using ITensors, ITensorMPS, HDF5
using ITensorMPOConstruction
const σx = [  0  1 ;  1  0 ]
const σz = [  1  0 ;  0 -1 ]

function MyMPO(os, sites)
    operatorNames = [ "I", "C", "Cdag", "N" ]
    opCacheVec = [ [OpInfo(ITensors.Op(name, n), sites[n]) for name in operatorNames] for n in eachindex(sites)  ]
    return MPO_new(os, sites ; basis_op_cache_vec = opCacheVec)
end

nm = 12
nf = 2
no = nm * nf

path = "nm_$(nm)_tmp/"
mkpath(path)
```
Like the previous section, we first put in the terms for Hamiltonian and the QNDiags 
```julia
ps_pot = [4.75, 1.] ./ 2
tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, 2, ps_pot) - 
    GetDenIntTerms(nm, 2, ps_pot, σx) - 
    3.16 * GetPolTerms(nm, 2, σz)
)
qnd = [ 
    GetNeQNDiag(no), 
    GetLz2QNDiag(nm, nf), 
    GetZnfChargeQNDiag(nm, nf) 
]
```
The Sites and Hamiltonian MPO can be generated with the function `GetMPOSites`. 
```julia
hmt, sites = GetMPOSites("hmt", tms_hmt, qnd ; path, mpo_method = MyMPO)
```
To generate the initial MPS that respects the $\mathbb{Z}_2$ symmetry, we can use a direct product state. 
```julia
cfi_p = [ [1, 0][f] for m = 1 : nm for f = 1 : nf ]
sti_p = MPS(sites, string.(cfi_p))
cfi_m = [ m == 1 ? [0, 1][f] : [1, 0][f] for m = 1 : nm for f = 1 : nf ]
sti_m = MPS(sites, string.(cfi_m))
```
The lowest eigenenergies and the eigenstate MPSs $|0\rangle,|\sigma\rangle,|\epsilon\rangle$ can be easily generated by the function `EasySweep`.
```julia
E0, st0 = EasySweep("0", hmt, sti_p ; path)
Ee, ste = EasySweep("e", hmt, sti_p ; path, proj = ["0"])
Es, sts = EasySweep("s", hmt, sti_m ; path)
```
To measure the angular momentum $L^2$ of the ground state, we generate the MPO for $L^2$.
```julia
tms_l2 = GetL2Terms(nm, 2)
l2 = GetMPO("l2", tms_l2, sites ; path, mpo_method = MyMPO)
val_l20 = inner(st0', l2, st0)
```
Similarly, to measure the OPE coefficient $f_{\sigma\sigma\epsilon}=\langle \sigma|n^x_{00}|\epsilon\rangle/\langle \sigma|n^x_{00}|0\rangle$
```julia
obs_nx = GetDensityObs(nm, 2, σx)
tms_nx00 = SimplifyTerms(GetComponent(obs_nx, 0.0, 0.0))
nx00 = GetMPO("nx00", tms_nx00, sites ; path, mpo_method = MyMPO)
f_sse = abs(inner(sts', nx00, ste) / inner(sts', nx00, st0))
```

## List of Examples

We offer a series of other examples that reproduces various achievements of fuzzy sphere. For a more detailed summary of the background, see the [review of existing works](@ref Review-of-Existing-Works). 

* [`ising_spectrum.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/ising_spectrum.jl) calculates the spectrum of 3d Ising model on the fuzzy sphere at ``N_m = 12``. For each ``(P,Z,R)`` sector, 20 states are calculated. This example reproduces Table I and Figure 4 in [Zhu 2022](@ref Zhu2022).
* [`ising_phase_diagram.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/ising_phase_diagram.jl) calculates the phase diagram of fuzzy sphere Ising modelby calculating the order parameter ``\langle M^2\rangle``. This example reproduces Figure 3 in [Zhu 2022](@ref Zhu2022).
* [`ising_ope.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/ising_ope.jl) calculates various OPE coefficients at ``N_m = 12`` by taking overlaps between CFT states and density operators and composite. This example reproduces Figure 2 and Table I in [Hu 2023Mar](@ref Hu2023Mar).
* [`ising_correlator.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/ising_correlator.jl) calculates the ``σσ`` two-point function on sphere and the ``σσσσ`` four-point function on sphere, 0 and ``∞``. This example reproduces Figures 1c and 2a in [Han 2023Jun](@ref Han2023Jun).
* [`ising_optimisation.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/ising_optimisation.jl) defines a cost function as the square sum of the deviations of descendants and stress tensor to evaluate the conformal symmetry for Ising model and minimises this cost function to find the best parameter.
* [`ising_full_spectrum.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/ising_full_spectrum.jl) calculates the full spectrum of 3d Ising model on the fuzzy sphere at ``N_m = 10`` for sector ``(P,Z,R) = (1,1,1)``.
* [`ising_space_entangle.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/ising_space_entangle.jl) calculates the entanglement entropy of the Ising ground state along the real space cut of ``θ = 0.500π`` and ``0.499π`` respectively, and use these two data to extract finite size ``F``-function without sustracting the IQHE contribution. This example reproduces Figures 3 in [Hu 2024](@ref Hu2024).
* [`ising_orbital_entangle.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/ising_space_entangle.jl) calculates the entanglement entropy of the Ising ground state along the orbital space cut at ``m = 0``, and also the entanglement spectrum in the half-filled ``l_z = 0, 1`` and  both ``\mathbb{Z}_2`` sectors.
* [`ising_generator.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/ising_generator.jl) examines the quality of conformal symmetry at ``N_m = 12`` by examining the matrix elements of conformal generators ``P^z + K^z`` and compare the states ``(P^z + K^z)|Φ⟩`` with the CFT expectations. This example reproduces Figure 7 in [Fardelli 2024](@ref Fardelli2024).
* [`defect_spectrum.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/defect_spectrum.jl) calculates the spectrum of magnetic line defect in 3d Ising model in ``l_z = 0, P = ±1`` and ``l_z = 1`` sectors, calibrated by bulk ``T``. This example reproduces Table I in [Hu 2023Aug](@ref Hu2023Aug).
* [`defect_correlator.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/defect_correlator.jl) calculates the 1-pt function ``σ`` and 2-pt function ``σ\hat{ϕ}`` of magnetic line defect in 3d Ising model. The normalisation of the correlators require extra bulk data. This example reproduces Figure 4 in [Hu 2023Aug](@ref Hu2023Aug).
* [`defect_changing.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/defect_changing.jl) calculates the spectrum of the defect creation and changing operators of the magnetic line defect in 3d Ising model. This example reproduces Table 2 and Figure 5 in [Zhou 2024Jan](@ref Zhou2024Jan).
* [`defect_overlap.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/defect_overlap.jl) calculates the ``g``-function of magnetic line defect in 3d Ising model using the ovelaps between the bulk, defect ground state and the lowest defect-creation state. This example reproduces Figure 6 in [Zhou 2024Jan](@ref Zhou2024Jan).
* [`cusp_dim.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/cusp_dim.jl) calculates the scaling dimension of the cusp of the magnetic line defect in 3d Ising model as a function of the angle ``θ``. This example reproduces Table 2, upper panel in [Cuomo 2024](@ref Cuomo2024).
* [`surface_ordinary_spectrum.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/surface_ordinary_spectrum.jl) calculates the spectrum of ordinary surface CFT in 3d Ising model calibrated by surface displacement operator ``D`` in the orbital boundary scheme. This example reproduces Figures 3 and 4 in [Zhou 2024Jul](@ref Zhou2024Jul).
* [`surface_normal_spectrum.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/surface_normal_spectrum.jl) calculates the spectrum of normal surface CFT in 3d Ising model calibrated by surface displacement operator ``D`` in the orbital boundary scheme. This example reproduces Figure 5 in [Zhou 2024Jul](@ref Zhou2024Jul).
* [`hsb_2l_spectrum.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/hsb_2l_spectrum.jl) calculates the spectrum of ``\mathrm{O}(3)`` Wilson-Fisher CFT using the bilayer Heisenberg model. This example reproduces Table I and Figure 2 in [Han 2023Dec](@ref Han2023Dec).
* [`so5_spectrum.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/so5_spectrum.jl) calculates the spectrum of SO(5) DQCP on the fuzzy sphere. This example reproduces Table II in [Zhou 2023](@ref Zhou2023).
* [`sp3_spectrum.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/sp3_spectrum.jl) calculates the spectrum of Sp(3) CFT on the fuzzy sphere. This example reproduces Table I in [Zhou 2024Oct](@ref Zhou2024Oct).
* [`ising_frac_fermion.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/ising_frac_fermion.jl) calculates the spectrum of 3d Ising model on the fuzzy sphere for fermions at fractional filling ``ν = 1/3``. This example reproduces Figure 10 in [Voinea 2024](@ref Voinea2024).
* [`ising_frac_boson.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/ising_frac_boson.jl) calculates the spectrum of 3d Ising model on the fuzzy sphere for bosons at fractional filling ``ν = 1/2`` with the module _Fuzzifino_. This example reproduces Figure 12a,b in [Voinea 2024](@ref Voinea2024).
* [`potts_spectrum.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/potts_spectrum.jl) calculates the spectrum of 3d Potts model on the fuzzy sphere. This example reproduces Table I and Figure 4 in [Yang 2025Jan](@ref Yang2025Jan).
* [`circle_ising.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/circle_ising.jl) calculates the spectrum of 2d Ising CFT on a fuzzy thin torus. This example reproduces Figure 4 and Tables I—III in [Han 2025](@ref Han2025).
* [`yang_lee_spectrum.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/yang_lee_spectrum.jl) calculates the spectrum of 3d non-unitary Yang-Lee CFT on the fuzzy sphere. This example reproduces the finite size data of Figure 11 [Arguello Cruz 2025](@ref Fan2025)
* [`free_scalar_spectrum.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/free_scalar_spectrum.jl) calculates the spectrum of free real scalar. This example reproduces Figure 2 in [He 2025Jun](@ref He2025Jun).
* [`free_scalar_correlator.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/free_scalar_correlator.jl) calculates calculates the correlator of $n^i$ for $ϕ$, $ϕ^2$ and $π$ of free real scalar. This example reproduces Figure 3 in [He 2025Jun](@ref He2025Jun)
* [`o4_dqcp_spectrum.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/o4_dqcp_spectrum.jl) calculates the spectrum of $\mathrm{O}(4)$ DQCP on the fuzzy sphere. This example reproduces Figure 4 and Table I in [Yang 2025Jul](@ref Yang2025Jul).
* [`su2_1_scal_spectrum.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/su2_1_scal_spectrum.jl) calculates the spectrum of the $\mathrm{SU}(2)_1$ coupled to a complex scalar Chern-Simons matter CFT on the fuzzy sphere. This example reproduces Figs. 2d and 3 in [Zhou 2025Jul](@ref Zhou2025Jul).
* [`free_majorana_spectrum.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/free_majorana_spectrum.jl) calculates the spectrum of free Majorana fermion. This example reproduces Figure 5 in [Zhou 2025Sep](@ref Zhou2025Sep). 
* [`free_majorana_correlator.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/free_majorana_correlator.jl) calculates the 2-pt correlator of free Majorana fermion. This example reproduces Figure 7 in [Zhou 2025Sep](@ref Zhou2025Sep).
* [`bpf_220_spectrum.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/bpf_220_spectrum.jl) calculates the spectrum of the transition between bosonic Pfaffian and Halperin 220 described by gauged Majorana fermion CFT # This example reproduces Figure 4a in [Voinea 2025](@ref Voinea2025)
* [`o3_wf_spectrum.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/o2_wf_spectrum.jl) calculates the spectrum of $\mathrm{O}(3)$ Wilson-Fisher CFT. This example takes the model from [Dey 2025](@ref Dey2025) and reproduces partly Tables I and II, Figures 1 and 2.
* [`o2_wf_spectrum.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/o2_wf_spectrum.jl) calculates the spectrum of $\mathrm{O}(2)$ Wilson-Fisher CFT. This example reproduces Figure 3, upper panel in [Guo 2025](@ref Guo2025).
* [`o2_free_spectrum.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/o2_free_spectrum.jl) calculates the spectrum of $\mathrm{O}(2)$ free scalar CFT. This example reproduces Figure 4, upper-left panel in [Guo 2025](@ref Guo2025).
* [`o4_wf_spectrum.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/o4_wf_spectrum.jl) calculates the spectrum of $\mathrm{O}(4)$ Wilson-Fisher CFT. This example reproduces Figure 3, lower panel in [Guo 2025](@ref Guo2025).
* [`ising_spectrum_krylov.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/ising_spectrum_krylov.jl) calculates the spectrum of 3d Ising model on the fuzzy sphere by calling the eigsolve function in KrylovKit.jl instead of Arpack.
* [`ising_spectrum_cuda.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/ising_spectrum_cuda.jl) calculates the spectrum of 3d Ising model on the fuzzy sphere for one sector by performing the sparse matrix multiplication on CUDA.