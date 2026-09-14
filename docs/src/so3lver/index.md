# SO(3)lver 

SO(3)lver is an extension of the package [FuzzifiED](https://docs.fuzzified.world). It enables exact diagonalization (ED) on the fuzzy sphere with the implementation of full $\mathrm{SO}(3)$ spatial rotation symmetry and some non-Abelian flavour symmetry through partitioning the Hilbert space into segments and compose these segments. 

The implementation of $\mathrm{SO}(3)$ symmetry greatly reduces the demand for the memory and enables calculation for much larger size. _E. g._, with the $\mathrm{SO}(3)$ resolution, the Hilbert space of the ground-state sector for $N_m=22$ has a dimension $2.07\times 10^7$, and the Hamiltonian requires a memory of 76.5 Go ; as a comparison, without the $\mathrm{SO}(3)$ resolution, the Hilbert space for $N_m=18$ has a similar dimension $2.84\times 10^7$, and the Hamiltonian requires a memory of 63.9 Go.

SO(3)lver can apply to a lot of current fuzzy-sphere research. Its flexible design also makes it straightforward to adapt to new models. We provide [a collection of examples](@ref Examples-Using-SO(3)lver).

To use the package, include at the start of the Julia script
```julia
using FuzzifiED
using FuzzifiED.SO3lver
```

If this module is helpful in your research, please cite : 

> Implementing Full $\mathrm{SO}(3)$ Rotation Symmetry for Exact Diagonalization on the Fuzzy Sphere, Zheng Zhou and Yin-Chen He, arXiv:2509.xxxxx (_to appear_).

## Outline 

```@contents
Pages = [
    "formalism.md",
    "tutorial.md",
    "interface.md"
]
Depth = 2
```

## Index 

```@index
Pages = [
    "interface.md"
]
```
