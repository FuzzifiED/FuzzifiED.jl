# FuzzifiED.jl 

[![Version](https://img.shields.io/badge/Version-1.1.0-orange)](https://juliahub.com/ui/Packages/General/FuzzifiED/1.1.0)
[![Documentation online](https://img.shields.io/badge/Documentation-Online-8e8eff)](https://docs.fuzzified.world/)
[![Documentation PDF](https://img.shields.io/badge/Documentation-PDF-8e8eff)](https://docs.fuzzified.world/assets/FuzzifiED_Documentation.pdf)
[![Source GitHub](https://img.shields.io/badge/Source-GitHub-silver)](https://github.com/FuzzifiED/FuzzifiED.jl)
[![arXiv](https://img.shields.io/badge/arXiv-2503.00100-b31b1b)](https://arxiv.org/abs/2503.00100)
[![Contact](https://img.shields.io/badge/Contact-Zheng_Zhou_周正-2e63b8)](mailto:Zheng%20Zhou%20%E5%91%A8%E6%AD%A3%20<physics@zhengzhou.page>)

Since its proposal, the fuzzy sphere regularisation has made significant contributions to the study of 3d CFTs. The Julia package FuzzifiED aims at simplifying the numerical calculations on the fuzzy sphere. It supports exact diagonalisation (ED) calculations, as well as the density matrix renormalisation group (DMRG) using the ITensor library. FuzzifiED can also apply to generic fermionic and bosonic models. This package offers the following features : 

* Flexibility : FuzzifiED can reproduce nearly all the ED and DMRG results from fuzzy sphere research. Its flexible design also makes it straightforward to adapt to new models.
* Usability : The expressive Julia interface simplifies coding and comprehension. To help users get started, we provide [a collection of examples](@ref List-of-examples).
* Efficiency : FuzzifiED produces results on reasonable system sizes within minutes.
* Open source : The FuzzifiED codebase is freely available under the MIT License, welcoming reviews and contributions from the wider community.

A PDF version of the documentation is provided at [this link](https://docs.fuzzified.world/assets/FuzzifiED_Documentation.pdf). If you have any questions, please contact Zheng Zhou (周正) at [physics@zhengzhou.page](mailto:Zheng%20Zhou%20%E5%91%A8%E6%AD%A3%20<physics@zhengzhou.page>).

## Installation

To install the package, run the following command in the Julia REPL (read-eval-print loop) (To enter Julia REPL, simply type `julia` in the command line) 
```julia
using Pkg ; Pkg.add("FuzzifiED")
```
To use the package, include at the start of the Julia script
```julia
using FuzzifiED
```
To obtain the documentation for an interface, type `?` followed by the keyword in the Julia REPL, _e.g._, `?Confs`.

## Citation

If this package is helpful in your research, please cite the package as : 

> FuzzifiED : Julia package for numerics on the fuzzy sphere, Zheng Zhou, [arXiv:2503.00100](https://arxiv.org/abs/2503.00100).

We have also provided a BibTeX file that includes all the works on the fuzzy sphere works at [this link](https://docs.fuzzified.world/assets/bib_fuzzy.bib).

## Useful information

* Jupyter Notebook is highly recommended as it allows you to run Julia (and Python) just like running a Mathematica notebook.
* The package regisitry may have some delay. If you encounter trouble at installation, to bring the registry up to date, use `Pkg.Registry.update()`, or install from the GitHub repos.
```Julia
using Pkg
Pkg.add(url="https://github.com/FuzzifiED/FuzzifiED_jll.jl")
Pkg.add(url="https://github.com/FuzzifiED/FuzzifiED.jl")
```
* Please find following some useful links.
    - Installation for Julia : <https://julialang.org/downloads>
    - Homepage : <https://www.fuzzified.world>
    - Julia source code : <https://github.com/FuzzifiED/FuzzifiED.jl>
    - JLL wrapper : <https://github.com/FuzzifiED/FuzzifiED_jll.jl>
    - Fortran source code : <https://github.com/FuzzifiED/FuzzifiED_Fortran>
    - Registry of the package : <https://juliahub.com/ui/Packages/General/FuzzifiED>

## Outline 

```@contents
Pages = [
    "intro.md",
    "tutorial.md",
    "core.md",
    "models.md",
    "itensors.md",
    "extension.md",
    "fuzzifino.md",
    "manifolds.md",
    "releases.md"
]
Depth = 2
```

## References

* __[Zhu 2022]__ Uncovering conformal symmetry in the 3d Ising transition : state-operator correspondence from a quantum fuzzy sphere regularisation, Wei Zhu, Chao Han, Emilie Huffman, Johannes S. Hofmann, and Yin-Chen He, [arXiv:2210.13482](https://arxiv.org/abs/2210.13482), [Phys. Rev. X __13__, 021009 (2023)](https://doi.org/10.1103/PhysRevX.13.021009).
* __[Hu 2023Mar]__ Operator product expansion coefficients of the 3d Ising criticality via quantum fuzzy sphere, Liangdong Hu, Yin-Chen He, and Wei Zhu, [arXiv:2303.08844](https://arxiv.org/abs/2303.08844), [Phys. Rev. Lett __131__, 031601 (2023)](https://doi.org/10.1103/PhysRevLett.131.031601).
* __[Han 2023Jun]__ Conformal four-point correlators of the 3d Ising transition via the quantum fuzzy sphere, Chao Han, Liangdong Hu, Wei Zhu, and Yin-Chen He, [arXiv:2306.04681](https://arxiv.org/abs/2306.04681), [Phys. Rev. B __108__, 235123 (2023)](https://doi.org/10.1103/PhysRevB.108.235123).
* __[Zhou 2023]__ The ``\mathrm{SO}(5)`` deconfined phase transition under the fuzzy sphere microscope: approximate conformal symmetry, pseudo-criticality, and operator spectrum, Zheng Zhou, Liangdong Hu, Wei Zhu, and Yin-Chen He, [arXiv:2306.16435](https://arxiv.org/abs/2306.16435), [Phys. Rev. X __14__, 021044 (2024)](https://doi.org/10.1103/PhysRevX.14.021044).
* __[Hu 2023Aug]__ Solving conformal defects in 3d conformal field theory using fuzzy sphere regularisation, Liangdong Hu, Yin-Chen He, and Wei Zhu, [arXiv:2308.01903](https://arxiv.org/abs/2308.01903), [Nat. Commun. __15__, 3659 (2024)](https://doi.org/10.1038/s41467-024-47978-y).
* __[Hofmann 2024]__ Quantum Monte Carlo simulation of the 3d Ising transition on the fuzzy sphere, Johannes S. Hofmann, Florian Goth, Wei Zhu, Yin-Chen He, and Emilie Huffman, [arXiv:2310.19880](https://arxiv.org/abs/2310.19880), [SciPost Phys. Core __7__, 028 (2024)](https://doi.org/10.21468/SciPostPhysCore.7.2.028).
* __[Han 2023Dec]__ Conformal operator content of the Wilson-Fisher transition on fuzzy sphere bilayers, Chao Han, Liangdong Hu, and Wei Zhu, [arXiv:2312.04047](https://arxiv.org/abs/2312.04047), [Phys. Rev. B __110__, 115113 (2024)](https://doi.org/10.1103/PhysRevB.110.115113).
* __[Zhou 2024Jan]__ The ``g``-function and defect changing operators from wavefunction overlap on a fuzzy sphere, Zheng Zhou, Davide Gaiotto, Yin-Chen He, Yijian Zou, [arXiv:2401.00039](https://arxiv.org/abs/2401.00039), [SciPost Phys. __17__, 021 (2024)](https://doi.org/10.21468/SciPostPhys.17.1.021).
* __[Hu 2024]__ Entropic ``F``-function of 3d Ising conformal field theory via the fuzzy sphere regularisation, Liangdong Hu, Wei Zhu, and Yin-Chen He, [arXiv:2401.17362](https://arxiv.org/abs/2401.17362), [Phys. Rev. B __111__, 155151 (2025)](https://doi.org/10.1103/PhysRevB.111.155151).
* __[Cuomo 2024]__ Impurities with a cusp : general theory and 3d Ising, Gabriel Cuomo, Yin-Chen He, Zohar Komargodski, [arXiv:2406.10186](https://arxiv.org/abs/2406.10186), [JHEP __11__ (2024) 061](https://doi.org/10.1007/JHEP11(2024)061). 
* __[Zhou 2024Jul]__ Studying the 3d Ising surface CFTs on the fuzzy sphere, Zheng Zhou, and Yijian Zou, [arXiv:2407.15914](https://arxiv.org/abs/2407.15914), [SciPost Phys. __18__, 031 (2025)](https://doi.org/10.21468/SciPostPhys.18.1.031).
* __[Dedushenko 2024]__ Ising BCFTs from the fuzzy hemisphere, Mykola Dedushenko, [arXiv:2407.15948](https://arxiv.org/abs/2407.15948).
* __[Fardelli 2024]__ Constructing the infrared conformal generators on the fuzzy sphere, Giulia Fardelli, A. Liam Fitzpatrick, and Emanuel Katz, [arXiv:2409.02998](https://arxiv.org/abs/2409.02998), [SciPost Phys. __18__, 086 (2025)](https://doi.org/10.21468/SciPostPhys.18.3.086).
* __[Fan 2024]__ Note on explicit construction of conformal generators on the fuzzy sphere, Ruihua Fan, [arXiv:2409.08257](https://arxiv.org/abs/2409.08257).
* __[Zhou 2024Oct]__ 3D conformal field theories with ``\mathrm{Sp}(N)`` global symmetry on fuzzy sphere, Zheng Zhou, and Yin-Chen He, [arXiv:2410.00087](https://arxiv.org/abs/2410.00087), [Phys. Rev. Lett. __135__, 026504 (2025)](https://doi.org/10.1103/xstj-xvcy).
* __[Voinea 2024]__ Regularising 3d conformal field theories via anyons on the fuzzy sphere, [arXiv:2411.15299](https://arxiv.org/abs/2411.15299), [Phys. Rev. X __15__, 031007 (2025)](https://doi.org/10.1103/bf4k-phl9).
* __[Han 2025]__ Quantum phase transitions on the noncommutative circle, Chao Han, and Wei Zhu, [Phys. Rev. B __111__, 085113 (2025)](https://doi.org/10.1103/PhysRevB.111.085113).
* __[Yang 2025Jan]__ Microscopic study of 3d Potts phase transition via fuzzy sphere regularisation, Shuai Yang, Yan-Guang Yue, Yin Tang, Chao Han, Wei Zhu, and Yan Chen, [arXiv:2501.14320](https://arxiv.org/abs/2501.14320), [Phys. Rev. B __112__, 024436 (2025)](https://doi.org/10.1103/x1qn-x6xb).
* __[Läuchli 2025]__ Exact diagonalization, matrix product states and conformal perturbation theory study of a 3d Ising fuzzy sphere model, Andreas M. Läuchli, Loïc Herviou, Patrick H. Wilhelm, and Slava Rychkov, [arXiv:2504.00842](https://arxiv.org/abs/2504.00842).
* __[Fan 2025]__ Simulating the non-unitary Yang-Lee conformal field theory on the fuzzy sphere, Ruihua Fan, Junkai Dong, and Ashvin Vishwanath [arXiv:2505.06342](https://arxiv.org/abs/2505.06342).
* __[Arguello Cruz 2025]__ Yang-Lee quantum criticality in various dimensions, Erick Arguello Cruz, Igor R. Klebanov, Grigory Tarnopolsky, and Yuan Xin, [arXiv:2505.06369](https://arxiv.org/abs/2505.06369).
* __[Elias Miro 2025]__ Flowing from the Ising model on the fuzzy sphere to the 3d Lee-Yang CFT, Joan Elias Miro, Olivier Delouche, [arXiv:2505.07655](https://arxiv.org/abs/2505.07655).
* __[He 2025Jun]__ Free real scalar CFT on fuzzy sphere : Spectrum, algebra and wavefunction ansatz, [arXiv:2506.14904](https://arxiv.org/abs/2506.14904).
* __[Taylor 2025]__ Conformal scalar field theory from Ising tricriticality on the fuzzy sphere, Joseph Taylor, Cristian Voinea, Zlatko Papić, Ruihua Fan, [arXiv:2506.22539](https://arxiv.org/abs/2506.22539).
* __[Yang 2025Jul]__ Conformal Operator Flows of the Deconfined Quantum Criticality from $\mathrm{SO}(5)$ to $\mathrm{O}(4)$, Shuai Yang, Liang-Dong Hu, Chao Han, Wei Zhu, Yan Chen, [arXiv:2507.01322](https://arxiv.org/abs/2507.01322).
* __[Zhou 2025Jul]__ Chern-Simons-matter conformal field theory on fuzzy sphere: Confinement transition of Kalmeyer-Laughlin chiral spin liquid, Zheng Zhou, Chong Wang, Yin-Chen He, [arXiv:2507.19580](https://arxiv.org/abs/2507.19580).

## Index 

```@index
Pages = [
    "core.md",
    "models.md",
    "itensors.md",
    "extension.md",
    "fuzzifino.md",
    "manifolds.md"
]
```
