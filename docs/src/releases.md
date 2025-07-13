# Release notes 

## Version 1.x

### Version 1.1

#### Version 1.1.0 (14 July, 2025)

- Add the removal and relabelling of orbitals for Terms. 
- Include SphereObs and AngModes for Fuzzifino.
- Add Contract for AngModes. 
- Fix the normalisation of Casimir normalisation.
- Update the reference and examples. 

### Version 1.0

#### Version 1.0.4 (15 May, 2025)

- Add example of Yang-Lee CFT.
- Update compat for ITensor.

#### Version 1.0.3 (30 April, 2025)

- Relax constraint on negative QNDiag and particle conservation. 
- Add GetTorusTranslQNOffd for torus. (We acknowledge Shuai Yang for the suggestion.)

#### Version 1.0.2 (28 March, 2025)

- Enable the normalisation of the observables. 
- Add uniform spatial integral and component filter of the observables. 
- Add conversion from Term to STerm.

#### Version 1.0.1 (11 March, 2025)

- Implement the pure bosonic models.

#### Version 1.0.0 (1 March, 2025)

- Release together with the documentation.
- Add support for torus. 
- Add example for Potts. 
- Fix bugs.

## Version 0.x

### Version 0.10

#### Version 0.10.7 (9 December, 2024)

- Add implementation of transformations.

#### Version 0.10.6 (6 December, 2024)

- Add entanglement calculations in Fuzzifino. 

#### Version 0.10.5 (4 December, 2024)

- Fix the dependency in FuzzifiED_jll
- Move KrylovKit and SparseArrays to extensions. 
- Fix minor bugs. 

#### Version 0.10.4 (30 November, 2024)

- Split ITensorsExt and EasySweepExt.
- Define a new SiteType "FuzzyFermion" for ITensors to avoid overwriting "Fermion" in ITensorsExt.
- Change the names of Density and Electron to avoid conflict with other packages.
- Create the alias Terms and STerms.
- Fix bugs in converting OpSum to Terms.
- Improve the interface for ITensor extension.

#### Version 0.10.3 (28 November, 2024)

- Add options for DMRG. 
- Optimise the implementation of confs in Fuzzifino.
- Add the examples of fractional filling.

#### Version 0.10.1 (25 November, 2024)

- Add CUDA extension. 
- Add interface with KrylovKit.

#### Version 0.10.0 (24 November, 2024)

- Add module Fuzzifino for boson-fermion systems.

### Version 0.9

#### Version 0.9.3 (23 November, 2024)

- Move HDF5 and ITensor to extensions
- Modify ITensor extension interfaces in alignment with the update of ITensor.

#### Version 0.9.2 (16 September, 2024)

- Add the example of Ising generators and Sp(3) CFT.

#### Version 0.9.1 (13 September, 2024)

- Allow input of initial vectors for diagonalisation. (We acknowledge Andrew Fitzpatrick for the suggestion.)

#### Version 0.9.0 (11 September, 2024)

- Add feature of angular modes observables.
- Fix typos and bugs.

### Version 0.8

#### Version 0.8.2 (28 July, 2024)

- Fix bugs in GetDenIntTerms and multiplication in SphereObs. 

#### Version 0.8.0 (26 July, 2024)

- Improve the performance of SimplifyTerms. 
- Add file operation. 

### Version 0.7

#### Version 0.7.2 (24 July, 2024)

- Add example and support for surface CFTs. 
- Add example of Ising cusp.

#### Version 0.7.1 (11 June, 2024)

- Add some new interfaces for built-in operators. 
- Add new examples. 
- Fix bugs

#### Version 0.7.0 (9 June, 2024)

- Revise the implementation of diagonal and off-diagonal quantum number. 

### Version 0.6

#### Version 0.6.3 (8 June, 2024)

- Add support for calculating entanglement spectrum. 
- Add global parameters to control the number of threads, the output and the path of the dynamic library. 
- Fix bugs. 

#### Version 0.6.0 (5 June, 2024)

- Add support for full diagonalisation. 
- Fix bugs and typos.
- Ready for formal release !

### Version 0.5

#### Version 0.5.8 (3 June, 2024)

- Change the binary dependence to Julia Binary Builder. 

#### Version 0.5.0 (30 May, 2024)

- Enable simplification of terms.
- Add general observables and built-in electrons and density operators. 
- Reorganise the realisations of built-in models.
- Cancel ITensorMPOConstruction dependence. 

### Version 0.4

#### Version 0.4.3 (29 May, 2024)

- Add QNU truncation for ITensor use.
- Change the Fortran code to be robust against QNU breaking terms
- Add built-in density-density interaction. 
- Add built-in 3-state Potts model.
- Add built-in Ising model with magnetic line defect. 

#### Version 0.4.0 (28 May, 2024)

- Add support for DMRG.
- Add convertion from diagonal QNs to sites. 
- Add the support of ``\mathbb{Z}_n`` diagonal quantum numbers in `Confs`.
- Merge the submodules to the main package. 
- Add Ising model in X basis.  
- Add functions in built-in models to export diagonal QNs. 

### Version 0.3

#### Version 0.3.0 (27 May, 2024)

- Add the support for Hamiltonians with real elements. 
- Add conversion with SparseMatrixCSC. 
- Add the conversion from terms to OpSum.
- Add the look-up of configurations. 
- For the built-in Ising model, add density operator.
- Add built-in ``\mathrm{Sp}(N)`` model. 

### Version 0.2

#### Version 0.2.0 (26 May, 2024)

- Add operations of terms.
- Add built-in Ising model. 