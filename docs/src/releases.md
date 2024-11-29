# Release notes 

## Version 0.x

### Version 0.10

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

### Version 0.7.0 (9 June, 2024)

- Revise the implementation of diagonal and off-diagonal quantum number. 

### Version 0.6

### Version 0.6.3 (8 June, 2024)

- Add support for calculating entanglement spectrum. 
- Add global parameters to control the number of threads, the output and the path of the dynamic library. 
- Fix bugs. 

### Version 0.6.0 (5 June, 2024)

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