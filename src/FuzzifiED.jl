module FuzzifiED 

using LinearAlgebra
using SparseArrays
using Requires
using WignerSymbols
using SphericalHarmonics
using FuzzifiED_jll
using HDF5
import Base.:+
import Base.:-
import Base.:*
import Base.:/
import Base.:÷
import Base.:^
import Base.zero
import Base.adjoint

include("core/param.jl")
export NumThreads
export SilentStd
export Libpath
export ElementType

include("core/qn.jl")
export QNDiag
export QNOffd

include("core/confs.jl")
export Confs
export GetConfId

include("core/basis.jl")
export Basis
export GetConfWeight

include("core/term.jl")
export Term
export ParticleHole
export NormalOrder
export SimplifyTerms
export SimplifyTermsOld

include("core/operator.jl")
export Operator

include("core/opmat.jl")
export OpMat
export GetEigensystem
export SparseMatrixCSCFromOpMat
export MatrixFromOpMat

include("core/entangle.jl")
export StateDecompMat
export GetEntSpec

include("models/qndiag.jl")
export GetNeQNDiag
export GetLz2QNDiag
export GetFlavQNDiag
export GetZnfChargeQNDiag
export GetPinOrbQNDiag

include("models/qnoffd.jl")
export GetParityQNOffd
export GetFlavPermQNOffd
export GetRotyQNOffd

include("models/opterms.jl")
export GetIntMatrix
export GetDenIntTerms
export GetPairIntTerms
export GetPolTerms
export GetIsingIntTerms
export GetL2Terms
export GetC2Terms

include("models/sphere_obs.jl")
export SphereObs
export StoreComps!
export StoreComps
export Laplacian
export GetComponent
export GetPointValue
export Electron
export Density
export Pairing
export PairObs

include("archieve/ar_core.jl")

include("archieve/ar_models.jl")
export GetLzQnu
export GetLzZnQnu
export GetLzConfs 
export GetLzZnConfs
export GetIsingQnz
export GetIsingBasis
export GetSnBasis
export GetXPolTerms
export GetZPolTerms

function __init__()
    FuzzifiED.NumThreads = Threads.nthreads()

    @require ITensors = "9136182c-28ba-11e9-034c-db9fb085ebd5" begin
        using ITensors

        include("itensors_support/itensors_format.jl")
        export QNDiagFromSites
        export ConfsFromSites
        export TermsFromOpSum
        export OpSumFromTerms
        export SitesFromQNDiag
        export TruncateQNDiag

        include("itensors_support/easy_sweep.jl")
        export SweepOne
        export EasySweep
        export GetMPOSites
        export GetMPO

        include("archieve/ar_itensor.jl")
        export TruncateQnu
        export SitesFromQnu
    end
end

end