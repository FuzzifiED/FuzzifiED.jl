module ITensorsExt

using ITensors 
using ITensorMPS
using FuzzifiED
using HDF5
using LinearAlgebra

import FuzzifiED.QNDiagFromSites
import FuzzifiED.ConfsFromSites
import FuzzifiED.TermsFromOpSum
import FuzzifiED.OpSumFromTerms
import FuzzifiED.SitesFromQNDiag
import FuzzifiED.TruncateQNDiag
import FuzzifiED.SweepOne
import FuzzifiED.EasySweep
import FuzzifiED.GetMPOSites
import FuzzifiED.GetMPO
import FuzzifiED.TruncateQnu
import FuzzifiED.SitesFromQnu

include("itensors_format.jl")
include("easy_sweep.jl")
include("ar_itensor.jl")

function __init__()
    # Define the space method at runtime
    @eval ITensors.space( :: SiteType"Fermion"; o :: Int, qnd :: Vector{QNDiag}) = [
        QN(
            [ (qndi.name, qndi.charge[o] * n, qndi.modul) for qndi in qnd ]...
        ) => 1 for n = 0 : 1
    ]
    BLAS.set_num_threads(1);
    NDTensors.Strided.disable_threads();
    ITensors.enable_threaded_blocksparse();
end

end 