module FuzzifiEDITensorsExt

using ITensors 
using ITensorMPS
using FuzzifiED
using LinearAlgebra

import Base.Vector
import FuzzifiED.QNDiagFromSites
import FuzzifiED.ConfsFromSites
import FuzzifiED.TermsFromOpSum
import FuzzifiED.OpSumFromTerms
import FuzzifiED.SitesFromQNDiag
import FuzzifiED.GetSites
import FuzzifiED.TruncateQNDiag
import FuzzifiED.SilentStd

include("fermion_type.jl")
include("itensors_format.jl")

function __init__()
    BLAS.set_num_threads(1);
    NDTensors.Strided.disable_threads();
    ITensors.enable_threaded_blocksparse();
end

end 