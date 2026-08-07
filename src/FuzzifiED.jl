module FuzzifiED 

using LinearAlgebra
using WignerSymbols
using SphericalHarmonics
using FuzzifiED_jll

include("core/param.jl")
include("core/qn.jl")
include("core/confs.jl")
include("core/basis.jl")
include("core/term.jl")
include("core/operator.jl")
include("core/opmat.jl")
include("core/transf.jl")
include("core/entangle.jl")

include("models/qndiag.jl")
include("models/qnoffd.jl")
include("models/opterms.jl")
include("models/sphere_obs.jl")
include("models/ang_modes.jl")

include("fuzzifino/fuzzifino.jl")
include("manifolds/manifolds.jl")
include("jack/jack.jl")

export GetSites, GetQNDiags, TruncateQNDiag
function GetSites end
function GetQNDiags end
function TruncateQNDiag end

export SweepOne, EasySweep, GetMPOSites, GetMPO
function SweepOne end
function EasySweep end
function GetMPOSites end
function GetMPO end

export GetEigensystemKrylov
function GetEigensystemKrylov end

export GetEigensystemCuda
function GetEigensystemCuda end

function __init__()
    FuzzifiED.NumThreads = Threads.nthreads()
end

end
