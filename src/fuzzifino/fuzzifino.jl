export Fuzzifino

module Fuzzifino

using FuzzifiED_jll
using FuzzifiED
using LinearAlgebra

import FuzzifiED: NumThreads, SilentStd, ElementType

"""
    Fuzzifino.Libpathino :: String = FuzzifiED_jll.LibpathFuzzifino

define path of the Fortran library `libfuzzifino.so`. You do not need to modify that by yourself. However, if you compile the Fortran codes by yourself, you need to point this to your compiled library. 
"""
Libpathino :: String = FuzzifiED_jll.LibpathFuzzifino

include("core/sqn.jl")
include("core/sconfs.jl")
include("core/sbasis.jl")
include("core/sterm.jl")
include("core/soperator.jl")
include("core/sopmat.jl")
include("core/stransf.jl")
include("core/sentangle.jl")

include("model/boson_sqn.jl")
include("model/boson_opsterms.jl")

end
