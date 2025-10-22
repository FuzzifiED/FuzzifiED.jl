
"""
    FuzzifiED.SilentStd :: Bool = false 

a flag to determine whether logs of the FuzzifiED functions should be turned off. False by default. If you want to evaluate without log, put `FuzzifiED.SilentStd = true`. This parameter can be defined for each process separately. 
"""
SilentStd :: Bool = false

"""
    FuzzifiED.Libpath :: String = FuzzifiED_jll.LibpathFuzzifiED

define path of the Fortran library `libfuzzified.so`. You do not need to modify that by yourself. However, if you compile the Fortran codes by yourself, you need to point this to your compiled library. 
"""
Libpath :: String = FuzzifiED_jll.LibpathFuzzifiED

"""
    FuzzifiED.NumThreads :: Int = Threads.nthreads()

an integer to define how many threads OpenMP uses. By default, it is the same as the number of threads in Julia. If you use Jupyter notebooks, which by default uses one core only, you may need to define this by hand, _e. g._, `FuzzifiED.NumThreads = 8`. This parameter can be defined for each process separately. 
"""
NumThreads :: Int = Threads.nthreads()

"""
    FuzzifiED.ElementType :: DataType = ComplexF64

set the default type of the operator elements, either `ComplexF64` or `Float64`. `ComplexF64` by default. 
"""
ElementType :: Union{Type{ComplexF64}, Type{Float64}} = ComplexF64

"""
    FuzzifiED.OpenHelp!()

A shortcut to open the link for documentation [docs.fuzzified.world](https://docs.fuzzified.world) in the system browser. 
"""
function OpenHelp!()
    url = "https://docs.fuzzified.world"
    if Sys.iswindows()
        run(`start $url` ; wait = false)
    elseif Sys.isapple()
        run(`open $url` ; wait = false)   # macOS
    elseif Sys.isunix()
        run(`xdg-open $url` ; wait = false)  # Linux
    end
end