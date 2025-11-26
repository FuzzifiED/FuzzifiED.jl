#=
julia -O3 --color=yes make.jl && rm -r publish/*/* && mv -f build/* publish && rm -r build && cd publish && git commit -a -m "a" && git push && cd ..
=#
push!(LOAD_PATH,"../src/")
push!(LOAD_PATH,"../ext/")

using Documenter
using WignerSymbols
using SparseArrays
using ITensors
using ITensorMPS
using CUDA
using LinearAlgebra
using HDF5
using KrylovKit
using FuzzifiED
using FuzzifiED.Fuzzifino
using FuzzifiED.FuzzyManifolds

makedocs(sitename = "FuzzifiED.jl", 
    pages = ["Home" => "index.md", 
        "Introduction" => "intro.md",
        "Tutorial" => "tutorial.md",
        "Core Functions" => "core.md",
        "Built-in Models" => "models.md",
        "ITensor Extension" => "itensors.md",
        "Other Extensions" => "extension.md", 
        "Fuzzifino" => "fuzzifino.md",
        "Fuzzy Manifolds" => "manifolds.md",
        "Releases" => "releases.md"],
    format = Documenter.HTML(
        assets = ["assets/serif.css", "assets/favicon.ico"], 
        repolink = "https://github.com/FuzzifiED/FuzzifiED.jl",
        footer = "Powered by [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) and the [Julia Programming Language](https://julialang.org/). Copyright (c) 2025 Zheng Zhou (周正) and contributors."
    )
)
