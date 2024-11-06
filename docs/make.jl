using PiccoloQuantumObjects
using Documenter

DocMeta.setdocmeta!(PiccoloQuantumObjects, :DocTestSetup, :(using PiccoloQuantumObjects); recursive=true)

makedocs(;
    modules=[PiccoloQuantumObjects],
    authors="Aaron Trowbridge <aaron.j.trowbridge@gmail.com> and contributors",
    sitename="PiccoloQuantumObjects.jl",
    format=Documenter.HTML(;
        canonical="https://aarontrowbridge.github.io/PiccoloQuantumObjects.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/aarontrowbridge/PiccoloQuantumObjects.jl",
    devbranch="main",
)
