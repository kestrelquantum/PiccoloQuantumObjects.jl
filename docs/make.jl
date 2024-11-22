using PiccoloQuantumObjects
using Documenter
using Literate

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

pages = [
    "Home" => "index.md",
    "Library" => "lib.md",
]

format = Documenter.HTML(;
    prettyurls=get(ENV, "CI", "false") == "true",
    canonical="https://kestrelquantum.github.io/PiccoloQuantumObjects.jl",
    edit_link="main",
    assets=String[],
    mathengine = MathJax3(Dict(
        :loader => Dict("load" => ["[tex]/physics"]),
        :tex => Dict(
            "inlineMath" => [["\$","\$"], ["\\(","\\)"]],
            "tags" => "ams",
            "packages" => [
                "base",
                "ams",
                "autoload",
                "physics"
            ],
        ),
    )),
)

src = joinpath(@__DIR__, "src")
lit = joinpath(@__DIR__, "literate")

lit_output = joinpath(src, "generated")

for (root, _, files) ∈ walkdir(lit), file ∈ files
    splitext(file)[2] == ".jl" || continue
    ipath = joinpath(root, file)
    opath = splitdir(replace(ipath, lit=>lit_output))[1]
    Literate.markdown(ipath, opath)
end

makedocs(;
    modules=[PiccoloQuantumObjects],
    authors="Aaron Trowbridge <aaron.j.trowbridge@gmail.com> and contributors",
    sitename="PiccoloQuantumObjects.jl",
    format=format,
    pages=pages,
    warnonly=true,
)

deploydocs(;
    repo="github.com/kestrelquantum/PiccoloQuantumObjects.jl.git",
    devbranch="main",
)
