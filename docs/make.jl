using CytoZoo
using Documenter

DocMeta.setdocmeta!(CytoZoo, :DocTestSetup, :(using CytoZoo); recursive=true)

makedocs(;
    modules=[CytoZoo],
    authors="Kyle Beggs (beggskw@gmail.com) and contributors",
    sitename="CytoZoo.jl",
    format=Documenter.HTML(;
        canonical="https://kylebeggs.github.io/CytoZoo.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kylebeggs/CytoZoo.jl",
    devbranch="main",
)
