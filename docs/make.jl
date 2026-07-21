using CytoZoo
using Documenter
using DocumenterVitepress

DocMeta.setdocmeta!(CytoZoo, :DocTestSetup, :(using CytoZoo); recursive = true)

makedocs(;
    modules = [CytoZoo],
    authors = "Kyle Beggs (beggskw@gmail.com) and contributors",
    sitename = "CytoZoo.jl",
    repo = Documenter.Remotes.GitHub("DerangedIons", "CytoZoo.jl"),
    format = DocumenterVitepress.MarkdownVitepress(;
        repo = "https://github.com/DerangedIons/CytoZoo.jl",
        devbranch = "main",
        devurl = "dev",
        build_vitepress = (!haskey(ENV, "VITEPRESS_DEV")),
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Guides" => [
            "The Cell Model Interface" => "guides/interface.md",
            "Stimulus" => "guides/stimulus.md",
            "Spatial Heterogeneity" => "guides/spatial.md",
            "Rush-Larsen Integration" => "guides/rush_larsen.md",
            "Coupling" => [
                "Overview" => "guides/coupling/index.md",
                "Share Edges" => "guides/coupling/share.md",
                "Connect Edges" => "guides/coupling/connect.md",
                "Patterns Cookbook" => "guides/coupling/patterns.md",
                "Limitations" => "guides/coupling/limitations.md",
            ],
            "Derived Observables" => "guides/monitors.md",
            "Implementing a Model" => "guides/implementing_a_model.md",
            "Integrations" => "guides/integrations.md",
            "Quick Reference" => "guides/quickref.md",
        ],
        "Reference" => [
            "Model Catalog" => "reference/models.md",
            "Design Notes" => "reference/design.md",
            "Coupling Internals" => "reference/internals.md",
            "API" => "reference/api.md",
        ],
    ],
    checkdocs = :exports,
    clean = false,
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/DerangedIons/CytoZoo.jl",
    target = joinpath(@__DIR__, "build"),
    branch = "gh-pages",
    devbranch = "main",
    push_preview = true,
)
