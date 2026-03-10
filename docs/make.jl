using CrystalUtils
using Documenter

DocMeta.setdocmeta!(CrystalUtils, :DocTestSetup, :(using CrystalUtils); recursive = true)

makedocs(;
    modules = [CrystalUtils],
    authors = "Junfeng Qiao <qiaojunfeng@outlook.com> and contributors",
    sitename = "CrystalUtils.jl",
    pages = [
        "Home" => "index.md",
        "API" => [
            "Type" => "api/type.md",
            "Lattice" => "api/lattice.md",
            "Kpoint Path" => "api/kpath.md",
        ],
    ],
)

deploydocs(; repo = "github.com/atomology/CrystalUtils.jl", devbranch = "main")
