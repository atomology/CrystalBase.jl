using CrystalBase
using Documenter

DocMeta.setdocmeta!(CrystalBase, :DocTestSetup, :(using CrystalBase); recursive = true)

makedocs(;
    modules = [CrystalBase],
    authors = "Junfeng Qiao <qiaojunfeng@outlook.com> and contributors",
    sitename = "CrystalBase.jl",
    pages = [
        "Home" => "index.md",
        "API" => [
            "Type" => "api/type.md",
            "Lattice" => "api/lattice.md",
            "Atom" => "api/atom.md",
            "Kpoint Path" => "api/kpath.md",
        ],
    ],
)

deploydocs(; repo = "github.com/atomology/CrystalBase.jl", devbranch = "main")
