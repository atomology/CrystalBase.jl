using CrystalUtils
using Documenter

DocMeta.setdocmeta!(CrystalUtils, :DocTestSetup, :(using CrystalUtils); recursive=true)

makedocs(;
    modules=[CrystalUtils],
    authors="Junfeng Qiao <qiaojunfeng@outlook.com> and contributors",
    sitename="CrystalUtils.jl",
    format=Documenter.HTML(;
        canonical="https://atomology.github.io/CrystalUtils.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/atomology/CrystalUtils.jl",
    devbranch="main",
)
