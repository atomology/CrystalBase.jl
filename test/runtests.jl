using CrystalBase
using Documenter

DocMeta.setdocmeta!(
    CrystalBase,
    :DocTestSetup, :(using CrystalBase);
    recursive=true,
)
doctest(
    CrystalBase,
    fix=true,  # update all the output in `jldoctest`
)

using TestItemRunner

@run_package_tests verbose = true
