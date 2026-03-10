using CrystalUtils
using Documenter

DocMeta.setdocmeta!(
    CrystalUtils,
    :DocTestSetup, :(using CrystalUtils);
    recursive=true,
)
doctest(
    CrystalUtils,
    fix=true,  # update all the output in `jldoctest`
)

using TestItemRunner

@run_package_tests verbose = true
