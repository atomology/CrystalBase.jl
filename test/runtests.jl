using TestItemRunner
using CrystalBase
using Documenter

# Get filter string from command line arguments if provided
# Usage: julia --project test/runtests.jl "Name of Test"
filter_name = length(ARGS) > 0 ? ARGS[1] : nothing

if isnothing(filter_name)
    println("Running all tests...")

    DocMeta.setdocmeta!(CrystalBase, :DocTestSetup, :(using CrystalBase); recursive=true)
    doctest(
        CrystalBase,
        # fix=true,  # update all the output in `jldoctest`
    )

    @run_package_tests verbose = true
else
    println("Running specific test: $filter_name")

    @run_package_tests verbose = true filter = ti -> endswith(ti.filename, filter_name)
end
