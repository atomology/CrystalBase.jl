using TestItemRunner
using CrystalBase
using Documenter

# Get filter strings from command line arguments if provided
# Usage: julia --project test/runtests.jl "kpath.jl" "type.jl"
filter_names = isempty(ARGS) ? nothing : ARGS

if isnothing(filter_names)
    println("Running all tests...")

    @run_package_tests verbose = true

    DocMeta.setdocmeta!(CrystalBase, :DocTestSetup, :(using CrystalBase); recursive=true)
    doctest(
        CrystalBase,
        # fix=true,  # update all the output in `jldoctest`
    )
else
    println("Running specific tests: $(join(filter_names, ", "))")

    @run_package_tests verbose = true filter =
        ti -> any(name -> endswith(ti.filename, name), filter_names)
end
