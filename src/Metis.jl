__precompile__(true)

module Metis

# Load libmetis with BinaryProvider
__init__() = check_deps()
let depsfile = joinpath(@__DIR__, "..", "deps", "deps.jl")
    if isfile(depsfile)
        include(depsfile)
    else
        error("$(depsfile) does not exist, Please re-run Pkg.build(\"Metis\"), and restart Julia.")
    end
end

end # module
