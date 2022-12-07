using Clang.Generators
using METIS_jll

cd(@__DIR__)

metis_include_dir = normpath(METIS_jll.artifact_dir, "include")

options = load_options(joinpath(@__DIR__, "generator.toml"))

args = get_default_args()
push!(args, "-I$(metis_include_dir)")

# Compiler flags from Yggdrasil (??)
# https://github.com/JuliaPackaging/Yggdrasil/blob/ed9cc4d0cec8d3bbc973a47cb227c42180261782/M/METIS/METIS%405/build_tarballs.jl
push!(args, "-DSHARED=1")
# push!(args, "-DIDXTYPEWIDTH=32") # 32 and 64
# push!(args, "-DREALTYPEWIDTH=32") # 32 and 64

headers = joinpath.(metis_include_dir, [
    "metis.h",
])

ctx = create_context(headers, args, options)

build!(ctx)
