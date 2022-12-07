module LibMetis

# Clang.jl auto-generated bindings
include("../lib/LibMetis.jl")

# Error checking
struct MetisError <: Exception
    fn::Symbol
    code::Cint
end

function Base.showerror(io::IO, me::MetisError)
    print(io, "MetisError: LibMETIS.$(me.fn) returned error code $(me.code): ")
    if me.code == METIS_ERROR_INPUT
        print(io, "input error")
    elseif me.code == METIS_ERROR_MEMORY
        print(io, "could not allocate the required memory")
    else
        print(io, "unknown error")
    end
    return
end

# Macro for checking return codes of ccalls
macro check(arg)
    Meta.isexpr(arg, :call) || throw(ArgumentError("wrong usage of @check"))
    return quote
        r = $(esc(arg))
        if r != METIS_OK
            throw(MetisError($(QuoteNode(arg.args[1])), r))
        end
        r
    end
end

# Export everything with METIS_ prefix
for name in names(@__MODULE__; all=true)
    if startswith(string(name), "METIS_")
        @eval export $name
    end
end

end # module LibMetis
