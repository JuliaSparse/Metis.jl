# convenienct function to flatten Vector of Vector
function flatten{T}(x::Vector{Vector{T}})
    out = T[]
    for y in x
        append!(out, y)
    end
    out
end

## Create the 0-based column pointers and row indices of the adjacency matrix
## as required by the Metis functions
function metisform(al::Graphs.GenericAdjacencyList)
    !Graphs.is_directed(al) || error("Metis functions require undirected graphs")
    @compat isa(al.vertices,UnitRange) && first(al.vertices) == 1 || error("Vertices must be numbered from 1")
    @compat length(al.adjlist), round.(Int32, cumsum(vcat(0, map(length, al.adjlist)))),
        round.(Int32, flatten(al.adjlist)) .- one(Int32)
end

## Create the 0-based column pointers and row indices of a symmetric sparse matrix
## after eliminating self-edges.  This is the form required by Metis functions.
function metisform(m::SparseMatrixCSC)
    issymmetric(m) || ishermitian(m) || error("m must be symmetric or Hermitian")

    ## copy m.rowval and m.colptr to Int32 vectors dropping diagonal elements
    adjncy = @compat sizehint!(Cint[],nnz(m))
    xadj = zeros(Cint,m.n+1)
    for j in 1:m.n
        count = 0
        for k in m.colptr[j] : (m.colptr[j+1] - 1)
            i = m.rowval[k]
            if i != j
                count += 1
                push!(adjncy,i-1)
            end
        end
        xadj[j+1] = xadj[j] + count
    end
    convert(Int32,m.n),xadj,adjncy
end

function metisform(g::LightGraphs.Graph)
    n = nv(g)
    xadj = zeros(Int32,n + 1)
    adjncy = sizehint!(Int32[],2*ne(g))
    for i in 1:n
        ein = [convert(Int32,x-1) for x in LightGraphs.neighbors(g, i)]
        xadj[i + 1] = xadj[i] + length(ein)
        append!(adjncy,ein)
    end
    n, xadj, adjncy
end

function testgraph(nm::String)
    pathnm = joinpath(dirname(@__FILE__), "..", "graphs", string(nm, ".graph"))
    ff = open(pathnm, "r")
    nvert, nedge = map(t -> parse(Int, t), split(readline(ff)))
    adjlist = Array{Vector{Int32}}(nvert)
    for i in 1:nvert adjlist[i] = map(t -> parse(Int32, t), split(readline(ff))) end
    close(ff)
    @compat GenericAdjacencyList{Int32,UnitRange{Int32},Vector{Vector{Int32}}}(false,
                                                                    map(Int32, 1:nvert),
                                                                    nedge,
                                                                    adjlist)
end
