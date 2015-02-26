## Create the 0-based column pointers and row indices of the adjacency matrix
## as required by the Metis functions
function metisform(al::Graphs.GenericAdjacencyList)
    !Graphs.is_directed(al) || error("Metis functions require undirected graphs")
    isa(al.vertices,Range1) && first(al.vertices) == 1 || error("Vertices must be numbered from 1")
    length(al.adjlist), int32(cumsum(vcat(0, map(length, al.adjlist)))),
        int32(vcat(al.adjlist...)) .- one(Int32)
end

## Create the 0-based column pointers and row indices of a symmetric sparse matrix
## after eliminating self-edges.  This is the form required by Metis functions.
function metisform(m::SparseMatrixCSC)
    issym(m) || ishermitian(m) || error("m must be symmetric or Hermitian")

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
        ein = [convert(Int32,dst(x)-1) for x in g.finclist[i]]
        xadj[i + 1] = xadj[i] + length(ein)
        append!(adjncy,ein)
    end
    n, xadj, adjncy
end

function testgraph(nm::ASCIIString)
    pathnm = joinpath(dirname(@__FILE__), "..", "graphs", string(nm, ".graph"))
    ff = open(pathnm, "r")
    nvert, nedge = map(int, split(readline(ff)))
    adjlist = Array(Vector{Int32}, nvert)
    for i in 1:nvert adjlist[i] = map(int32, split(readline(ff))) end
    close(ff)
    GenericAdjacencyList{Int32,Range1{Int32},Vector{Vector{Int32}}}(false,
                                                                    int32(1:nvert),
                                                                    nedge,
                                                                    adjlist)
end
