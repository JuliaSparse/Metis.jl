## Create a 0-based CSR redundant representation of an undirected graph or hermitian sparse matrix

function mkadjCSR(al::GenericAdjacencyList)
    !is_directed(al) || error("Metis functions require undirected graphs")
    isa(al.vertices,Range1) && first(al.vertices) == 1 || error("Vertices must be numbered from 1")
    length(al.adjlist), int32(cumsum(vcat(0, map(length, al.adjlist)))),
        int32(vcat(al.adjlist...)) .- one(Int32)
end

function mkadjCSR(m::SparseMatrixCSC)
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


function testgraph(nm::ASCIIString)
    pathnm = Pkg.dir("Metis", "graphs", string(nm, ".graph"))
    ff = open(pathnm, "r")
    nvert, nedge = map(int, split(readline(ff)))
    adjlist = Array(Vector{Int32}, nvert)
    for i in 1:nvert adjlist[i] = map(int32, split(readline(ff))) end
    GenericAdjacencyList{Int32,Range1{Int32},Vector{Vector{Int32}}}(false,
                                                                    int32(1:nvert),
                                                                    nedge,
                                                                    adjlist)
end

