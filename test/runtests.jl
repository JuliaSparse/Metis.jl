using Graphs
using LightGraphs
using Metis
using Compat
using Base.Test

# ==========================================================
# Initialize test graphs
# ==========================================================

"""
```
laplacian2d(nx, ny) -> L
```

Construct a 2D discrete Laplacian operator on an `nx` x `ny` grid. `L`
is a sparse matrix.
"""
function laplacian2d{T<:Integer}(nx::T,ny::T)

    # Initialize storage for graph
    n = nx*ny
    nzest = 5*n
    I = @compat sizehint!(Int32[],nzest)
    J = @compat sizehint!(Int32[],nzest)
    V = @compat sizehint!(Float64[],nzest)

    # Add edge to graph
    function addEdge(i,j,v)
        push!(I,i)
        push!(J,j)
        push!(V,v)
    end

    # Construct discrete Laplacian
    for y in 1:ny
        for x in 1:nx
            s = x + (y-1)*nx
            addEdge(s,s,4)
            x > 1  && addEdge(s-1,s,-1)
            x < nx && addEdge(s+1,s,-1)
            y > 1  && addEdge(s-nx,s,-1)
            y < ny && addEdge(s+nx,s,-1)
        end
    end
    return sparse(I,J,V,n,n)
end

# Construct SparseMatrixCSC examples
srand(12321)
const a = convert(SparseMatrixCSC{Float64,Int32},
                  sprand(1000,100,0.01))
const ata = a'a
const laplace = laplacian2d(123,89)
const path_CSC = spdiagm((ones(Int,7),ones(Int,7)), (-1,1), 8, 8)
const wheel = sparse([0   100 100 1   1   1   1   1  ;
                      100 0   100 0   0   0   0   1  ;
                      100 100 0   1   0   0   0   0  ;
                      1   0   1   0   100 0   0   0  ;
                      1   0   0   100 0   100 0   0  ;
                      1   0   0   0   100 0   1   0  ;
                      1   0   0   0   0   1   0   100;
                      1   1   0   0   0   0   100 0  ])

# Construct Graphs.GenericAdjacencyList examples
const copter2     = Metis.testgraph("copter2");
const mdual       = Metis.testgraph("mdual")
const path_Graphs = Metis.testgraph("path")

# Construct LightGraphs.Graph examples
const tutte            = LightGraphs.TutteGraph()
const path_LightGraphs = LightGraphs.Graph(LightGraphs.PathDiGraph(8))

# Vertex weights for path graph and wheel graph
const vwgt_path  = [4,3,4,3,4,1,1,1]
const vwgt_wheel = [1,1,1,1,1,1,1,2]

# ==========================================================
# Tests with nodeND
# ==========================================================

# copter2
perm, iperm = nodeND(copter2)
@test all(invperm(perm) .== iperm)

# tutte
perm, iperm = nodeND(tutte)
@test all(invperm(perm) .== iperm)

# ata
perm, iperm = nodeND(ata)
@test all(invperm(perm) .== iperm)

# ==========================================================
# Tests with partGraphKway
# ==========================================================

"""
```
interface(g, part) -> Int
```

Compute the number of vertices that lie on partition interfaces. `g`
can be a Graphs package adjacency list or a LightGraphs package graph
and `part` is a partition vector.
"""
# Compute interface size on Graphs.GenericAdjacencyList
function interface(g::Graphs.GenericAdjacencyList, part::Vector)
    n = length(part)
    (n==length(g.vertices)) || error("partition length != # of vertices")
    interfaceSize = 0
    for i in 1:n
        if any(j -> (part[j]!=part[i]), g.adjlist[i])
            interfaceSize += 1
        end
    end
    return interfaceSize
end

# Compute interface size on LightGraphs.Graph
function interface(g::LightGraphs.Graph, part::Vector)
    n = nv(g)
    (n==length(part)) || error("partition length != # of vertices")
    interfaceSize = 0
    for i in 1:n
        if any(j -> (part[j]!=part[i]), map(dst,g.finclist[i]))
            interfaceSize += 1
        end
    end
    return interfaceSize
end

# copter2
objval, part = partGraphKway(copter2, 6)
@test extrema(part) == (1,6)
@test all(i->findfirst(part,i) != 0,1:6)
@test interface(copter2,part) < 6000    # 5907 on an Ubuntu system

# mdual
objval, part = partGraphKway(mdual, 10)
@test extrema(part) == (1,10)
@test all(i->findfirst(part,i) != 0,1:10)
@test interface(mdual,part) < 19000     # 18263 on an Ubuntu system

# tutte
objval, part = partGraphKway(tutte, 6)
@test extrema(part) == (1,6)
@test length(part) == LightGraphs.nv(tutte)
@test all(i-> findfirst(part,i) != 0,1:6)

# path_CSC
objval, part = partGraphKway(path_CSC, 3, vwgt=vwgt_path)
@test objval == 2
@test ((part.==1) == [1,1,0,0,0,0,0,0]
       || (part.==1) == [0,0,1,1,0,0,0,0]
       || (part.==1) == [0,0,0,0,1,1,1,1])
@test ((part.==2) == [1,1,0,0,0,0,0,0]
       || (part.==2) == [0,0,1,1,0,0,0,0]
       || (part.==2) == [0,0,0,0,1,1,1,1])
@test ((part.==3) == [1,1,0,0,0,0,0,0]
       || (part.==3) == [0,0,1,1,0,0,0,0]
       || (part.==3) == [0,0,0,0,1,1,1,1])

# path_Graphs
objval, part = partGraphKway(path_Graphs, 3, vwgt=vwgt_path)
@test objval == 2
@test ((part.==1) == [1,1,0,0,0,0,0,0]
       || (part.==1) == [0,0,1,1,0,0,0,0]
       || (part.==1) == [0,0,0,0,1,1,1,1])
@test ((part.==2) == [1,1,0,0,0,0,0,0]
       || (part.==2) == [0,0,1,1,0,0,0,0]
       || (part.==2) == [0,0,0,0,1,1,1,1])
@test ((part.==3) == [1,1,0,0,0,0,0,0]
       || (part.==3) == [0,0,1,1,0,0,0,0]
       || (part.==3) == [0,0,0,0,1,1,1,1])

# path_LightGraphs
objval, part = partGraphKway(path_LightGraphs, 3, vwgt=vwgt_path)
@test objval == 2
@test ((part.==1) == [1,1,0,0,0,0,0,0]
       || (part.==1) == [0,0,1,1,0,0,0,0]
       || (part.==1) == [0,0,0,0,1,1,1,1])
@test ((part.==2) == [1,1,0,0,0,0,0,0]
       || (part.==2) == [0,0,1,1,0,0,0,0]
       || (part.==2) == [0,0,0,0,1,1,1,1])
@test ((part.==3) == [1,1,0,0,0,0,0,0]
       || (part.==3) == [0,0,1,1,0,0,0,0]
       || (part.==3) == [0,0,0,0,1,1,1,1])

# wheel
objval, part = partGraphKway(wheel, 3, adjwgt=true, vwgt=vwgt_wheel)
@test objval == 8
@test ((part.==1) == [1,1,1,0,0,0,0,0]
       || (part.==1) == [0,0,0,1,1,1,0,0]
       || (part.==1) == [0,0,0,0,0,0,1,1])
@test ((part.==2) == [1,1,1,0,0,0,0,0]
       || (part.==2) == [0,0,0,1,1,1,0,0]
       || (part.==2) == [0,0,0,0,0,0,1,1])
@test ((part.==3) == [1,1,1,0,0,0,0,0]
       || (part.==3) == [0,0,0,1,1,1,0,0]
       || (part.==3) == [0,0,0,0,0,0,1,1])

# ==========================================================
# Tests with partGraphRecursive
# ==========================================================

# copter2
objval, part = partGraphRecursive(copter2,6)
@test extrema(part) == (1,6)
@test all(i->findfirst(part,i) != 0,1:6)
@test interface(copter2,part) < 6000    # 5907 on an Ubuntu system

# mdual
objval, part = partGraphKway(mdual, 10)
@test extrema(part) == (1,10)
@test all(i->findfirst(part,i) != 0,1:10)
@test interface(mdual,part) < 19000     # 18263 on an Ubuntu system

# tutte
objval, part = partGraphKway(tutte, 6)
@test extrema(part) == (1,6)
@test length(part) == LightGraphs.nv(tutte)
@test all(i-> findfirst(part,i) != 0,1:6)

# ==========================================================
# Tests with vertexSep
# ==========================================================

"""
```
testVertexSep(g, part) -> Bool
```

Test whether a graph partition describes a vertex separator. `g` can
be a Graphs package adjacency list, LightGraphs graph, or adjacency
matrix in CSC format. `part` is a partition vector with three
parts. The third part is being evaluated as a vertex separator between
the first two parts.
"""
# Test a vertex separator on Graphs.GenericAdjacencyList
function testVertexSep(g::Graphs.GenericAdjacencyList, part::Vector)
    validPart = true
    for i in 1:length(part)
        partVal = part[i]
        if partVal == 0
            for j in g.adjlist[i]
                if part[j] == 1
                    println("Edge ($i,$j) connects sets 0 and 1")
                    validPart = false
                end
            end
        elseif partVal == 1
            for j in g.adjlist[i]
                if part[j] == 0
                    println("Edge ($i,$j) connects sets 1 and 0")
                    validPart = false
                end
            end
        elseif partVal == 2
            continue
        else
            println("Vertex $i assigned to set $partVal")
            validPart = false
        end
    end
    return validPart
end

# Test a vertex separator on LightGraphs.Graph
function testVertexSep(g::LightGraphs.Graph, part::Vector)
    validPart = true
    for i in 1:length(part)
        partVal = part[i]
        if partVal == 0
            for j in g.fadjlist[i]
                if part[j] == 1
                    println("Edge ($i,$j) connects sets 0 and 1")
                    validPart = false
                end
            end
        elseif partVal == 1
            for j in g.fadjlist[i]
                if part[j] == 0
                    println("Edge ($i,$j) connects sets 1 and 0")
                    validPart = false
                end
            end
        elseif partVal == 2
            continue
        else
            println("Vertex $i assigned to set $partVal")
            validPart = false
        end
    end
    return validPart
end

# Test a vertex separator on SparseMatrixCSC
function testVertexSep(m::SparseMatrixCSC, part::Vector)
    validPart = true
    for i in 1:length(part)
        partVal = part[i]
        if partVal == 0
            for k in m.colptr[i]:m.colptr[i+1]-1
                j = m.rowval[k]
                if part[j] == 1
                    println("Edge ($i,$j) connects sets 0 and 1")
                    validPart = false
                end
            end
        elseif partVal == 1
            for k in m.colptr[i]:m.colptr[i+1]-1
                j = m.rowval[k]
                if part[j] == 0
                    println("Edge ($i,$j) connects sets 1 and 0")
                    validPart = false
                end
            end
        elseif partVal == 2
            continue
        else
            println("Vertex $i assigned to set $partVal")
            validPart = false
        end
    end
    return validPart
end

# copter2
sizes, part = vertexSep(copter2)
@test sizes[1]!=0 && sizes[2]!=0
@test testVertexSep(copter2,part)

# tutte
sizes, part = vertexSep(tutte)
@test sizes[1]!=0 && sizes[2]!=0
@test testVertexSep(tutte,part)

# laplace
sizes, part = vertexSep(laplace)
@test sizes[1]!=0 && sizes[2]!=0
@test testVertexSep(laplace,part)
