using Graphs
using Metis
using Compat
using Base.Test

const copter2 = Metis.testgraph("copter2");
perm, iperm = nodeND(copter2)
@test all(invperm(perm) .== iperm)

srand(12321)
const a = convert(SparseMatrixCSC{Float64,Int32},sprand(1000,100,0.01))
const ata = a'a
perm, iperm = nodeND(ata,3)
@test all(invperm(perm) .== iperm)
perm, iperm = nodeND(ata,0)             # reset verbosity level

function counts{T<:Integer}(v::Vector{T})
    ans = zeros(Int,maximum(v))
    for vv in v
        ans[vv] += 1
    end
    ans
end

objval, part = partGraphKway(copter2, 6)
@test extrema(part) == (1,6)            # subset values are in the correct range
@test all(i->bool(findfirst(part,i)),1:6) # the 6 subsets are non-empty

# the sizes of partition subsets are not consistent across OS's
#@test counts(part) == [9076,9374,9384,9523,8978,9141]

function interface(part::Vector, g::Graphs.GenericAdjacencyList)
    (nv = length(part)) == length(g.vertices) || error("partition length != # of vertices")
    conn = falses(nv)                   # vertex connected to another subset?
    for i in 1:nv
        pp = part[i]
        conn[i] = any(j -> part[j] != pp, g.adjlist[i])
    end
    countnz(conn)
end

@test interface(part,copter2) < 6000    # 5907 on an Ubuntu system

objval, part = partGraphRecursive(copter2,6)
@test extrema(part) == (1,6)
@test all(i->bool(findfirst(part,i)),1:6)
@test interface(part,copter2) < 6000

const mdual = Metis.testgraph("mdual")
objval, part = partGraphKway(mdual, 10)
@test extrema(part) == (1,10)
@test all(i->bool(findfirst(part,i)),1:10)
@test interface(part,mdual) < 19000     # 18263 on an Ubuntu system

sizes, copterPart = vertexSep(copter2)

function testGraphPart(g,part)
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
        elseif partVal != 2
            println("Vertex $i assigned to set $partVal")
            validPart = false
        end
    end
    validPart
end

@test testGraphPart(copter2,copterPart)

function appendel(I,J,V,i,j,v)
    push!(I,i)
    push!(J,j)
    push!(V,v)
end

function laplacian2d{T<:Integer}(nx::T,ny::T)
    n = nx*ny
    nzest = 5n
    I = @compat sizehint!(Int32[],nzest)
    J = @compat sizehint!(Int32[],nzest)
    V = @compat sizehint!(Float64[],nzest)
    for x in 1:nx
        for y in 1:ny
            s = x + (y-1)*nx
            appendel(I,J,V,s,s,2)
            x > 1 && appendel(I,J,V,s,s-1,-1)
            y > 1 && appendel(I,J,V,s,s-nx,-1)
        end
    end
    A = sparse(I,J,V,n,n)
    A + A'
end

const A = laplacian2d(100,110)

sizes, matPart = vertexSep(A)

function testMatPart(m,part)
    validPart = true
    for i in 1:length(part)
        partVal = part[i]
        if partVal == 0
            for k in m.colptr[i]:m.colptr[i+1]-1
                j = m.rowval[k]
                if part[j] == 1
                    println("Nonzero ($i,$j) connects sets 0 and 1")
                    validPart = false
                end
            end
        elseif partVal == 1
            for k in m.colptr[i]:m.colptr[i+1]-1
                j = m.rowval[k]
                if part[j] == 0
                    println("Nonzero ($i,$j) connects sets 1 and 0")
                    validPart = false
                end
            end
        elseif partVal != 2
            println("Vertex $i assigned to set $partVal")
            validPart = false
        end
    end
    validPart
end

@test testMatPart(A,matPart)

g = LightGraphs.TutteGraph()
x, y = partGraphKway(g, 6)
@test maximum(y) == 6
@test length(y) == LightGraphs.nv(g)
