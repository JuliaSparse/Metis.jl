using Graphs
using Metis
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
@test counts(part) == [9076,9374,9384,9523,8978,9141]

objval, part = partGraphRecursive(copter2,6)
@test counts(part) == [9076,9374,9384,9523,8978,9141]

const mdual = Metis.testgraph("mdual")
objval, part = partGraphKway(mdual, 10)
@test counts(part) == [25789,25731,25790,25998,25728,25724,25722,26061,25995,26031]

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
    I = sizehint(Int32[],nzest)
    J = sizehint(Int32[],nzest)
    V = sizehint(Float64[],nzest)
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
