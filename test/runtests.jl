using Graphs
using Metis
using Base.Test

copter2 = Metis.testgraph("copter2");
perm, iperm = nodeND(copter2)
@test all(invperm(perm) .== iperm)

srand(12321)
a = convert(SparseMatrixCSC{Float64,Int32},sprand(1000,100,0.01))
ata = a'a
perm, iperm = nodeND(ata,3)
@test all(invperm(perm) .== iperm)

perm, iperm = nodeND(ata,0)             # reset the verbosity
@test all(invperm(perm) .== iperm)

function counts{T<:Integer}(v::Vector{T},k::Int)
    ans = zeros(Int,k)
    for i in 1:length(v)
        ans[v[i]] += 1
    end
    ans
end


objval, part = partGraphKway(copter2, 6)
@test counts(part,6) == [9076,9374,9384,9523,8978,9141]

objval, part = partGraphRecursive(copter2,6)
@test counts(part,6) == [9076,9374,9384,9523,8978,9141]

mdual = Metis.testgraph("mdual")
objval, part = partGraphKway(mdual, 10)
@test counts(part,10) == [25789,25731,25790,25998,25728,25724,25722,26061,25995,26031]

copter2 = Metis.testgraph("copter2")
# There does not seem to be a simple way to copy a graph...
copter2Copy = Metis.testgraph("copter2")
sepSize, copterPart = vertexSep(copter2Copy)

function testGraphPart(g,part)
  n = length(copter2.vertices)
  validPart = true
  for i=1:n
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

nx = 100
ny = 110
A = speye(nx*ny,nx*ny)
ACopy = speye(nx*ny,nx*ny)
# NOTE: At this time, copy(A) does not produce a deep copy
for x=1:nx
  for y=1:ny
    s = x + (y-1)*nx
    A[s,s] = 4
    ACopy[s,s] = 4
    if x > 1
      A[s,s-1] = -1
      ACopy[s,s-1] = -1
    end
    if x < nx
      A[s,s+1] = -1
      ACopy[s,s+1] = -1
    end
    if y > 1
      A[s,s-nx] = -1
      ACopy[s,s-nx] = -1
    end
    if y < ny
      A[s,s+nx] = -1
      ACopy[s,s+nx] = -1
    end
  end
end

sepSize, matPart = vertexSep(A) 

function testMatPart(m,part)
  validPart = true
  n = m.n
  validPart = true
  for i=1:n
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

@test testMatPart(ACopy,matPart)
