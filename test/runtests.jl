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
