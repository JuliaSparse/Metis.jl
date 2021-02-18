using Metis
using Random
using Test
using SparseArrays
import LightGraphs
import SimpleWeightedGraphs

@testset "Metis.permutation" begin
    Random.seed!(0)
    S = sprand(10, 10, 0.5); S = S + S'; fill(S.nzval, 1)
    perm, iperm = Metis.permutation(S)
    @test isperm(perm)
    @test isperm(iperm)
    @test S == S[perm,perm][iperm,iperm]
end

@testset "Metis.partition" begin
    Random.seed!(0)
    S = sprand(10, 10, 0.5); S = S + S'; fill(S.nzval, 1)
    T = LightGraphs.smallgraph(:tutte)
    for G in (S, T), alg in (:RECURSIVE, :KWAY), nparts in (3, 4)
        partition = Metis.partition(G, nparts, alg = alg)
        @test extrema(partition) == (1, nparts)
        @test all(x -> findfirst(==(x), partition) !== nothing, 1:nparts)
    end
end

@testset "Metis.partition with edge weights" begin
    T = LightGraphs.smallgraph(:tutte)
    G = SimpleWeightedGraphs.SimpleWeightedGraph(LightGraphs.nv(T))
    for i in 1:LightGraphs.nv(T)
        for j in LightGraphs.neighbors(T, i)
            SimpleWeightedGraphs.add_edge!(G, i, j, 1)
        end
    end
    for alg in (:RECURSIVE, :KWAY), nparts in (3, 4)
        partition = Metis.partition(G, nparts, alg = alg)
        @test extrema(partition) == (1, nparts)
        @test all(x -> findfirst(==(x), partition) !== nothing, 1:nparts)
    end
end

@testset "Metis.separator" begin
    Random.seed!(0)
    S = sprand(10, 10, 0.5); S = S + S'; fill(S.nzval, 1)
    T = LightGraphs.smallgraph(:tutte)
    for G in (S, T)
        parts = Metis.separator(G)
        @test extrema(parts) == (1, 3)
    end
end
