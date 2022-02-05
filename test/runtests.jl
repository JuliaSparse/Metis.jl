using Metis
using Random
using Test
using SparseArrays
import LightGraphs, Graphs

@testset "Metis.permutation" begin
    rng = MersenneTwister(0)
    S = sprand(rng, 10, 10, 0.5); S = S + S'; fill(S.nzval, 1)
    perm, iperm = Metis.permutation(S)
    @test isperm(perm)
    @test isperm(iperm)
    @test S == S[perm,perm][iperm,iperm]
end

@testset "Metis.partition" begin
    rng = MersenneTwister(0)
    S = sprand(rng, 10, 10, 0.5); S = S + S'; fill(S.nzval, 1)
    T = Graphs.smallgraph(:tutte)
    TL = LightGraphs.smallgraph(:tutte)
    for G in (S, T, TL), alg in (:RECURSIVE, :KWAY), nparts in (3, 4)
        partition = Metis.partition(G, nparts, alg = alg)
        @test extrema(partition) == (1, nparts)
        @test all(x -> findfirst(==(x), partition) !== nothing, 1:nparts)
    end
end

@testset "Metis.separator" begin
    rng = MersenneTwister(0)
    S = sprand(rng, 10, 10, 0.5); S = S + S'; fill(S.nzval, 1)
    T = Graphs.smallgraph(:tutte)
    TL = LightGraphs.smallgraph(:tutte)
    for G in (S, T, TL)
        parts = Metis.separator(G)
        @test extrema(parts) == (1, 3)
    end
end
