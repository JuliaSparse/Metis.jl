using Metis
using Compat.Test
using Compat.Random
import LightGraphs

@testset "Metis.permutation" begin
    srand(0)
    S = sprand(10, 10, 0.5); S = S + S'; fill(S.nzval, 1)
    perm, iperm = Metis.permutation(S)
    @test isperm(perm)
    @test isperm(iperm)
    @test S == S[perm,perm][iperm,iperm]
end

@testset "Metis.partition" begin
    srand(0)
    S = sprand(10, 10, 0.5); S = S + S'; fill(S.nzval, 1)
    T = LightGraphs.TutteGraph()
    for G in (S, T), alg in (:RECURSIVE, :KWAY), nparts in (3, 4)
        partition = Metis.partition(G, nparts, alg = alg)
        @test extrema(partition) == (1, nparts)
        @test all(x -> Compat.findfirst(==(x), partition) !== nothing, 1:nparts)
    end
end

@testset "Metis.separator" begin
    srand(0)
    S = sprand(10, 10, 0.5); S = S + S'; fill(S.nzval, 1)
    T = LightGraphs.TutteGraph()
    for G in (S, T)
        parts = Metis.separator(G)
        @test extrema(parts) == (1, 3)
    end
end