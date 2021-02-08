using Metis
using Random
using Test
using SparseArrays
import LightGraphs

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

@testset "Metis.separator" begin
    Random.seed!(0)
    S = sprand(10, 10, 0.5); S = S + S'; fill(S.nzval, 1)
    T = LightGraphs.smallgraph(:tutte)
    for G in (S, T)
        parts = Metis.separator(G)
        @test extrema(parts) == (1, 3)
    end
end

@testset "Metis.MeshToDual" begin
    # Connectivity matrix of 2D mesh of 4 triangles.
    cells = [1 2 2 3;
             2 5 3 6;
             4 4 5 5]
    # Tell Metis to consider a cell connected if it shares two vertices.
    ncommon = 2

    # Get a Graph object. The number of graph verticies is the number of mesh
    # elements.
    G = Metis.graph(cells, ncommon)

    # Partition the graph.
    parts = Metis.partition(G, 2)
    # First two elements should be in the same partition.
    @test parts[1] == parts[2]
    # Last two elements should be in the same partition.
    @test parts[3] == parts[4]
end
