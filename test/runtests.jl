using Metis
using Random
using Test
using SparseArrays
using LinearAlgebra: Symmetric, Hermitian
import LightGraphs, Graphs

@testset "Metis.graph(::SparseMatrixCSC)" begin
    rng = MersenneTwister(0)
    S = sprand(rng, Int, 10, 10, 0.2); S = S + S'; fill!(S.nzval, 1)
    foreach(i -> S[i, i] = 0, 1:10); dropzeros!(S)
    g = Metis.graph(S)
    gw = Metis.graph(S; weights=true)
    @test g.xadj == gw.xadj
    @test g.adjncy == gw.adjncy
    @test g.adjwgt == C_NULL
    @test gw.adjwgt == ones(Int, length(gw.adjncy))
    G  = SparseMatrixCSC(10, 10, g.xadj, g.adjncy, ones(Int, length(g.adjncy)))
    GW = SparseMatrixCSC(10, 10, gw.xadj, gw.adjncy, gw.adjwgt)
    @test iszero(S - G)
    @test iszero(S - GW)
end

@testset "Metis.graph(::Union{Hermitian, Symmetric})" begin
    rng = MersenneTwister(0)
    for T in (Symmetric, Hermitian), uplo in (:U, :L)
        S = sprand(rng, Int, 10, 10, 0.2); fill!(S.nzval, 1)
        TS = T(S, uplo)
        CSCS = SparseMatrixCSC(TS)
        @test TS == CSCS
        g1 = Metis.graph(TS)
        g2 = Metis.graph(CSCS)
        @test g1.nvtxs == g2.nvtxs
        @test g1.xadj == g2.xadj
        @test g1.adjncy == g2.adjncy
        @test g1.vwgt == g2.vwgt == C_NULL
        @test g1.adjwgt == g2.adjwgt == C_NULL
        @test_throws ArgumentError Metis.graph(TS; weights = true)
    end
end

@testset "Metis.graph(::SimpleWeightedDiGraph)" begin
    if isdefined(Base, :get_extension)
        import SimpleWeightedGraphs
        # Issue #56: the weight must be that of the edge from j to i
        g = Metis.graph(SimpleWeightedGraphs.SimpleWeightedDiGraph([1], [2], [1]))
        @test g.xadj == [1, 2, 2]
        @test g.adjncy == [2]
        @test g.adjwgt == [1]
        g = Metis.graph(SimpleWeightedGraphs.SimpleWeightedDiGraph([1, 2], [2, 1], [3, 5]))
        @test g.xadj == [1, 2, 3]
        @test g.adjncy == [2, 1]
        @test g.adjwgt == [3, 5]
    end
end

@testset "Edge weight validation" begin
    # Stored zeros are kept by `sparse`, so they are seen as edges
    for w in (0, 0.0, -5, -5.0, 0.5, 1 + 2im)
        S = sparse([1, 2], [2, 1], [w, conj(w)])
        @test_throws ArgumentError Metis.graph(S; weights=true)
        # The weights are only validated when they are used
        @test Metis.graph(S).adjwgt == C_NULL
    end
    # The message names the offending entry of the matrix
    A = sparse([2, 1], [1, 2], [-5, 3])
    err = try
        Metis.graph(A; weights=true, check_hermitian=false)
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("got weight -5 for edge (2, 1)", err.msg)
    # Real integers stored as other types are valid
    for w in (2, 2.0, 2 + 0im, 2.0 + 0.0im, 2//1)
        @test Metis.graph(sparse([1, 2], [2, 1], [w, w]); weights=true).adjwgt == [2, 2]
    end
    # Self edges are skipped, so diagonal entries are not validated
    S = sparse([1, 2, 1, 2], [1, 2, 2, 1], [0, -1, 3, 3])
    g = Metis.graph(S; weights=true)
    @test g.adjncy == [2, 1]
    @test g.adjwgt == [3, 3]
    # The graph extensions use the same validation
    if isdefined(Base, :get_extension)
        import SimpleWeightedGraphs
        for w in (-5, 0.5)
            G = SimpleWeightedGraphs.SimpleWeightedGraph([1], [2], [w])
            @test_throws ArgumentError Metis.graph(G)
        end
        @test Metis.graph(SimpleWeightedGraphs.SimpleWeightedGraph([1], [2], [2.0])).adjwgt == [2, 2]
    end
end

@testset "Metis.permutation" begin
    rng = MersenneTwister(0)
    S = sprand(rng, 10, 10, 0.5); S = S + S'; fill!(S.nzval, 1)
    perm, iperm = Metis.permutation(S)
    @test isperm(perm)
    @test isperm(iperm)
    @test S == S[perm,perm][iperm,iperm]
end

@testset "Metis.partition" begin
    T = Graphs.smallgraph(:tutte)
    S = Graphs.adjacency_matrix(T)
    SG = Metis.graph(S; weights=true)
    TL = LightGraphs.smallgraph(:tutte)
    for G in (S, T, TL, SG), alg in (:RECURSIVE, :KWAY), nparts in (1, 3, 4)
        partition = Metis.partition(G, nparts, alg = alg)
        @test extrema(partition) == (1, nparts)
        @test all(x -> findfirst(==(x), partition) !== nothing, 1:nparts)
    end
    # With weights
    if isdefined(Base, :get_extension)
        import SimpleWeightedGraphs
        G1 = SimpleWeightedGraphs.SimpleWeightedGraph(Graphs.nv(T))
        G2 = SimpleWeightedGraphs.SimpleWeightedGraph(Graphs.nv(T))
        for edge in Graphs.edges(T)
            i, j = Tuple(edge)
            SimpleWeightedGraphs.add_edge!(G1, i, j, 1)
            SimpleWeightedGraphs.add_edge!(G2, i, j, i+j)
        end
        for alg in (:RECURSIVE, :KWAY), nparts in (3, 4)
            unwpartition = Metis.partition(T, nparts, alg = alg)
            partition = Metis.partition(G1, nparts, alg = alg)
            @test partition == unwpartition
            partition = Metis.partition(G2, nparts, alg = alg)
            @test partition != unwpartition
            @test extrema(partition) == (1, nparts)
            @test all(x -> findfirst(==(x), partition) !== nothing, 1:nparts)
        end
    end
end

@testset "Metis.separator" begin
    T = Graphs.smallgraph(:tutte)
    S = Graphs.adjacency_matrix(T)
    SG = Metis.graph(S; weights=true)
    T = Graphs.smallgraph(:tutte)
    TL = LightGraphs.smallgraph(:tutte)
    for G in (S, T, TL, SG)
        parts = Metis.separator(G)
        @test extrema(parts) == (1, 3)
    end
end
