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

@testset "Metis.graph(::SparseMatrixCSC) rejects non-square matrices" begin
    # Compare the full exception: a non-square matrix already fails the Hermitian check.
    for check_hermitian in (true, false)
        tall = sparse([10], [5], [1], 10, 5)
        @test_throws ArgumentError("matrix must be square, got size (10, 5)") Metis.graph(tall; check_hermitian)
        wide = sparse([1], [9], [1], 3, 10)
        @test_throws ArgumentError("matrix must be square, got size (3, 10)") Metis.graph(wide; check_hermitian)
    end
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

@testset "Metis.graph rejects unknown keywords for graph types" begin
    GG = Graphs.smallgraph(:house)
    LG = LightGraphs.smallgraph(:house)
    Gs = Any[GG, LG]
    if isdefined(Base, :get_extension)
        import SimpleWeightedGraphs
        push!(Gs, SimpleWeightedGraphs.SimpleWeightedGraph(GG))
    end
    for G in Gs
        @test_throws MethodError Metis.graph(G; weight=true)
        @test_throws MethodError Metis.graph(G; check_hermitian=false)
        g = Metis.graph(G; weights=false)
        @test g.adjwgt == C_NULL
        @test length(g.adjncy) == 2 * Graphs.ne(GG)
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

@testset "Metis.graph(::AbstractSimpleGraph; weights=true)" begin
    # Built by hand since `LightGraphs.path_graph` only exists in newer LightGraphs versions
    LG = LightGraphs.SimpleGraph(3)
    LightGraphs.add_edge!(LG, 1, 2)
    LightGraphs.add_edge!(LG, 2, 3)
    for (G, T) in ((Graphs.path_graph(3), Graphs.smallgraph(:tutte)),
                   (LG, LightGraphs.smallgraph(:tutte)))
        g = Metis.graph(G)
        gw = Metis.graph(G; weights=true)
        @test g.adjwgt == C_NULL
        @test gw.adjwgt == [1, 1, 1, 1]
        @test gw.xadj == g.xadj
        @test gw.adjncy == g.adjncy
        for alg in (:RECURSIVE, :KWAY)
            # Unit weights partition like no weights. Only the larger graph is checked for
            # non-empty parts: :KWAY puts this path graph in one part, with or without weights.
            partition = Metis.partition(gw, 2, alg = alg)
            @test all(in(1:2), partition)
            @test partition == Metis.partition(G, 2, alg = alg)
            partition = Metis.partition(Metis.graph(T; weights=true), 3, alg = alg)
            @test extrema(partition) == (1, 3)
            @test partition == Metis.partition(T, 3, alg = alg)
        end
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

@testset "Metis.partition argument validation" begin
    S = Graphs.adjacency_matrix(Graphs.cycle_graph(8))
    for alg in (:RECURSIVE, :KWAY), nparts in (0, -1)
        @test_throws ArgumentError Metis.partition(S, nparts; alg=alg)
    end
    for nparts in (1, 2)
        @test_throws ArgumentError Metis.partition(S, nparts; alg=:BOGUS)
    end
    # More parts than vertices is an error, as many parts as vertices is fine
    C4 = Graphs.adjacency_matrix(Graphs.cycle_graph(4))
    for alg in (:RECURSIVE, :KWAY)
        @test_throws ArgumentError Metis.partition(C4, 5; alg=alg)
        part = Metis.partition(C4, 4; alg=alg)
        @test length(part) == 4
        @test all(in(1:4), part)
    end
    # A graph without vertices can only be partitioned into one part
    E = spzeros(Int, 0, 0)
    for alg in (:RECURSIVE, :KWAY)
        part = Metis.partition(E, 1; alg=alg)
        @test part isa Vector
        @test isempty(part)
        @test_throws ArgumentError Metis.partition(E, 2; alg=alg)
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

@testset "Edge weights of SimpleWeightedGraphs that are invalid for Metis" begin
    if isdefined(Base, :get_extension)
        import SimpleWeightedGraphs
        g = SimpleWeightedGraphs.SimpleWeightedGraph([1, 2, 3], [2, 3, 4], [0.5, 1.5, 2.5])
        perm, iperm = Metis.permutation(g)
        @test isperm(perm)
        @test isperm(iperm)
        @test length(perm) == 4
        parts = Metis.separator(g)
        @test length(parts) == 4
        @test all(in(1:3), parts)
        # Edge weights are ignored, so the results are those of the unweighted graph
        T = Graphs.smallgraph(:tutte)
        G = SimpleWeightedGraphs.SimpleWeightedGraph(Graphs.nv(T))
        for edge in Graphs.edges(T)
            SimpleWeightedGraphs.add_edge!(G, Graphs.src(edge), Graphs.dst(edge), 0.5)
        end
        @test Metis.permutation(G) == Metis.permutation(T)
        @test Metis.separator(G) == Metis.separator(T)
        # METIS_PartGraph* do use edge weights, so they must still be rejected
        @test_throws ArgumentError Metis.partition(g, 2)
        @test_throws ArgumentError Metis.partition(G, 3)
    end
end
