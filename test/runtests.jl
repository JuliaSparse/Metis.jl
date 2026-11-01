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

@testset "Graph structure validation" begin
    using LinearAlgebra: ishermitian, tril, triu
    idx = Metis.LibMetis.idx_t
    # The functions that pass a graph to Metis
    consumers = (
        G -> Metis.partition(G, 2), G -> Metis.partition(G, 2; alg=:RECURSIVE),
        G -> Metis.partition(G, 1), Metis.permutation, Metis.separator,
    )
    # Inputs known to crash or hang Metis go to `check_graph` only, which never calls Metis.
    # The inputs given to `consumers` are small ones that Metis gets through without the check.

    # Directed graphs
    onewaypath = Graphs.SimpleDiGraph(4)
    foreach(i -> Graphs.add_edge!(onewaypath, i, i + 1), 1:3)
    @test Metis.graph(onewaypath).adjncy == [2, 3, 4] # Metis.graph is a raw converter
    for f in consumers
        @test_throws ArgumentError f(onewaypath)
    end
    onewaypathlg = LightGraphs.SimpleDiGraph(4)
    foreach(i -> LightGraphs.add_edge!(onewaypathlg, i, i + 1), 1:3)
    @test_throws ArgumentError Metis.partition(onewaypathlg, 2)
    randomdigraph = Graphs.SimpleDiGraph(2000, 6000; rng=MersenneTwister(1))
    @test_throws ArgumentError Metis.check_graph(Metis.graph(randomdigraph))
    # A directed graph with both directions of each edge works
    T = Graphs.smallgraph(:tutte)
    D = Graphs.SimpleDiGraph(T)
    @test Metis.partition(D, 3) == Metis.partition(T, 3)
    @test Metis.permutation(D) == Metis.permutation(T)
    @test Metis.separator(D) == Metis.separator(T)
    if isdefined(Base, :get_extension)
        import SimpleWeightedGraphs
        @test_throws ArgumentError Metis.partition(SimpleWeightedGraphs.SimpleWeightedDiGraph([1], [2], [1]), 2)
        # Only the structure must be symmetric, not the weights
        W = SimpleWeightedGraphs.SimpleWeightedDiGraph([1, 2, 2, 3], [2, 1, 3, 2], [3, 5, 7, 9])
        @test length(Metis.partition(W, 2)) == 3
    end

    # One-sided stored zeros
    onesided = SparseMatrixCSC(4, 4, [1, 2, 4, 5, 6], Int[2, 1, 3, 2, 1], [1, 1, 1, 1, 0])
    @test ishermitian(onesided)
    @test Metis.graph(onesided).adjncy == [2, 1, 3, 2, 1] # Metis.graph is a raw converter
    for f in consumers
        @test_throws ArgumentError f(onesided)
    end
    dropped = dropzeros(onesided)
    @test Metis.graph(dropped).adjncy == [2, 1, 3, 2]
    @test length(Metis.partition(dropped, 2)) == 4
    @test isperm(Metis.permutation(dropped)[1])
    @test length(Metis.separator(dropped)) == 4
    allzero = sprand(MersenneTwister(1), 2000, 2000, 0.003)
    allzero.nzval .= 0
    @test_throws ArgumentError Metis.check_graph(Metis.graph(allzero))
    @test isperm(Metis.permutation(dropzeros(allzero))[1])

    # Stored zeros on both sides of the diagonal are edges
    A = Graphs.adjacency_matrix(Graphs.grid([20, 20]))
    Z = copy(A)
    Z.nzval .= 0
    @test nnz(Z) == nnz(A)
    @test Metis.graph(Z).adjncy == Metis.graph(A).adjncy
    @test Metis.partition(Z, 2) == Metis.partition(A, 2)
    @test Metis.permutation(Z) == Metis.permutation(A)
    @test Metis.separator(Z) == Metis.separator(A)

    # Triangular input with check_hermitian=false
    P = Graphs.adjacency_matrix(Graphs.path_graph(4))
    for tri in (triu(P), tril(P))
        G = Metis.graph(tri; check_hermitian=false)
        for f in consumers
            @test_throws ArgumentError f(G)
        end
    end
    @test_throws ArgumentError Metis.check_graph(Metis.graph(triu(A); check_hermitian=false))
    @test_throws ArgumentError Metis.partition(Metis.graph(triu(A); check_hermitian=false), 2)
    @test Metis.partition(Symmetric(triu(A)), 2) == Metis.partition(A, 2)
    # check_hermitian=false still allows numerically nonsymmetric values
    N = copy(A)
    N.nzval .= 1:nnz(N)
    @test_throws ArgumentError Metis.graph(N)
    @test Metis.partition(Metis.graph(N; check_hermitian=false), 2) == Metis.partition(A, 2)

    # Malformed Metis.Graph: a path with 3 vertices
    n = idx(3)
    xadj = idx[1, 2, 4, 5]
    adjncy = idx[2, 1, 3, 2]
    @test Metis.check_graph(Metis.Graph(n, xadj, adjncy)) === nothing
    for G in (
        Metis.Graph(n, xadj .- idx(1), adjncy),                 # 0-based xadj
        Metis.Graph(n, idx[1, 2, 4], adjncy),                   # xadj too short
        Metis.Graph(n, idx[1, 2, 4, 5, 5], adjncy),             # xadj too long
        Metis.Graph(n, idx[], adjncy),                          # empty xadj
        Metis.Graph(n, idx[1, 4, 2, 5], adjncy),                # decreasing xadj
        Metis.Graph(idx(-1), xadj, adjncy),                     # negative nvtxs
        Metis.Graph(n, xadj, C_NULL),                           # adjncy not a vector
        Metis.Graph(n, xadj, adjncy[1:3]),                      # adjncy too short
        Metis.Graph(n, idx[1, 2, 4, typemax(idx)], adjncy),     # adjncy far too short
        Metis.Graph(n, xadj, adjncy .- idx(1)),                 # 0-based adjncy
        Metis.Graph(n, xadj, idx[2, 1, 4, 2]),                  # vertex out of range
        Metis.Graph(n, xadj, idx[2, 2, 3, 2]),                  # self edge
        Metis.Graph(n, xadj, adjncy, ones(idx, 2)),             # vwgt too short
        Metis.Graph(n, xadj, adjncy, ones(idx, 4)),             # vwgt too long
        Metis.Graph(n, xadj, adjncy, C_NULL, ones(idx, 3)),     # adjwgt too short
        Metis.Graph(n, idx[1, 3, 4, 4], idx[2, 2, 1]),          # 1 → 2 stored twice, 2 → 1 once
    )
        @test_throws ArgumentError Metis.check_graph(G)
    end
    for f in consumers
        @test_throws ArgumentError f(Metis.Graph(n, xadj, adjncy, ones(idx, 2)))
    end
    # Entries past xadj[end] are not used
    @test Metis.check_graph(Metis.Graph(n, xadj, vcat(adjncy, idx(99)), C_NULL, ones(idx, 5))) === nothing
    # Graphs without edges or vertices are valid
    @test Metis.check_graph(Metis.Graph(n, ones(idx, 4), idx[])) === nothing
    @test Metis.check_graph(Metis.Graph(idx(0), idx[1], idx[])) === nothing
    # adjncy can be C_NULL for a graph without edges
    edgeless = Metis.Graph(n, ones(idx, 4), C_NULL)
    @test length(Metis.partition(edgeless, 2)) == 3
    @test length(Metis.partition(edgeless, 2; alg=:RECURSIVE)) == 3
    @test isperm(Metis.permutation(edgeless)[1])
    @test length(Metis.separator(edgeless)) == 3

    # The error explains the two common causes
    err = try
        Metis.check_graph(Metis.graph(onewaypath))
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("directed graph", err.msg)
    @test occursin("SparseArrays.dropzeros", err.msg)
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
