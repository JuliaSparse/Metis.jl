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
    # Compare the full exception: a non-square matrix would also fail the symmetry check.
    for check in (true, false)
        tall = sparse([10], [5], [1], 10, 5)
        @test_throws ArgumentError("matrix must be square, got size (10, 5)") Metis.graph(tall; check)
        wide = sparse([1], [9], [1], 3, 10)
        @test_throws ArgumentError("matrix must be square, got size (3, 10)") Metis.graph(wide; check)
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
        # Metis requires edges in both directions
        @test_throws ArgumentError Metis.graph(SimpleWeightedGraphs.SimpleWeightedDiGraph([1], [2], [1]))
        # Issue #56: the weight must be that of the edge from j to i. Metis requires equal
        # weights, so this can only be seen with check=false.
        g = Metis.graph(SimpleWeightedGraphs.SimpleWeightedDiGraph([1, 2], [2, 1], [3, 5]); check=false)
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
        for check in (true, false)
            g = Metis.graph(G; weights=false, check)
            @test g.adjwgt == C_NULL
            @test length(g.adjncy) == 2 * Graphs.ne(GG)
        end
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
        Metis.graph(A; weights=true, check=false)
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
    using LinearAlgebra: ishermitian, triu
    # Directed graphs need both directions of each edge
    onewaypath = Graphs.SimpleDiGraph(4)
    foreach(i -> Graphs.add_edge!(onewaypath, i, i + 1), 1:3)
    @test_throws ArgumentError Metis.graph(onewaypath)
    @test_throws ArgumentError Metis.partition(onewaypath, 2)
    @test_throws ArgumentError Metis.permutation(onewaypath)
    @test_throws ArgumentError Metis.separator(onewaypath)
    err = try
        Metis.graph(onewaypath)
    catch e
        e
    end
    @test occursin("edge 1 → 2 without an edge 2 → 1", err.msg)
    # check=false skips the check, and the graph is converted as is
    @test Metis.graph(onewaypath; check=false).adjncy == [2, 3, 4]
    onewaypathlg = LightGraphs.SimpleDiGraph(4)
    foreach(i -> LightGraphs.add_edge!(onewaypathlg, i, i + 1), 1:3)
    @test_throws ArgumentError Metis.graph(onewaypathlg)
    @test Metis.graph(onewaypathlg; check=false).adjncy == [2, 3, 4]
    # A directed graph with both directions of each edge works
    T = Graphs.smallgraph(:tutte)
    D = Graphs.SimpleDiGraph(T)
    @test Metis.partition(D, 3) == Metis.partition(T, 3)
    @test Metis.permutation(D) == Metis.permutation(T)
    @test Metis.separator(D) == Metis.separator(T)
    if isdefined(Base, :get_extension)
        import SimpleWeightedGraphs
        oneway = SimpleWeightedGraphs.SimpleWeightedDiGraph([1], [2], [1])
        @test_throws ArgumentError Metis.graph(oneway)
        @test Metis.graph(oneway; check=false).adjwgt == [1]
        # permutation and separator ignore the weights, but not the structure
        @test_throws ArgumentError Metis.permutation(oneway)
        @test_throws ArgumentError Metis.separator(oneway)
        # The weights must be equal in both directions
        W = SimpleWeightedGraphs.SimpleWeightedDiGraph([1, 2, 2, 3], [2, 1, 3, 2], [3, 5, 7, 9])
        @test_throws ArgumentError Metis.graph(W)
        @test_throws ArgumentError Metis.partition(W, 2)
        err = try
            Metis.graph(W)
        catch e
            e
        end
        @test occursin("weight 3 for edge 1 → 2 and weight 5 for edge 2 → 1", err.msg)
        @test Metis.graph(W; check=false).adjwgt == [3, 5, 7, 9]
        # The weights aren't used, and aren't checked, by permutation and separator
        @test isperm(Metis.permutation(W)[1])
        @test length(Metis.separator(W)) == 3
        # Equal weights in both directions work
        E = SimpleWeightedGraphs.SimpleWeightedDiGraph([1, 2, 2, 3], [2, 1, 3, 2], [3, 3, 7, 7])
        @test Metis.graph(E).adjwgt == [3, 3, 7, 7]
        @test length(Metis.partition(E, 2)) == 3
    end

    # One-sided stored zeros pass ishermitian, which compares values, but not the pattern check
    onesided = SparseMatrixCSC(4, 4, [1, 2, 4, 5, 6], Int[2, 1, 3, 2, 1], [1, 1, 1, 1, 0])
    @test ishermitian(onesided)
    @test_throws ArgumentError Metis.graph(onesided)
    err = try
        Metis.graph(onesided)
    catch e
        e
    end
    @test occursin("entry (1, 4) is stored and entry (4, 1) isn't", err.msg)
    @test_throws ArgumentError Metis.partition(onesided, 2)
    @test Metis.graph(dropzeros(onesided)).adjncy == [2, 1, 3, 2]
    @test length(Metis.partition(dropzeros(onesided), 2)) == 4
    # check=false skips the checks, and the matrix is converted as is
    @test Metis.graph(onesided; check=false).adjncy == [2, 1, 3, 2, 1]
    # check_hermitian is a deprecated alias of check
    @test (@test_deprecated Metis.graph(onesided; check_hermitian=false)).adjncy == [2, 1, 3, 2, 1]

    # Stored zeros on both sides of the diagonal are edges
    A = Graphs.adjacency_matrix(Graphs.grid([20, 20]))
    Z = copy(A)
    Z.nzval .= 0
    @test Metis.graph(Z).adjncy == Metis.graph(A).adjncy
    @test Metis.partition(Z, 2) == Metis.partition(A, 2)
    @test Metis.permutation(Z) == Metis.permutation(A)
    @test Metis.separator(Z) == Metis.separator(A)

    # Triangular matrices must be wrapped
    @test_throws ArgumentError Metis.graph(triu(A))
    @test Metis.partition(Symmetric(triu(A)), 2) == Metis.partition(A, 2)
    # Without weights only the pattern must be symmetric, not the values
    N = copy(A)
    N.nzval .= 1:nnz(N)
    @test Metis.graph(N).adjncy == Metis.graph(A).adjncy
    @test Metis.partition(N, 2) == Metis.partition(A, 2)
    @test Metis.permutation(N) == Metis.permutation(A)
    @test Metis.separator(N) == Metis.separator(A)
    # With weights the values must be symmetric too
    @test_throws ArgumentError Metis.graph(N; weights=true)
    W = sparse([2, 1], [1, 2], [3, 5])
    err = try
        Metis.graph(W; weights=true)
    catch e
        e
    end
    @test occursin("entry (2, 1) is 3 and entry (1, 2) is 5", err.msg)
    @test Metis.graph(W; weights=true, check=false).adjwgt == [3, 5]
    @test Metis.graph(A + A'; weights=true).adjwgt == fill(2, nnz(A))

    # The check itself
    @test Metis.find_asymmetry(A, true) === nothing
    @test Metis.find_asymmetry(triu(A), false)[3] === :pattern
    @test Metis.find_asymmetry(N, false) === nothing
    @test Metis.find_asymmetry(N, true)[3] === :value
    @test Metis.find_asymmetry(spdiagm(0 => 1:3), true) === nothing # the diagonal is ignored
    @test Metis.find_asymmetry(spzeros(0, 0), true) === nothing
    rng = MersenneTwister(1)
    for _ in 1:50
        S = sprand(rng, 30, 30, 0.1)
        P = SparseMatrixCSC(size(S)..., S.colptr, S.rowval, ones(nnz(S)))
        @test (Metis.find_asymmetry(S, false) === nothing) == (P == permutedims(P))
        @test Metis.find_asymmetry(S + S', true) === nothing
        S2 = S + permutedims(P) # symmetric pattern, values usually not symmetric
        @test Metis.find_asymmetry(S2, false) === nothing
        @test (Metis.find_asymmetry(S2, true) === nothing) == (S2 == permutedims(S2))
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
