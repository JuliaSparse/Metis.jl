using Metis
using Random
using Test
using SparseArrays
import TriangleMesh
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

@testset "Metis.mesh_partition" begin
    poly = TriangleMesh.polygon_Lshape()
    L_mesh = TriangleMesh.create_mesh(poly, info_str="L_mesh", voronoi=true, delaunay=true)
    L_mesh = TriangleMesh.refine(L_mesh, divide_cell_into=1_000, keep_edges=true)
    poly = TriangleMesh.polygon_unitSquareWithHole()
    S_mesh = TriangleMesh.create_mesh(poly, info_str="S_mesh", voronoi=true, delaunay=true)
    S_mesh = TriangleMesh.refine(S_mesh, divide_cell_into=1_000, keep_edges=true)    
    for mesh in (L_mesh, S_mesh), alg in (:DUAL, :NODAL), nparts in (20, 50)
        epart, npart = Metis.mesh_partition(mesh, nparts, alg = alg)
        @test extrema(epart) == (1, nparts)
        @test extrema(npart) == (1, nparts)
        @test all(x -> findfirst(==(x), epart) !== nothing, 1:nparts)
        @test all(x -> findfirst(==(x), npart) !== nothing, 1:nparts)
    end
end