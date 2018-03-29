using Metis
using Compat.Test

@testset "Metis.permutation" begin
    S = sprand(10, 10, 0.5); S = S + S'; fill(S.nzval, 1)
    perm, iperm = Metis.permutation(S)
    @test isperm(perm)
    @test isperm(iperm)
    @test S == S[perm,perm][iperm,iperm]
end
