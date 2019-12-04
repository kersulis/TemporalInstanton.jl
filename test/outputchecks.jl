@testset "Output checks" begin
    atol = 1e-4
w1 = check_generator_outputs(nd, i, o)
    @test length(w1) == 8
    @test w1[1]["output"] ≈ -0.01774 atol=atol
    w2 = check_decision_variables(i, o)
    @test length(w2) == 3
    @test w2[1]["bus"] == 103

    D = TemporalInstanton.return_demand(i)
    @test D[1][1] ≈ 0.36067 atol=atol
    P = TemporalInstanton.return_injections(nd, o)
    @test P[1][1] ≈ -0.36067 atol=atol
    R = TemporalInstanton.return_decicion_variables(i, o)
    @test R[1][12] ≈ 1.7477 atol=atol
end
