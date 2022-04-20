@testset "logrank test" begin

    # Simple test with no ties
    km1 = fit(KaplanMeier, [1, 3, 5], [1, 0, 1])
    km2 = fit(KaplanMeier, [2, 4, 6], [1, 1, 0])
    kt1 = logrank_test(km1, km2)
    kt2 = logrank_test([1, 3, 5, 2, 4, 6], [1, 0, 1, 1, 1, 0], [1, 1, 1, 2, 2, 2])
    for kt in [kt1, kt2]
        @test isapprox(kt.stat, 0.073903, atol=1e-5, rtol=1e-5)
        @test isapprox(kt.dof, 1)
        @test isapprox(kt.expected, [1.733333, 2.266667], atol=1e-5, rtol=1e-5)
        @test isapprox(kt.var, [0.9622222 -0.9622222; -0.9622222 0.9622222], atol=1e-5, rtol=1e-5)
    end

    # Simple test with ties
    km1 = fit(KaplanMeier, [1, 2, 1], [1, 0, 1])
    km2 = fit(KaplanMeier, [1, 2, 2], [1, 1, 0])
    kt = logrank_test(km1, km2)
    @test isapprox(kt.stat, 0.04132231, atol=1e-5, rtol=1e-5)
    @test isapprox(kt.dof, 1)
    @test isapprox(kt.expected, [1.833333, 2.1666667], atol=1e-5, rtol=1e-5)
    @test isapprox(kt.var, [0.6722222 -0.6722222; -0.6722222 0.6722222], atol=1e-5, rtol=1e-5)

    # Three group test with different sample sizes
    km1 = fit(KaplanMeier, [1, 2, 1], [1, 0, 1])
    km2 = fit(KaplanMeier, [1, 2, 2, 3], [1, 1, 0, 1])
    km3 = fit(KaplanMeier, [1, 2, 2, 4, 5], [1, 1, 0, 0, 1])
    kt = logrank_test(km1, km2, km3)
    @test isapprox(kt.stat, 1.507214, atol=1e-5, rtol=1e-5)
    @test isapprox(kt.dof, 2)
    @test isapprox(kt.expected, [1.25, 2.41666667, 4.33333], atol=1e-5, rtol=1e-5)
    @test isapprox(kt.var, [0.7329545 -0.3227814 -0.4101732;
                            -0.3227814 1.2704726 -0.9476912;
                            -0.4101732 -0.9476912 1.3578644], atol=1e-5, rtol=1e-5)

    # Test with stratification.
    # Column 1 is time, column 2 is status, column 3 is group,
    # column 4 is strata.
    da = [1 1 1 1; 3 0 1 1; 5 1 1 1; 2 1 2 1; 4 1 2 1; 6 0 2 1
          8 1 1 1; 5 1 1 1; 2 0 1 2; 1 0 2 2; 3 1 2 2; 4 1 2 2]
    km1 = Vector{Vector{KaplanMeier}}()
    for j in 1:2 # stratum
        kx = KaplanMeier[]
        for i in 1:2 # group
            ii = (da[:, 3] .== i) .& (da[:, 4] .== j)
            push!(kx, fit(KaplanMeier, da[ii, 1], da[ii, 2]))
        end
        push!(km1, kx)
    end
    kt1 = logrank_test(km1)
    kt2 = logrank_test(da[:,1], da[:,2], da[:,3], da[:,4])
    for kt in [kt1, kt2]
        @test isapprox(kt.stat, 0.09065547, atol=1e-5, rtol=1e-5)
        @test isapprox(kt.dof, 1)
        @test isapprox(kt.expected, [4.3, 3.7], atol=0.1, rtol=0.1)
    end

    # Test with Wilcoxon weights
    # Stata seems to report the unweighted expected values
    # whereas we are reporing the weighted expected values,
    # so the expected values are not tested here.
    kt = logrank_test(da[:,1], da[:,2], da[:,3]; method=:WBG)
    @test isapprox(kt.stat, 0.51362, atol=1e-5, rtol=1e-5)
    @test isapprox(kt.dof, 1)

    # Test with Wilcoxon weights and strata
    kt = logrank_test(da[:,1], da[:,2], da[:,3], da[:,4]; method=:WBG)
    @test isapprox(kt.stat, 0.10810811, atol=1e-5, rtol=1e-5)
    @test isapprox(kt.dof, 1)

    # Test with Tarone-Ware weights and strata
    kt = logrank_test(da[:,1], da[:,2], da[:,3], da[:,4]; method=:TW)
    @test isapprox(kt.stat, 0.10857866, atol=1e-5, rtol=1e-5)
    @test isapprox(kt.dof, 1)

    # Test with Fleming-Harrington weights
    kt = logrank_test(da[:,1], da[:,2], da[:,3]; method=:FH)
    @test isapprox(kt.stat, 0.9047277, atol=1e-5, rtol=1e-5)
    @test isapprox(kt.dof, 1)

    # Test with Fleming-Harrington weights and strata
    kt = logrank_test(da[:,1], da[:,2], da[:,3], da[:,4]; method=:FH)
    @test isapprox(kt.stat, 0.11739737, atol=1e-5, rtol=1e-5)
    @test isapprox(kt.dof, 1)

    # Another test with Fleming-Harrington weights and strata
    kt = logrank_test(da[:,1], da[:,2], da[:,3], da[:,4]; method=:FH, fh=[0.6, 0.2])
    @test isapprox(kt.stat, 0.63748839, atol=1e-5, rtol=1e-5)
    @test isapprox(kt.dof, 1)

	# Trend test
    # Column 1 is time, column 2 is status, column 3 is group,
    # column 4 is strata.
    da = [1 1 1 1; 3 0 1 1; 5 1 1 1; 2 1 2 1; 4 1 2 1; 6 0 2 1
          8 1 3 1; 5 1 3 1; 2 0 3 2; 1 0 1 2; 3 1 1 2; 4 1 1 2;
          6 1 2 2; 7 0 2 2; 1 1 2 2; 1 1 3 2; 3 0 3 2; 5 0 3 2]
	trend = Dict(1=>1, 2=>2, 3=>3)
    # Trend test without stratification.
    kt = logrank_test(da[:,1], da[:,2], da[:,3]; trend=trend)
	@assert isapprox(kt.stat, 1.93, atol=1e-2, rtol=1e-2)
	@assert isapprox(kt.dof, 1)
    # Trend test with stratification.
    kt = logrank_test(da[:,1], da[:,2], da[:,3], da[:, 4]; trend=trend)
	@assert isapprox(kt.stat, 1.52, atol=1e-2, rtol=1e-2)
	@assert isapprox(kt.dof, 1)
end
