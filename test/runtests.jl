using Survival
using Base.Test

@testset "Event times" begin
    @testset "EventTime" begin
        @test isevent(EventTime{Int}(44, true))
        @test !isevent(EventTime(3.2, false))
        @test !isevent(EventTime(2.5f0))

        @test iscensored(EventTime(3))
        @test !iscensored(EventTime(2.1, true))

        @test eltype(EventTime(1)) == Int

        @test sprint(show, EventTime(1)) == "1+"

        @test convert(Int, EventTime(1)) == 1
        @test convert(EventTime, 1) == EventTime(1)
    end

    @testset "EventTimeVector" begin
        let x = EventTimeVector{Float32}([3.4f0, 2.1f0], falses(2))
            @test eltype(x) == EventTime{Float32}
            @test length(x) == 2
            @test endof(x) == 2
            @test size(x) == (2,)
            @test eachindex(x) == Base.OneTo(2)
            @test indices(x) == (Base.OneTo(2),)
            @test !isempty(x)
            @test Base.iteratorsize(x) == Base.HasShape()
            @test Base.iteratoreltype(x) == Base.HasEltype()
            @test Base.IndexStyle(x) == Base.IndexLinear()
            @test Base.IndexStyle(typeof(x)) == Base.IndexLinear()
            @test start(x) == 1
            @test next(x, 1) == (EventTime(3.4f0, false), 2)
            @test !done(x, 1)
            @test x[1] == EventTime(3.4f0, false)
            @test x[1:2] == x

            xc = x[:]
            @test xc == x && xc !== x
            cx = copy(x)
            @test cx == x && cx !== x

            x[1] = EventTime(3.5f0, true)
            @test x == EventTimeVector([3.5f0, 2.1f0], [true, false])
            x[2:2] = EventTime(6.7f0)
            @test x == EventTimeVector([3.5f0, 6.7f0], [true, false])

            @test any(isevent, x)
            @test !all(isevent, x)
            @test any(iscensored, x)
            @test !all(iscensored, x)

            @test isevent.(x) == [true, false]
            @test iscensored.(x) == [false, true]

            @test sort(x) == x
        end

        @test_throws DimensionMismatch EventTimeVector([1,2,3], [true])

        @test EventTimeVector([EventTime(1), EventTime(2)]) == EventTimeVector([1,2], falses(2))

        let x = EventTimeVector([1, 2])
            push!(x, EventTime(3))
            @test x == EventTimeVector([1, 2, 3])

            append!(x, EventTimeVector([4, 5]))
            @test x == EventTimeVector([1, 2, 3, 4, 5])

            prepend!(x, EventTimeVector([6]))
            @test x == EventTimeVector([6; 1:5])

            sort!(x)
            @test x == EventTimeVector(collect(1:6))

            y = EventTimeVector([7])
            @test [x; y] == EventTimeVector(collect(1:7))
            @test [y; y; y] == EventTimeVector([7, 7, 7])
        end
    end
end

@testset "Kaplan-Meier" begin
    t = [
        310, 361, 654, 728,  61,  81, 520, 473, 107, 122, 965, 731, 153, 433, 145,  95, 765,
        735,   5, 687, 345, 444,  60, 208, 821, 305, 226, 426, 705, 363, 167, 641, 740, 245,
        588, 166, 559, 450, 529, 351, 201, 524, 199, 550, 551, 543, 293, 511, 511, 371, 201,
         62, 356, 340, 315, 182, 364, 376, 384, 268, 266, 194, 348, 382, 296, 186, 145, 269,
        350, 272, 292, 332, 285, 243, 276,  79, 240, 202, 235, 224, 239, 173, 252,  92, 192,
        211, 175, 203, 105, 177,
    ]
    s = [
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1,
        1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1,
        0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0,
    ]

    # Results generated from R 3.3.2 using survival 2.41
    # survfit(Surv(t, s) ~ 1)
    # $surv
    r_surv = [
        0.9888889, 0.9777778, 0.9666667, 0.9555556, 0.9444444,  0.9333333,  0.9333333,  0.9220884,
        0.9220884, 0.9107045, 0.8993207, 0.8765531, 0.8651693,  0.8537855,  0.8424017,  0.8424017,
        0.8424017, 0.8424017, 0.8305369, 0.8186721, 0.8186721,  0.8066328,  0.7945935,  0.7705149,
        0.7705149, 0.7705149, 0.7580872, 0.7580872, 0.7580872,  0.7452383,  0.7452383,  0.7321639,
        0.7321639, 0.7321639, 0.7186054, 0.7186054, 0.7186054,  0.7045151,  0.7045151,  0.7045151,
        0.7045151, 0.6895254, 0.6895254, 0.6742026, 0.6742026,  0.6585235,  0.6428443,  0.6428443,
        0.6428443, 0.6263611, 0.609878,  0.5933948, 0.5769116,  0.5604284,  0.5604284,  0.5434457,
        0.526463,  0.526463,  0.5089143, 0.5089143, 0.5089143,  0.5089143,  0.4893406,  0.469767,
        0.4501934, 0.4306198, 0.4110461, 0.4110461, 0.3894121,  0.3677781,  0.3677781,  0.3677781,
        0.3432596, 0.3432596, 0.3432596, 0.3432596, 0.3120542,  0.2808487,  0.2496433,  0.2184379,
        0.1872325, 0.1560271, 0.1248217, 0.1248217, 0.08321444, 0.08321444, 0.08321444,
    ]
    # $n.risk
    r_risk = [
        90, 89, 88, 87, 86, 85, 84, 83, 82, 81, 80, 79, 77, 76, 75, 74, 73, 72, 71, 70, 69, 68,
        67, 66, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45,
        44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23,
        22, 21, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1
    ]

    km = fit(KaplanMeier, t, s)
    jl_surv = km.survival

    @test length(jl_surv) == length(r_surv)
    @test all(isapprox.(jl_surv, r_surv, atol=1e-6))

    @test km.times == sort!(unique(t))
    @test km.natrisk == r_risk
    @test km.nevents == [sum(s[t .== tᵢ]) for tᵢ in sort!(unique(t))]
    @test km.ncensor == [sum(iszero, s[t .== tᵢ]) for tᵢ in sort!(unique(t))]

    km_ev = fit(KaplanMeier, EventTimeVector(t, BitVector(s)))
    @test km_ev.times == km.times
    @test km_ev.natrisk == km.natrisk
    @test km_ev.nevents == km.nevents
    @test km_ev.ncensor == km.ncensor
    @test km_ev.survival ≈ jl_surv

    @test_throws DimensionMismatch fit(KaplanMeier, [1, 2], [true])
    @test_throws ArgumentError fit(KaplanMeier, Float64[], Bool[])
end
