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
