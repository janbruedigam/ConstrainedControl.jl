using Test
using SafeTestsets

@safetestset "Example Tests" begin
    include("example_test.jl")
end
