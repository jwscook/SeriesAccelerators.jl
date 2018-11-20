using SeriesAccelerators
using Test, BenchmarkTools
import SpecialFunctions: factorial

@testset "Shanks" begin
  summand(x, i) = Float64(x^i / factorial(BigInt(i)))
  accelerator = SeriesAccelerators.shanks
  @test exp(0.0) ≈ accelerator(i->summand(0.0, i), 0, 1)[1] rtol=sqrt(eps())
  @test exp(0.0) ≈ accelerator(i->summand(0.0, i), 0, 2)[1] rtol=sqrt(eps())
  @test exp(1.0) ≈ accelerator(i->summand(1.0, i), 0, 10)[1] rtol=sqrt(eps())
  @test exp(1.0) ≈ accelerator(i->summand(1.0, i), 1, 10)[1] rtol=sqrt(eps())
  @test exp(1.0) ≈ accelerator(i->summand(1.0, i))[1] rtol=sqrt(eps())
  x = -1.0
  @test exp(x) ≈ accelerator(i->summand(x, i), 0, 10)[1] rtol=sqrt(eps())
  @test exp(x) ≈ accelerator(i->summand(x, i))[1] rtol=sqrt(eps())
  x = -2.0
  @test exp(x) ≈ accelerator(i->summand(x, i), 1, 20)[1] rtol=sqrt(eps())
  @test exp(x) ≈ accelerator(i->summand(x, i))[1] rtol=sqrt(eps())
  x = 2.0
  @test exp(x) ≈ accelerator(i->summand(x, i), 0, 15)[1] rtol=sqrt(eps())
  @test exp(x) ≈ accelerator(i->summand(x, i))[1] rtol=sqrt(eps())
end

@testset "van Wijngaarden" begin
  summand(x, i) = x^i / factorial(Float64(i))
  accelerator = SeriesAccelerators.vanwijngaarden
  for x = [0.0, 1.0, -1.0], c ∈ 1:5, i ∈ 1:5
    j = c + i + 3
    result = accelerator(i->summand(x, i), i, j)[1]
  end
  @test exp(0.0) ≈ accelerator(i->summand(0.0, i), 5, 7)[1] rtol=sqrt(eps())
  @test exp(1.0) ≈ accelerator(i->summand(1.0, i), 8, 16)[1] rtol=sqrt(eps())
  @test exp(1.0) ≈ accelerator(i->summand(1.0, i))[1] rtol=sqrt(eps())
  x = -1.0
  @test exp(x) ≈ accelerator(i->summand(x, i), 10, 25)[1] rtol=sqrt(eps())
  @test exp(x) ≈ accelerator(i->summand(x, i))[1] rtol=sqrt(eps())
  x = -2.0
  @test exp(x) ≈ accelerator(i->summand(x, i), 10, 25)[1] rtol=sqrt(eps())
  @test exp(x) ≈ accelerator(i->summand(x, i))[1] rtol=sqrt(eps())
  x = 2.0
  @test exp(x) ≈ accelerator(i->summand(x, i), 10, 25)[1] rtol=sqrt(eps())
  @test exp(x) ≈ accelerator(i->summand(x, i))[1] rtol=sqrt(eps())
end

@testset "Vector results" begin
  summand(x, i) = Float64(x^i / factorial(BigInt(i)))
  accelerator = SeriesAccelerators.shanks
  result = accelerator(i->[summand(0.0, i), 2*summand(0.0, i)], 0, 1)[1]
  @test [exp(0.0), 2*exp(0.0)] ≈ result rtol=sqrt(eps())
  result = accelerator(i->[summand(0.0, i), 2*summand(0.0, i)])[1]
  @test [exp(0.0), 2*exp(0.0)] ≈ result rtol=sqrt(eps())
end

function mFnnaive(a::AbstractVector{S}, b::AbstractVector{U}, z::V
                 ) where {S<:Number, U<:Number, V<:Number}
  terms = Dict()
  proda(n) = isempty(a) ? 1 : prod(a .+ n)
  prodb(n) = isempty(b) ? 1 : prod(b .+ n)
  function summand(n)
    haskey(terms, n) || (terms[n] = z / (n + 1) * proda(n) / prodb(n))
    return prod(terms[i] for i in 0:n)
  end
  (value, isconverged) = shanks(summand, 8)
  return 1 + value
end


@testset "sum_limit test" begin
a, b, z = [2.36231, -0.901808, 1.00893], [1.39985, -2.34968], 0.970806084927347
  @test mFnnaive(a, b, z) != nothing
end
