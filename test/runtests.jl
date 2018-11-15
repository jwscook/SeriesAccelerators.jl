using SeriesAccelerators
using Test, BenchmarkTools

@testset "Shanks" begin
  summand(x, i) = Float64(x^i / factorial(BigInt(i)))
  accelerator = SeriesAccelerators.shanks
  @test exp(0.0) ≈ accelerator(i->summand(0.0, i), 0, 1)[1] rtol=sqrt(eps())
  @test exp(0.0) ≈ accelerator(i->summand(0.0, i), 0, 2)[1] rtol=sqrt(eps())
  @test exp(1.0) ≈ accelerator(i->summand(1.0, i), 0, 10)[1] rtol=sqrt(eps())
  @test exp(1.0) ≈ accelerator(i->summand(1.0, i), 1, 10)[1] rtol=sqrt(eps())
  x = -1.0
  @test exp(x) ≈ accelerator(i->summand(x, i), 0, 10)[1] rtol=sqrt(eps())
  x = -2.0
  @test exp(x) ≈ accelerator(i->summand(x, i), 1, 20)[1] rtol=sqrt(eps())
  x = 2.0
  @test exp(x) ≈ accelerator(i->summand(x, i), 0, 15)[1] rtol=sqrt(eps())
end

@testset "van Wijngaarden" begin
  summand(x, i) = Float64(x^i / factorial(BigInt(i)))
  accelerator = SeriesAccelerators.vanwijngaarden
  for x = [0.0, 1.0, -1.0], c ∈ 1:5, i ∈ 1:5
    j = c + i
    result = accelerator(i->summand(x, i), i, j)[1]
    #@show x, i, j, exp(x), result
  end

  @test exp(0.0) ≈ accelerator(i->summand(0.0, i), 5, 6)[1] rtol=sqrt(eps())
  @test exp(0.0) ≈ accelerator(i->summand(0.0, i), 5, 7)[1] rtol=sqrt(eps())
  @test exp(1.0) ≈ accelerator(i->summand(1.0, i), 8, 16)[1] rtol=sqrt(eps())
  x = -1.0
  @test exp(x) ≈ accelerator(i->summand(x, i), 8, 16)[1] rtol=sqrt(eps())
  @test exp(x) ≈ accelerator(i->summand(x, i), 16, 32)[1] rtol=sqrt(eps())
  x = -2.0
  @test exp(x) ≈ accelerator(i->summand(x, i), 10, 20)[1] rtol=sqrt(eps())
  @test exp(x) ≈ accelerator(i->summand(x, i), 20, 40)[1] rtol=sqrt(eps())
  x = 2.0
  @test exp(x) ≈ accelerator(i->summand(x, i), 10, 20)[1] rtol=sqrt(eps())
  @test exp(x) ≈ accelerator(i->summand(x, i), 20, 40)[1] rtol=sqrt(eps())
end

# time comparison
summand(x, i) = Float64(x^i / factorial(BigInt(i)))
@btime [SeriesAccelerators.shanks(i->summand(x, i), 2, 20)[1] for x ∈ (-5.0, -2.0, 0.0, 2.0, 5.0)]
@btime [SeriesAccelerators.vanwijngaarden(i->summand(x, i), 10, 20)[1] for x ∈ (-5.0, -2.0, 0.0, 2.0, 5.0)]
for x ∈ range(-10, stop=10, length=21)
  ta = @elapsed a = SeriesAccelerators.shanks(i->summand(x, i), 2, 20)[1]
  tb = @elapsed b = SeriesAccelerators.vanwijngaarden(i->summand(x, i), 10, 20)[1]
end


@testset "Vector results" begin
  summand(x, i) = Float64(x^i / factorial(BigInt(i)))
  accelerator = SeriesAccelerators.shanks
  result = accelerator(i->[summand(0.0, i), 2*summand(0.0, i)], 0, 1)[1]
  @test [exp(0.0), 2*exp(0.0)] ≈ result rtol=sqrt(eps())
end
