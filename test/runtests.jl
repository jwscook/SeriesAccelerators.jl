using SeriesAccelerators
using Test, BenchmarkTools

@testset "Shanks" begin
  summand(x, i) = Float64(x^i / factorial(BigInt(i)))
  accelerator = SeriesAccelerators.shanks
  @test exp(0.0) ≈ accelerator(i->summand(0.0, i), 1, 0)[1] rtol=sqrt(eps())
  @test exp(0.0) ≈ accelerator(i->summand(0.0, i), 2, 0)[1] rtol=sqrt(eps())
  @test exp(1.0) ≈ accelerator(i->summand(1.0, i), 10, 0)[1] rtol=sqrt(eps())
  @test exp(1.0) ≈ accelerator(i->summand(1.0, i), 10, 1)[1] rtol=sqrt(eps())
  x = -1.0
  @test exp(x) ≈ accelerator(i->summand(x, i), 10, 0)[1] rtol=sqrt(eps())
  @test exp(x) ≈ accelerator(i->summand(x, i), 10, 0)[1] rtol=sqrt(eps())
  x = -2.0
  @test exp(x) ≈ accelerator(i->summand(x, i), 20, 1)[1] rtol=sqrt(eps())
  x = 2.0
  @test exp(x) ≈ accelerator(i->summand(x, i), 15, 0)[1] rtol=sqrt(eps())
end

@testset "van Wijngaarden" begin
  summand(x, i) = Float64(x^i / factorial(BigInt(i)))
  accelerator = SeriesAccelerators.vanwijngaarden
  for x = [0.0, 1.0, -1.0], c ∈ 1:5, j ∈ 1:5
    i = c + j 
    result = accelerator(i->summand(x, i), i, j)[1]
    #@show x, i, j, exp(x), result
  end

  @test exp(0.0) ≈ accelerator(i->summand(0.0, i), 6, 5)[1] rtol=sqrt(eps())
  @test exp(0.0) ≈ accelerator(i->summand(0.0, i), 7, 5)[1] rtol=sqrt(eps())
  @test exp(1.0) ≈ accelerator(i->summand(1.0, i), 16, 8)[1] rtol=sqrt(eps())
  x = -1.0
  @test exp(x) ≈ accelerator(i->summand(x, i), 16, 8)[1] rtol=sqrt(eps())
  @test exp(x) ≈ accelerator(i->summand(x, i), 32, 16)[1] rtol=sqrt(eps())
  x = -2.0
  @test exp(x) ≈ accelerator(i->summand(x, i), 20, 10)[1] rtol=sqrt(eps())
  @test exp(x) ≈ accelerator(i->summand(x, i), 40, 20)[1] rtol=sqrt(eps())
  x = 2.0
  @test exp(x) ≈ accelerator(i->summand(x, i), 20, 10)[1] rtol=sqrt(eps())
  @test exp(x) ≈ accelerator(i->summand(x, i), 40, 20)[1] rtol=sqrt(eps())
end

# time comparison
summand(x, i) = Float64(x^i / factorial(BigInt(i)))
@btime [SeriesAccelerators.shanks(i->summand(x, i), 20, 2)[1] for x ∈ (-5.0, -2.0, 0.0, 2.0, 5.0)]
@btime [SeriesAccelerators.vanwijngaarden(i->summand(x, i), 20, 10)[1] for x ∈ (-5.0, -2.0, 0.0, 2.0, 5.0)]
for x ∈ range(-10, stop=10, length=21)
  ta = @elapsed a = SeriesAccelerators.shanks(i->summand(x, i), 20, 2)[1]
  tb = @elapsed b = SeriesAccelerators.vanwijngaarden(i->summand(x, i), 20, 10)[1]
  #r = exp(x)
  #@show x, r, (a-r)/r,(b-r)/r, ta / tb
end
