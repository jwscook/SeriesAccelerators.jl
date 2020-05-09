const default_atol = eps()
const default_rtol = sqrt(eps())
const default_initial_iterate = 5
const default_sum_limit = 1_000_000
const default_recursion = 1

function _seriesaccelerator(accelerator::T, series::U,
    recursion::Int, sum_limit::V, rtol, atol
    ) where {T<:Function, U<:Function, V<:Int, W<:Int}
  iterate = max(default_initial_iterate, recursion + 1)
  old_value = series(iterate)
  isconverged = false
  while !isconverged && iterate < sum_limit
    new_value = accelerator(series, recursion, iterate)
    any(isfinite.(new_value)) || break
    isconverged = isapprox(new_value, old_value, rtol=rtol, atol=atol)
    old_value = deepcopy(new_value)
    iterate += 1
  end
  return old_value, isconverged
end

function _memoise(f::T) where {T<:Function}
  data = Dict(0 => f(0))
  function fmemoised(i)
    !haskey(data, i) && (data[i] = f(i))
    return data[i]
  end
  return fmemoised, data
end

"""
    shanks(series, recursion::Int=1, sum_limit::U=1_000_000; rtol=sqrt(eps()), atol=eps()) 

Shanks transformation series accelerator.
https://en.wikipedia.org/wiki/Shanks_transformation

Arguments:
  series (optional, Function) : a function that accepts an argument, n::Int,
    and returns the nth value in the series
  sum_limit (optional, Int) : the value at which to stop the series
    (default is 1,000,000)
  rtol (optional, Number) : relative stopping tolerance
  atol (optional, Number) : absolute stopping tolerance
"""
function shanks(series::T,
    recursion::Int=default_recursion,
    sum_limit::U=default_sum_limit;
    rtol=default_rtol, atol=default_atol) where {T<:Function, U<:Int}
  @assert 0 <= recursion <= sum_limit "$recursion, $sum_limit"
  memoisedseries, data = _memoise(series)
  sum_limit -= recursion
  f(n) = mapreduce(memoisedseries, +, 0:n)
  return _seriesaccelerator(_shanks, f, recursion, sum_limit, rtol, atol)
end

function _shanks(f::T, recursion::Int, termindex::U) where {T<:Function, U<:Int}
  function _shanks_value(An1, An, An_1)
    denominator = ((An1 .- An) .- (An .- An_1))
    any(iszero.(denominator)) && return An1 # then it's converged
    return An1 .- (An1 .- An).^2 ./ denominator
  end
  if recursion == 0
    An1 = f(termindex + 1)
    An = f(termindex + 0)
    An_1 = f(termindex - 1)
    return _shanks_value(An1, An, An_1)
  else
    An1 = _shanks(f, recursion - 1, termindex + 1)
    An = _shanks(f, recursion - 1, termindex)
    An_1 = _shanks(f, recursion - 1, termindex - 1)
    return _shanks_value(An1, An, An_1)
  end
  throw("Shouldn't be able to reach here")
end

"""
    vanwijngaarden(series, recursion::Int=1, sum_limit::U=1_000_000; rtol=sqrt(eps()), atol=eps()) 
van Wijngaarden transformation series accelerator.
https://en.wikipedia.org/wiki/Van_Wijngaarden_transformation

Arguments:
  series (optional, Function) : a function that accepts an argument, n::Int,
    and returns the nth value in the series
  sum_limit (optional, Int) : the value at which to stop the series
    (default is 1,000,000)
  rtol (optional, Number) : relative stopping tolerance
  atol (optional, Number) : absolute stopping tolerance
"""
function vanwijngaarden(series::T,
    recursion::Int=default_recursion,
    sum_limit::U=default_sum_limit;
    rtol=default_rtol, atol=default_atol) where {T<:Function, U<:Int}
  @assert 0 <= recursion "$recursion"
  @assert 0 < sum_limit "$sum_limit"
  memoisedseries, data = _memoise(series)
  return _seriesaccelerator(_vanwijngaarden, memoisedseries, recursion,
    sum_limit, rtol, atol)
end

function _vanwijngaarden(f::T, recursion::V, sum_limit::U
    ) where {T<:Function, U<:Int, V<:Int}
  if recursion == 0
    return mapreduce(n -> f(n), +, 0:sum_limit)
  else
    return 0.5 * mapreduce(k -> _vanwijngaarden(f, recursion-1, k),
      +, [1, -1] .+ sum_limit)
  end
  throw("Shouldn't be able to reach here")
end

