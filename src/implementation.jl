const default_atol = eps()
const default_rtol = sqrt(eps())
const default_initial_iterate = 5
const default_sum_limit = 1_000_000
const default_recursion = 1

function checkconvergence(a::Number, b::Number, atol, rtol)
  return isapprox(a, b, rtol=rtol, atol=atol)
end
function checkconvergence(a, b, atol, rtol)
  return all(i->checkconvergence(a[i], b[i], atol, rtol), eachindex(a, b))
end

function _seriesaccelerator(accelerator::T, series::U,
    recursion::Int, sum_limit::V, rtol, atol
    ) where {T, U, V<:Int, W<:Int}
  iterate = max(default_initial_iterate, recursion + 1)
  old_value = series(iterate)
  isconverged = false
  while !isconverged && iterate < sum_limit
    new_value = accelerator(series, recursion, iterate)
    any(isfinite.(new_value)) || break
    isconverged = checkconvergence(new_value, old_value, atol, rtol)
    old_value = deepcopy(new_value)
    iterate += 1
  end
  return old_value, isconverged
end

# Returns S(n) = sum(series(i) for i in 0:n), cached so each new partial sum
# costs one term evaluation and one addition.
function _memoisedpartialsums(series::T) where {T}
  sums = [series(0)]
  function partialsum(n::Int)
    while length(sums) <= n
      push!(sums, sums[end] + series(length(sums)))
    end
    return sums[n + 1]
  end
  return partialsum, sums
end

"""
    shanks(series, recursion::Int=1, sum_limit::U=1_000_000; rtol=sqrt(eps()), atol=eps()) 

Shanks transformation series accelerator.
https://en.wikipedia.org/wiki/Shanks_transformation

Arguments:
  series (behaves like a Function) : a function that accepts an argument, n::Int,
    and returns the nth value in the series
  sum_limit (optional, Int) : the value at which to stop the series
    (default is 1,000,000)
  rtol (optional, Number) : relative stopping tolerance
  atol (optional, Number) : absolute stopping tolerance
"""
function shanks(series::T,
    recursion::Int=default_recursion,
    sum_limit::U=default_sum_limit;
    rtol=default_rtol, atol=default_atol) where {T, U<:Int}
  @assert 0 <= recursion <= sum_limit (recursion, sum_limit)
  f, _ = _memoisedpartialsums(series)
  sum_limit -= recursion
  return _seriesaccelerator(_shanks, f, recursion, sum_limit, rtol, atol)
end

function _shanks(f::T, recursion::Int, termindex::U) where {T, U<:Int}
  function _shanks_value(An1, An, An_1)
    denominator = ((An1 .- An) .- (An .- An_1))
    any(iszero, denominator) && return An1 # then it's converged
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
  series (behaves like a Function) : a function that accepts an argument, n::Int,
    and returns the nth value in the series
  sum_limit (optional, Int) : the value at which to stop the series
    (default is 1,000,000)
  rtol (optional, Number) : relative stopping tolerance
  atol (optional, Number) : absolute stopping tolerance
"""
function vanwijngaarden(series::T,
    recursion::Int=default_recursion,
    sum_limit::U=default_sum_limit;
    rtol=default_rtol, atol=default_atol) where {T, U<:Int}
  @assert 0 <= recursion recursion
  @assert 0 < sum_limit sum_limit
  partialsum, _ = _memoisedpartialsums(series)
  return _seriesaccelerator(_vanwijngaarden, partialsum, recursion,
    sum_limit, rtol, atol)
end

# Averaging adjacent partial sums, s_{j+1}(n) = (s_j(n) + s_j(n+1)) / 2,
# `recursion` times is Pascal's triangle over S(n), ..., S(n + recursion).
# Collapse it level by level: O(recursion^2) additions, not 2^recursion calls.
function _vanwijngaarden(f::T, recursion::V, n::U
    ) where {T, U<:Int, V<:Int}
  s = [f(n + k) for k in 0:recursion]
  for level in recursion:-1:1, k in 1:level
    s[k] = (s[k] + s[k + 1]) / 2
  end
  return s[1]
end

