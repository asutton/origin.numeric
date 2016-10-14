// Copyright (c) 2008-2016 Andrew Sutton
// All rights reserved

#ifndef ORIGIN_NUMERIC_FP_NATURAL_HPP
#define ORIGIN_NUMERIC_FP_NATURAL_HPP

#include <origin/numeric/number.hpp>


namespace origin {
namespace numeric {

// -------------------------------------------------------------------------- //
// Fixed precision natural numbers                               [natural.fp] //

// A generic fixed point natural number. This is a fixed-length sequence of 
// P digits with radix R.
//
// Unlike unsigned integers, natural numbers are not modular; overflow and
// underflow result in undefined behavior.
//
// TODO: Factor this class into two parts: one that provides storage and
// access and a second that adds all of the operators. We should define our
// algorithms only in terms of the data representation. Also, add allocator
// support.
//
// TODO: Store a pointer to the most significant digit so that we don't
// end up with unnecessarily expensive computations for relatively small
// numbers stored in greater storage.
// 
// TODO: Support conversions between numbers in other bases? That would make
// sense.
//
// TODO: Factor multiplication and division algorithms into a generic library
// so they can be reused by different numeric facilities.
//
// TODO: An ideal "natural" type would have radix 32 or 64 as the underlying
// digit type, but self-advertise as having radix 2 since the underlying
// value is in fact stored and manipulated in binary.
template<Nonzero P, Radix R>
struct fp_natural
{
  using digit_type = digit<R>;
  using storage_type = std::array<digit_type, P>;
  using iterator = typename storage_type::iterator;
  using const_iterator = typename storage_type::const_iterator;
  using reverse_iterator = typename storage_type::reverse_iterator;
  using const_reverse_iterator = typename storage_type::const_reverse_iterator;

  fp_natural();

  // Copy
  fp_natural(fp_natural const&);
  fp_natural& operator=(fp_natural const&);

  // Move
  fp_natural(fp_natural&&);
  fp_natural& operator=(fp_natural&&);
  
  // Implicit widening.
  template<Nonzero P2>
    requires P2 < P
  fp_natural(fp_natural<P2, R> const&);

  // Value initialization.
  fp_natural(std::uintmax_t n);

  ~fp_natural();

  // Digit access.
  digit_type operator[](std::size_t n) const { return (*digs_)[n]; }
  digit_type& operator[](std::size_t n) { return (*digs_)[n]; }

  // Observers
  constexpr std::size_t size() const { return P; }

  // Iterators
  const_iterator begin() const { return digs_->begin(); }
  const_iterator end() const { return digs_->end(); }
  iterator begin() { return digs_->begin(); }
  iterator end() { return digs_->end(); }
  
  const_reverse_iterator rbegin() const { return digs_->rbegin(); }
  const_reverse_iterator rend() const { return digs_->rend(); }
  reverse_iterator rbegin() { return digs_->rbegin(); }
  reverse_iterator rend() { return digs_->rend(); }
  
  storage_type* digs_;
};

// Initialize this value to 0.
template<Nonzero P, Radix R>
fp_natural<P, R>::fp_natural()
  : digs_(new storage_type())
{
  std::fill(begin(), end(), digit_type(0));
}

// Initialize this value as a copy of n.
template<Nonzero P, Radix R>
fp_natural<P, R>::fp_natural(fp_natural const& n)
  : digs_(new storage_type())
{
  std::copy(n.begin(), n.end(), begin());
}

// Make this value a copy of n.
template<Nonzero P, Radix R>
fp_natural<P, R>&
fp_natural<P, R>::operator=(fp_natural const& n)
{
  if (!digs_)
    digs_ = new storage_type();
  std::copy(n.begin(), n.end(), begin());
  return *this;
}

// Make this value that of n, and reset n.
//
// NOTE: This makes n invalid, but assignable and destructible. No other
// operations are valid.
template<Nonzero P, Radix R>
fp_natural<P, R>::fp_natural(fp_natural&& n)
  : digs_(n.digs_)
{
  n.digs_ = nullptr;
}

// Make this value that of n, and reset n.
//
// TODO: Consider stealing n's pointer.
template<Nonzero P, Radix R>
fp_natural<P, R>&
fp_natural<P, R>::operator=(fp_natural&& n)
{
  std::swap(digs_, n.digs_);
  return *this;
}

// Support conversion from numbers with the same base and fewer digits.
template<Nonzero P, Radix R>
template<Nonzero P2>
  requires P2 < P
fp_natural<P, R>::fp_natural(fp_natural<P2, R> const& n)
  : digs_(new storage_type())
{
  auto iter = std::copy(n.begin(), n.end(), begin()); 
  std::fill(iter, end(), digit_type(0));
}

// Support initialization from standard integer types.
template<Nonzero P, Radix R>
fp_natural<P, R>::fp_natural(std::uintmax_t n)
  : digs_(new storage_type())
{
  n = fill_converted(begin(), end(), n);
  assert(n == 0);
}

// Reclaim memory used by the number.
template<Nonzero P, Radix R>
fp_natural<P, R>::~fp_natural()
{
  delete digs_;
}


// -------------------------------------------------------------------------- //
// Comparison

template<Nonzero P, Radix R>
inline bool
operator==(fp_natural<P, R> const& a, fp_natural<P, R> const& b)
{
  return std::equal(a.begin(), a.end(), b.begin());
}


template<Nonzero P, Radix R>
inline bool
operator!=(fp_natural<P, R> const& a, fp_natural<P, R> const& b)
{
  return !(a == b);
}


template<Nonzero P, Radix R>
inline bool
operator<(fp_natural<P, R> const& a, fp_natural<P, R> const& b)
{
  return less_n(a.end(), b.end(), P);
}


template<Nonzero P, Radix R>
inline bool
operator>(fp_natural<P, R> const& a, fp_natural<P, R> const& b)
{
  return b < a;
}


template<Nonzero P, Radix R>
inline bool
operator<=(fp_natural<P, R> const& a, fp_natural<P, R> const& b)
{
  return !(b < a);
}


template<Nonzero P, Radix R>
inline bool
operator>=(fp_natural<P, R> const& a, fp_natural<P, R> const& b)
{
  return !(a < b);
}

// -------------------------------------------------------------------------- //
// Predicates

// Returns true if the number is odd.
template<Nonzero P, Radix R>
inline bool
is_odd(fp_natural<P, R> const& n)
{
  return is_odd(n[0]);
}


// Return true when the number is even.
template<Nonzero P, Radix R>
inline bool
is_even(fp_natural<P, R> const& n)
{
  return is_even(n[0]);
}

// -------------------------------------------------------------------------- //
// Addition and subtraction

template<Nonzero P, Radix R>
fp_natural<P, R>&
operator+=(fp_natural<P, R>& r, fp_natural<P, R> const& n)
{
  add_n(r.begin(), r.begin(), n.begin(), P);
  return r;
}


template<Nonzero P, Radix R>
fp_natural<P, R>
operator+(fp_natural<P, R> const& a, fp_natural<P, R> const& b)
{
  fp_natural<P, R> r = a;
  return r += b;
}


template<Nonzero P, Radix R>
fp_natural<P, R>&
operator-=(fp_natural<P, R>& r, fp_natural<P, R> const& n)
{
  sub_n(r.begin(), r.begin(), n.begin(), P);
  return r;
}


template<Nonzero P, Radix R>
fp_natural<P, R>
operator-(fp_natural<P, R> const& a, fp_natural<P, R> const& b)
{
  fp_natural<P, R> r = a;
  return r -= b;
}


// -------------------------------------------------------------------------- //
// Exponential scaling

// Compute r = n * 2.
template<Nonzero P, Radix R>
fp_natural<P, R>&
twice(fp_natural<P, R>& r, fp_natural<P, R> const& n)
{
  lsh_n(r.begin(), n.begin(), 1, P);
  return r;
}


// Computes n * 2.
template<Nonzero P, Radix R>
fp_natural<P, R>
twice(fp_natural<P, R> const& n)
{
  fp_natural<P, R> r;
  twice(r, n);
  return r;
}


// Compute r = n / 2.
template<Nonzero P, Radix R>
fp_natural<P, R>&
half(fp_natural<P, R>& r, fp_natural<P, R> const& n)
{
  rsh_n(r.begin(), n.begin(), 1, P);
  return r;
}


// Computes n / 2.
template<Nonzero P, Radix R>
fp_natural<P, R>
half(fp_natural<P, R> const& n)
{
  fp_natural<P, R> r;
  half(r, n);
  return r;
}


// Computes n = n * b^k where b is the base of n.
template<Nonzero P, Radix R>
fp_natural<P, R>&
operator<<=(fp_natural<P, R>& n, int k)
{
  lsh_n(n.begin(), n.begin(), k, P);
  return n;
}


// Computes n * b^k where b is the base of n.
template<Nonzero P, Radix R>
fp_natural<P, R>
operator<<(fp_natural<P, R>& n, int k)
{
  fp_natural<P, R> r = n;
  return r <<= k;
}


// Divide n = n / b^k where b is the base of n.
template<Nonzero P, Radix R>
fp_natural<P, R>&
operator>>=(fp_natural<P, R>& n, int k)
{
  rsh_n(n.begin(), n.begin(), k, P);
  return n;
}


// Computes n / b^k where b is the base of n.
template<Nonzero P, Radix R>
fp_natural<P, R>
operator>>(fp_natural<P, R>& n, int k)
{
  fp_natural<P, R> r = n;
  return r >>= k;
}


namespace detail
{

// Compute r = a * b.
//
// TODO: Rewrite this using digit range algorithms. 
template<Nonzero P1, Nonzero P2, Radix R>
  requires P2 >= 2 * P1
void
long_multiply_non_overflowing(fp_natural<P2, R>& r, 
                              fp_natural<P1, R> const& a, 
                              fp_natural<P1, R> const& b)
{
  for(std::size_t i = 0; i < b.size(); ++i) {
    // Construct the product of a multiplied by the ith digit of b,
    // scaled to the ith power of b. This is really just a careful
    // positioning of values within output vector.
    fp_natural<P2, R> p;
    auto iter = p.begin();
    std::fill(iter, iter + i, 0);
    Digit c = mul(iter + i, a.begin(), a.begin(), b[i]);
    iter += P1;
    *iter++ = c;
    std::fill(iter, iter + P1 - 1, 0);

    // Accumulate the result.
    r += p;
  }
}

// Perform peasant multiplication on two numbers that have been appropriately
// expanded to accommodate the result.
//
// TODO: Rewrite this as an application of the power algorithm in EoP. 
template<Nonzero P, Radix R>
fp_natural<P, R>
peasant_multiply_recursive(fp_natural<P, R> const& a, fp_natural<P, R> const& b)
{
  // requires a,b < 2^(P/2) - 1
  using T = fp_natural<P, R>;
  if (b == T(0))
    return 0;
  if (b == T(1))
    return a;
  fp_natural<P, R> r = peasant_multiply_recursive(twice(a), half(b));
  if (is_odd(b))
    r += a;
  return r;
}

// Peasant multiplication.
//
// Extend the values a and b so that overflow will never happen. The result
// will be stored in number twice as large.
template<Nonzero P, Radix R>
fp_natural<2 * P, R>
peasant_multiply(fp_natural<P, R> const& a, fp_natural<P, R> const& b)
{
  fp_natural<2 * P, R> x = a;
  fp_natural<2 * P, R> y = b;
  return peasant_multiply_recursive(x, y);
}


template<Nonzero P, Radix R>
std::pair<fp_natural<P, R>, fp_natural<P, R>>
slow_divide_recursive(fp_natural<P, R> const& a, fp_natural<P, R> const& b)
{
  using T = fp_natural<P, R>;
  if (a < b)
    return {T(0), a};
  auto qr = slow_divide_recursive(a - b, b);
  return {qr.first + T(1), qr.second};
}


template<Nonzero P, Radix R>
std::pair<fp_natural<P, R>, fp_natural<P, R>>
divide_recursive(fp_natural<P, R> const& a, fp_natural<P, R> const& b)
{
  using T = fp_natural<P, R>;
  assert(b != T(0));
  if (a < b)
    return {T(0), a};
  if (a - b < b)
    return {T(1), a - b};
  auto qr = divide_recursive(a, twice(b));
  T q = twice(qr.first);
  T r = qr.second;
  if (r < b) 
    return {q, r};
  else
    return {q + T(1), r - b};
}


} // namespace detail


template<Nonzero P, Radix R>
fp_natural<P, R>&
operator*=(fp_natural<P, R>& r, fp_natural<P, R> const& n)
{
  // Compute the full product and copy the least digits back into r.
  fp_natural<2 * P, R> x = detail::peasant_multiply(r, n);
  std::copy(x.begin(), x.begin() + P, r.begin());
  return r;
}


template<Nonzero P, Radix R>
fp_natural<P, R>
operator*(fp_natural<P, R> const& a, fp_natural<P, R> const& b)
{
  fp_natural<P, R> r = a;
  return r *= b;
}


template<Nonzero P, Radix R>
fp_natural<P, R>&
operator/=(fp_natural<P, R>& r, fp_natural<P, R> const& n)
{
  auto qr = detail::divide_recursive(r, n);
  return r = qr.first;
}


template<Nonzero P, Radix R>
fp_natural<P, R>
operator/(fp_natural<P, R> const& a, fp_natural<P, R> const& b)
{
  fp_natural<P, R> r = a;
  return r /= b;
}


template<Nonzero P, Radix R>
fp_natural<P, R>&
operator%=(fp_natural<P, R>& r, fp_natural<P, R> const& n)
{
  auto qr = detail::divide_recursive(r, n);
  return r = qr.second;
}


template<Nonzero P, Radix R>
fp_natural<P, R>
operator%(fp_natural<P, R> const& a, fp_natural<P, R> const& b)
{
  fp_natural<P, R> r = a;
  return r %= b;
}

// TODO: Don't print leading non-significant digits?
template<Nonzero P, Radix R>
std::ostream& 
operator<<(std::ostream& os, fp_natural<P, R> const& n)
{
  for (auto iter = n.rbegin(); iter != n.rend(); ++iter)
    os << *iter;
  return os;
}


} // namespace numeric
} // namespace origin


namespace std
{

// Specializations of numeric limits for digits.
//
// TODO: There are a bunch of other properties that could be filled out.
template<::origin::numeric::Nonzero P, ::origin::numeric::Radix R>
struct numeric_limits<::origin::numeric::fp_natural<P, R>>
{
  using T = ::origin::numeric::fp_natural<P, R>;

  static constexpr bool is_specialized = true;
  static constexpr bool is_signed = false;
  static constexpr bool is_integer = true;
  static constexpr bool is_exact = true;
  static constexpr bool is_modular = false;

  // The radix of a bit data type is R.
  static constexpr std::uintmax_t radix = R;

  // A fixed-point number has a precise number of digits.
  static constexpr int digits = P;

  static constexpr T min() { return 0; }
  static constexpr T max() { return -1; } // FIXME: Wrong!
};


} // namespace std


#endif
