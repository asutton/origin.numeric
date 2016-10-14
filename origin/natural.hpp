// Copyright (c) 2008-2016 Andrew Sutton
// All rights reserved

#ifndef ORIGIN_NUMERIC_NATURAL_HPP
#define ORIGIN_NUMERIC_NATURAL_HPP

// Defines fixed and arbitrary precision natural numbers and their associated
// operations.

#include <origin/digit.hpp>

#include <array>
#include <vector>


namespace origin {
namespace numeric {

// -------------------------------------------------------------------------- //
// Concepts

// True for all unsigned integers that are non-zero.
//
// TODO: Find a better module for this concept.
template<std::size_t N>
concept bool Nonzero()
{
  return N > 0;
}


// -------------------------------------------------------------------------- //
// Digit Range Algorithms                                   [digit.algorithm] //
//
// TODO: Rewrite all of these to use iterators instead of pointers. Or ranges?

// Returns the number of significant digits in a number. If sigdig(f, l) == n,
// then the position of the most significant digit is n - 1.
//
// TODO: This algorithm has the wrong name, and probably does the wrong thing.
// It might be better to simply return an iterator to the most significant
// digit.
template<Digit D>
int
sigdig(D* first, D* limit)
{
  while (limit != first) {
    --limit;
    if (*limit != D(0))
      return limit - first + 1;
  }
  return 0;
}


// Initialize the digits in [first, limit) from the value n. This returns the
// residual value of n computed from the initialization of digits. If the
// number cannot be stored in limit - first digits, then the residue will
// be non-zero.
// 
// If n is -1, then the resulting range will represent the maximum value for 
// those digits.
//
// FIXME: This results spurious overflow if last - first < sizeof(T) * CHAR_BIT 
// and n == -1 (unsigned)
template<Digit D, Unsigned T>
T
fill_converted(D* first, D* limit, T n)
{
  // Handle the max case separately.
  if (n == T(-1)) {
    std::fill(first, limit, max<D>());
    return 0;
  }

  while (first != limit && n != 0) {
    *first++ = n % radix<D>();
    n /= radix<D>();
  }
  std::fill(first, limit, D(0));
  return n;
}


// Add a digit to the number in [first, last).
template<Digit D>
D
add(D* out, D const* first, D const* limit, D n)
{
  D c = n;
  while (first != limit) {
    auto s = add(*first++, D(0), c);
    *out++ = s.first;
    c = s.second;
  }
  return c;
}


// Add the digits in [first1, first1 + n) to those in [first2, first2 + n) and 
// store the results in [out, out + n). Returns the carry digit.
template<Digit D>
D
add_n(D* out, D const* first1, D const* first2, int n)
{
  D c = 0;
  while (n != 0) {
    auto s = add(*first1++, *first2++, c);
    *out++ = s.first;
    c = s.second;
    --n;
  }
  return c;
}


// Subtract the digits in [first2, first2 + n) from those in 
// [first1, first1 + n) and store the results in [out, out + n). Returns a
// borrowed digit.
template<Digit D>
D
sub_n(D* out, D const* first1, D const* first2, int n)
{
  D b = 0;
  while (n != 0) {
    auto s = sub(*first1++, *first2++, b);
    *out++ = s.first;
    b = s.second;
    --n;
  }
  return b;
}


// Multiply each digit in [first, limit) by n. Returns the final carry.
template<Digit D>
D
mul(D* out, D const* first, D const* limit, D n)
{
  D c(0);
  while (first != limit) {
    auto p = mul(*first++, n);
    auto s = add(p.first, c);
    // TODO: I don't believe that the should ever actually carry. Prove it.
    *out++ = s.first;
    assert(s.second == 0); 
    c = p.second + s.second;
  }
  return c;
}


// Shift k digits to the left.
//
// FIXME: This algorithm requires k * n shifts, although k <= digits<D>().
// Still, it would be more effective to determine when digit copies are
// sufficient and when shifts must actually be performed.
//
// TODO: Provide an overflow buffer to accept bits shifted out of the
// representation.
template<Digit D>
void
lsh_n(D* out, D const* in, int k, int n)
{
  assert(k < n);
  while (k != 0) {
    D c = 0;
    while (n != 0) {
      auto x = twice(*in++, c);
      *out++ = x.first;
      c = x.second;
      --n;
    }
    --k;
  }
}


// Shift k digits to the right. This shifts each digit starting from the
// most significant and working toward the least. The remainder in each
// shift is propagated to the next digit.
//
// TODO: Provide an overflow buffer to accept bits shifted out of the
// representation.
template<Digit D>
void
rsh_n(D* out, D const* in, int k, int n)
{
  while (k != 0) {
    D const* limit = in + n;
    D* out_limit = out + n;
    D r = 0;
    while (n != 0) {
      auto x = half(*--limit, r);
      *--out_limit = x.first;
      r = x.second;
      --n;
    }
    --k;
  }
}


// Returns true if the digit sequence in [p - n, p) is lexicographically
// less than that in [q - n, q). Digits are compared from the most significant
// to the least.
template<Digit D>
bool
less_n(D const* p, D const* q, int n)
{
  while (n != 0) {
    --p;
    --q;
    if (*p < *q)
      return true;
    if (*q < *p)
      return false;
    --n;
  }
  return false;
}


// -------------------------------------------------------------------------- //
// Fixed precision natural number

// A generic fixed point natural number. This is a fixed-length sequence of 
// P digits with radix R.
//
// Unlike unsigned integers, natural numbers are not modular; overflow and
// underflow result in undefined behavior.
//
// TODO: Factor this class into two parts: one that provides storage and
// access and a second that adds all of the operators. We should define our
// algorithms only in terms of the data representation.
//
// TODO: Store a pointer to the most significant digit so that we don't
// end up with unnecessarily expensive computations for relatively small
// numbers stored in greater storage.
// 
// TODO: Support conversions between numbers in other bases? That would make
// sense.
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


// -------------------------------------------------------------------------- //
// Multiplication and division

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


// -------------------------------------------------------------------------- //
// Streaming

// TODO: Don't print leading non-significant digits?
template<Nonzero P, Radix R>
std::ostream& 
operator<<(std::ostream& os, fp_natural<P, R> const& n)
{
  for (auto iter = n.rbegin(); iter != n.rend(); ++iter)
    os << *iter;
  return os;
}


// -------------------------------------------------------------------------- //
// Arbitrary precision natural number

// A generic arbitrary precision natural number. This is a variable-length
// sequence of digits in base R.
//
// TODO: Small numbers should not allocate.
template<Radix R>
struct ap_natural
{
  using digit_type = digit<R>;
  using storage_type = std::vector<digit_type>;
  using iterator = digit_type*;
  using const_iterator = digit_type const*;

  ap_natural();
  ap_natural(std::uintmax_t);

  // Digit access.
  digit_type operator[](std::size_t n) const { return digs_[n]; }
  digit_type& operator[](std::size_t n) { return digs_[n]; }

  // Observers
  std::size_t size() const { return digs_.size(); }
  std::size_t digits() const { return size(); }
  std::size_t sig_digits() const { return sigdig(begin(), end()); }

  // Memory
  void extend(std::size_t n) { 
    if (n > size())
      digs_.resize(n, digit_type(0)); 
  }

  // Iterators
  const_iterator begin() const { return digs_.data(); }
  const_iterator end() const { return begin() + size(); }
  iterator begin() { return digs_.data(); }
  iterator end() { return begin() + size(); }
  
  storage_type digs_;
};

template<Radix R>
ap_natural<R>::ap_natural()
  : digs_(1)
{
  digs_[0] = 0;
}


// Support initialization from standard integer types.
template<Radix R>
ap_natural<R>::ap_natural(std::uintmax_t n)
  : digs_(algo::digits(n, R))
{
  fill_converted(begin(), end(), n);
}


// -------------------------------------------------------------------------- //
// Comparison

// Two ap numbers are equal iff the corresponding significant digits are the 
// same.
template<Radix R>
inline bool
operator==(ap_natural<R> const& a, ap_natural<R> const& b)
{
  std::size_t n1 = a.sig_digits();
  std::size_t n2 = a.sig_digits();
  if (n1 != n2)
    return false;
  else
    return std::equal(a.begin(), a.begin() + n1, b.begin());
}


template<Radix R>
inline bool
operator!=(ap_natural<R> const& a, ap_natural<R> const& b)
{
  return !(a == b);
}


// A number with fewer signficant digits is less than any number with more.
// For numbers of the same significance, lexicographically compare their
// digits.
template<Radix R>
inline bool
operator<(ap_natural<R> const& a, ap_natural<R> const& b)
{
  std::size_t n1 = a.sig_digits();
  std::size_t n2 = a.sig_digits();
  if (n1 < n2)
    return true;
  if (n2 < n1)
    return false;
  return less_n(a.end(), b.end(), n1);
}


template<Radix R>
inline bool
operator>(ap_natural<R> const& a, ap_natural<R> const& b)
{
  return b < a;
}


template<Radix R>
inline bool
operator<=(ap_natural<R> const& a, ap_natural<R> const& b)
{
  return !(b < a);
}


template<Radix R>
inline bool
operator>=(ap_natural<R> const& a, ap_natural<R> const& b)
{
  return !(a < b);
}

// -------------------------------------------------------------------------- //
// Predicates

// Returns true if the number is odd.
template<Radix R>
inline bool
is_odd(ap_natural<R> const& n)
{
  return is_odd(n[0]);
}


// Return true when the number is even.
template<Radix R>
inline bool
is_even(ap_natural<R> const& n)
{
  return is_even(n[0]);
}


// -------------------------------------------------------------------------- //
// Addition and subtraction

template<Radix R>
ap_natural<R>&
operator+=(ap_natural<R>& a, ap_natural<R> const& b)
{
  // FIXME: This definitely over-allocates memory. We should only extend
  // the number if the number of significant digits exceeds the size of a.
  std::size_t n = std::max(a.digits(), b.digits());
  a.extend(n + 1);
    Digit c = add_n(a.begin(), a.begin(), b.begin(), b.digits());
  auto iter = a.begin() + b.digits();
  add(iter, iter, a.end(), c);
  return a;
}


template<Radix R>
ap_natural<R>
operator+(ap_natural<R> const& a, ap_natural<R> const& b)
{
  ap_natural<R> r = a;
  return r += b;
}


template<Radix R>
ap_natural<R>&
operator-=(ap_natural<R>& a, ap_natural<R> const& b)
{
  // sub_n(a.begin(), a.begin(), b.begin(), P);
  return a;
}


template<Radix R>
ap_natural<R>
operator-(ap_natural<R> const& a, ap_natural<R> const& b)
{
  ap_natural<R> r = a;
  return r -= b;
}


// -------------------------------------------------------------------------- //
// Streaming

// TODO: Don't print leading non-significant digits?
template<Radix R>
std::ostream& 
operator<<(std::ostream& os, ap_natural<R> const& n)
{
  auto* p = n.end();
  while (p != n.begin()) {
    --p;
    os << *p;
  }
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


// Specializations of numeric limits for digits.
//
// TODO: There are a bunch of other properties that could be filled out.
template<::origin::numeric::Radix R>
struct numeric_limits<::origin::numeric::ap_natural<R>>
{
  using T = ::origin::numeric::ap_natural<R>;

  static constexpr bool is_specialized = true;
  static constexpr bool is_signed = false;
  static constexpr bool is_integer = true;
  static constexpr bool is_exact = true;
  static constexpr bool is_modular = false;

  // The radix of a bit data type is R.
  static constexpr std::uintmax_t radix = R;

  // A fixed-point number has a precise number of digits.
  // static constexpr int digits = P;

  static constexpr T min() { return 0; }
  // static constexpr T max() { return -1; } // FIXME: Wrong!
};



} // namespace std

#endif
