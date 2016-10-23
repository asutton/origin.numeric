// Copyright (c) 2008-2016 Andrew Sutton
// All rights reserved

#ifndef ORIGIN_NUMERIC_DIGIT_HPP
#define ORIGIN_NUMERIC_DIGIT_HPP

// Defines primitive data types and elementary operations for manipulating 
// individual digits within a number. Digits are natural numbers ranging from
// 0 to r - 1 where r is the digit's radix (or base). That value is determined
// by a template argument.
//
// The underlying storage for a digit is a standard integral type whose size 
// is large enough to contain the bits required for each digit. No effort is
// made to compress multiple digits into a single unit of storage.
// 
// This library is not efficient for digits whose base is not m + 1 where m is 
// the maximum value of an unsigned standard integer type.

#include <origin/numeric/concepts.hpp>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <climits>
#include <algorithm>
#include <limits>
#include <type_traits>
#include <iosfwd>

#include <iostream>
#include <bitset>


namespace origin {
namespace numeric {

// -------------------------------------------------------------------------- //
// System limits                                                     [limits] //

// A type is numeric if it specializes std::numeric_limits.
//
// TODO: This belongs in a separate library that defines lots of system
// and type limitations and their related concepts.
template<typename T>
concept bool Numeric() 
{
  return std::numeric_limits<T>::is_specialized;
}


// -------------------------------------------------------------------------- //
// Mathematical support for digits                            [digit.support] //

// Returns floor(log_b n).
template<Integral T>
constexpr T ilog(T b, T n) 
{
  if (n == 1)
    return 0;
  else
    return 1 + ilog(b, T(n / b));
}

// Returns floor(log2 n).
template<Integral T>
constexpr T ilog2(T n) 
{
  return ilog(T(2), n);
}

// Returns the number of digits required to encode the value n in base b.
template<Integral T>
constexpr T
digits(T n, T b)
{
  return ilog(b, n) + 1;
}

// Returns the number of bits needed to store the value n.
template<Integral T>
constexpr T
bits(T n)
{
  return digits(n, T(2));
}


// Returns the minimum digit value for T.
//
// TODO: This belongs in the limits library.
template<Numeric T>
constexpr T 
min() 
{
  return std::numeric_limits<T>::min();
}

// Returns the maximum digit value for T.
//
// TODO: This belongs in the limits library.
template<Numeric T>
constexpr T 
max()
{
  return std::numeric_limits<T>::max();
}

// Returns the base value of the type. This is the radix published by the
// numeric_limits specialization for T.
//
// TODO: This belongs in the limits library.
template<Numeric T>
constexpr std::uintmax_t radix()
{
  return std::numeric_limits<T>::radix;
}


// Returns the number of bits needed to store the value n.
//
// TODO: This does not work for non-digit types since numeric limits does
// not actually define the bits member. We really have a refined version of 
// numeric limits for this type.
template<Numeric T>
constexpr T
bits()
{
  return std::numeric_limits<T>::bits;
}


// -------------------------------------------------------------------------- //
// Digit concepts                                             [digit.concept] //

// The smallest component of a number in some base. A digit is a numeric
// type that is an unsigned integer. 
// 
// The values of all digits d range from 0 <= d < base<D>().
//
// TODO: Digits are numbers. They have arithmetic operators in their basis.
//
// TODO: This concept is a little broadly defined. 
template<typename T>
concept bool Digit() 
{
  return Numeric<T>() 
      && !std::numeric_limits<T>::is_signed
      && std::numeric_limits<T>::is_integer
      && std::numeric_limits<T>::is_exact
      && std::numeric_limits<T>::is_modular
      && std::numeric_limits<T>::digits == 1;
}


// An integer is a radix if it is greater than 2.
template<std::uintmax_t N>
concept bool Radix()
{
  return N >= 2;
}


// -------------------------------------------------------------------------- //
// Digit class                                                  [digit.class] //

namespace detail
{

// Used to select an underlying storage type for the digit. This is the 
// unsigned integer type that is large enough to contain the values of
// the digit: the number of bits it takes to encode the value.
//
// When B is larger than then largest integer size of the host architecture,
// the digit is not representable and the program will be ill-formed
//
// FIXME: When B exceeds system limits, we could construct one of our 
// own fixed-point integers to store the value. But who really needs base-256
// digits.
template<int B>
struct select_storage;

template<int B>
  requires 0 < B && B <= 8
struct select_storage<B> { using type = std::uint8_t; };

template<int B>
  requires 8 < B && B <= 16
struct select_storage<B> { using type = std::uint16_t; };

template<int B>
  requires 16 < B && B <= 32
struct select_storage<B> { using type = std::uint32_t; };

template<int B>
  requires 32 < B && B <= 64
struct select_storage<B> { using type = std::uint64_t; };

template<int B>
  requires 64 < B && B <= 128
struct select_storage<B> { using type = __uint128_t; };

template<Radix B>
using storage_type_t = typename select_storage<B>::type;

} // namespace detail


// The digit class can be used to store digits on a given base. This class
// should not be used for storing digits where 2 <= B <= 64 and B is a power 
// of two.
template<Radix B>
struct digit
{
  // The underlying storage for the digit.
  using storage_type = detail::storage_type_t<ilog2(B)>;

  // The radix of the digit.
  static constexpr std::uintmax_t radix = B;

  constexpr digit() = default;

  // Allow implicit conversion from the storage type.
  constexpr digit(storage_type x)
    : val(x)
  { 
    assert(val < B && "invalid digit"); 
  }

  // Support checked assignment from the storage type.
  constexpr digit& operator=(storage_type x)
  {
    assert(x < B && "invalid digit"); 
    val = x;
    return *this;
  }

  // Also allow implicit conversion to the storage type.
  constexpr operator storage_type() const { return val; }
  
  storage_type val;
};


// Digit steaming
//
// FIXME: We should actually be trying to assign glyphs to these digits.
// Otherwise, we're just printing their ordinal values, which is probably
// okay for large radix digits.
//
// TODO: Assign constraints for C and T.

template<typename C, typename T, Radix B>
std::basic_ostream<C, T>& 
operator<<(std::basic_ostream<C, T>& os, digit<B> n)
{
  return os << n.val;
}

// Specialization for small digits.
template<typename C, typename T, Radix B>
  requires sizeof(digit<B>) == 1
std::basic_ostream<C, T>& 
operator<<(std::basic_ostream<C, T>& os, digit<B> n)
{
  os << (unsigned)n;
}


template<typename C, typename T, Radix B>
std::basic_istream<C, T>& 
operator>>(std::basic_istream<C, T>& is, digit<B> n)
{
  return is >> n.val;
}

template<typename C, typename T, Radix B>
  requires sizeof(digit<B>) == 1
std::basic_istream<C, T>& 
operator>>(std::basic_istream<C, T>& is, digit<B> n)
{
  unsigned x;
  is >> x;
  n = x;
  return is;
}


// Common (or not-so-common) types of digits

// A binary digit
using bit = digit<2>;


// Common large radix values.
constexpr std::uintmax_t radix16 = 1ull << 16;
constexpr std::uintmax_t radix32 = 1ull << 32;



// -------------------------------------------------------------------------- //
// Digit conversions                                          [digit.convert] //
//
// TODO: These are almost certainly internal to digit algorithms and are
// not likely to be useful outside of this module. Make the private?


// Returns the promoted size of the stored value type.
template<Digit D>
using promoted_type_t = detail::storage_type_t<2 * ilog2(D::radix)>;


// Return a fundamental type that can store two the results of addition
// of three digits and the multiplication of two.
template<Digit D>
constexpr typename D::storage_type
native(D d)
{
  return d;
}


// Return a fundamental type that can store two the results of addition
// of three digits and the multiplication of two.
template<Digit D>
constexpr promoted_type_t<D>
promote(D d)
{
  return d.val;
}


// Returns a mask for the least significant digit in a promoted integer.
template<Digit D, Integral T = promoted_type_t<D>>
constexpr T
lsd_mask()
{
  return (1ull << bits((T)max<D>())) - 1;
}


// Returns the least significant digit of a promoted digit value.
//
// TODO: Actually get the constraint on this algorithm right.
template<Digit D, Integral T>
  // requires Same_as<T, promoted_type_t<D>>()
constexpr typename D::storage_type
lsd(T n)
{
  return n & lsd_mask<D>();
}


// Returns a mask for the most significant digit in a promoted integer.
template<Digit D, Integral T = promoted_type_t<D>>
constexpr T
msd_mask()
{
  return lsd_mask<D>() << bits((T)max<D>());
}

// Returns the most significant digit of a promoted digit value.
//
// TODO: Actually get the constraint on this algorithm right.
template<Digit D, Integral T>
  // requires Same_as<T, promoted_type_t<D>>()
constexpr typename D::storage_type
msd(T n)
{
  constexpr T log2 = ilog2(D::radix);
  return (n & msd_mask<D>()) >> log2;
}


// -------------------------------------------------------------------------- //
// Digit algorithms                                         [digit.algorithm] //
//
// TODO: When the digit's radix is m + 1 where m is the max value of an
// unsigned integer type, we can propagate carries and remainders using
// only shifts and we can eliminate masking operations.
//
// TODO: Find a way to avoid promotions when they are not necessary. There's 
// a runtime cost for extending integers. This probably negligible for small
// digits, but when r is 2^32 or 2^64, that can be very expensive.

// Zero-fill the digits in the range [f, l).
template<Digit D>
void 
zero_digits(D* f, D* l)
{
  std::fill(f, l, D(0));
}

// Computes a + b + c -> (s, c) where s is the sum and c is the carry.
template<Digit D>
std::pair<D, D>
add_digits(D a, D b, D c)
{
  assert(c <= 1);
  auto x = promote(a) + promote(b) + c;
  if (x > max<D>())
    return {x - radix<D>(), D(1)};
  else 
    return {x, D(0)};
}

// Add two digits, returning the sum and carry.
template<Digit D>
std::pair<D, D>
add_digits(D a, D b)
{
  return add_digits(a, b, D(0));
}

// Computes a - (b + c) -> (d, 0) if b + c < a. Otherwise, computes
// (r + a) - (b + c) -> (d, 1). In both cases d is the difference of
// subtraction.
//
// Note that this is implemented using the Austrian method for subtraction. 
template<Digit D>
std::pair<D, D>
subtract_digits(D a, D b, D c)
{
  assert(c <= D(1));
  auto p = promote(a);
  auto q = promote(b) + promote(c);
  if (p >= q) {
    auto x = p - q;
    return {x, D(0)};
  }
  else {
    auto x = p + radix<D>() - q;
    return {x, D(1)};
  }
}

// Equivalent to sub(a, b, 0).
template<Digit D>
std::pair<D, D>
subtract_digits(D a, D b)
{
  return subtract_digits(a, b, D(0));
}

// Computes a * b -> (p, c) where p is the least significant digit in the
// produce and c is the overflow digit. For example, multiplying the digits 
// 9 * 9 (base 10) results in the pair (1, 8).
template<Digit D>
std::pair<D, D>
mul(D a, D b)
{
  auto x = promote(a) * promote(b);
  return {
    lsd<D>(x), // The low digit
    msd<D>(x)  // The high digit
  };
}

// Computes a / b -> (q, r) where q is the quotient and r is the remainder
// of division. For example, div(5, 2) yields (2, 1).
template<Digit D>
std::pair<D, D>
div(D a, D b)
{
  return {a / b, a % b};
}

// Computes 2n + c0 -> (p, c1) where c0 an input carry (0 or 1), p is the
// least significant digit in the product) and c1 is the carry (0, 1).
//
// The expression twice(n, 0) is equal to mul(d, 2) but requires fewer 
// operations.
template<Digit D>
std::pair<D, D>
twice(D n, D c) 
{
  auto x = (promote(n) << 1) | promote(c);
  return {
    lsd<D>(x),     // The doubled value
    msd<D>(x) != 0 // The carry contains the shifted-off bit.
  };
}

// Computes (r0 * r) + n / 2 -> (q, r1) where r0 is an input remainder (0 or 1),
// r is the radix of the digit, q is the quotient, and r1 is the remainder of 
// division (0 or 1 since the divisor is 2).
//
// The expression half(n, 0) is equal to div(n, 2) but requires few system
// operations.
template<Digit D>
std::pair<D, D>
half(D n, D r) 
{
  constexpr int shift = ilog2((std::uintmax_t)max<D>());
  auto s = n & 1;
  auto x = (promote(n) + promote(r) * radix<D>()) / 2;
  return {
    x,   // The halved value
    s,   // True only if odd.
  };
}


// -------------------------------------------------------------------------- //
// Predicates

// Returns true if and only if the digit is odd. 
template<Digit D>
inline bool
is_odd(D n) 
{
  return n & 1;
}

// Returns true if and only if the digit is even.
template<Digit D>
inline bool
is_even(D n) 
{
  return !(n & 1);
}


// -------------------------------------------------------------------------- //
// Digit Range Algorithms                                   [digit.algorithm] //
//
// The algorithms in this section are defined on sequences of digits delimited
// by a pair of iterators. 
//
// TODO: Rewrite all of these to use iterators instead of pointers. Or ranges?

// Returns the number of significant digits in a number. If sigdig(f, l) == n,
// then the position of the most significant digit is n - 1.
template<Digit D>
D*
find_most_significant(D* first, D* limit)
{
  while (limit != first) {
    --limit;
    if (*limit != D(0))
      return limit;
  }
  return first;
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

// Add the digit n to the number represented in [first, last), storing
// the output in [out, out + (last - first)). Returns a pair containing
// the iterator out + (last - first) and the carry.
//
// This is equivalent to adding the number in [first, last) to an equally 
// long sequence of zeros, using n as the initial carry.
//
// This operation requires exactly limit - first calls to add_digit().
template<Digit D>
std::pair<D*, D>
add_overflow_digit(D* out, D const* f, D const* l, D n)
{
  std::pair<D, D> result{D(0), n};
  while (f != l) {
    result = add_digits(*f++, D(0), result.second);
    *out++ = result.first;
  }
  return {out, result.second};
}

// Add the significant digits in [fi1, li1) to those in [fi2, li2) and store 
// the results in the sequence of digits [out, out + max(li1 - fi1, li2 - fi2)). 
// Returns a pair containing the iterator past the last computed digit and
// the overflow carry.
//
// TODO: Develop precise aliasing requirements for fi1, fi2, and out.
template<Digit D>
std::pair<D*, D>
add_significant_digits(D* out,
                       D const* fi1, D const* li1, 
                       D const* fi2, D const* li2)
{
  // Add the shared number of significant digits.
  std::pair<D, D> result{D(0), D(0)};
  while (fi1 != li1 && fi2 != li2) {
    result = add_digits(*fi1++, *fi2++, result.second);
    *out++ = result.first;
  }

  // Add the remaining digits in the larger sequence.
  if (fi1 != li1)
    return add_overflow_digit(out, fi1, li1, result.second);
  else if (fi2 != li2)
    return add_overflow_digit(out, fi2, li2, result.second);
  else
    return {out, result.second};
}

// Subtract from the digits in [f, l) an equally long sequence of 0 digits,
// with n being the initial borrow. The subtracted digits are stored in
// [out, out + (f - l)). Returns a pair containing the iterator past the
// last compute out iterator and an overflow borrow digit.
//
// If x is the number in [f, l), this computes x - 0 borrow n.
template<Digit D>
std::pair<D*, D>
subtract_left_overflow_digit(D* out, D const* f, D const* l, D n)
{
  std::pair<D, D> result{D(0), n};
  while (f != l) {
    result = subtract_digits(*f++, D(0), result.second);
    *out++ = result.first;
  }
  return {out, result.second};
}

// Subtract the digits in [f, l) from an equally long sequence of 0 digits
// with n being an initial borrow. The subtracted digits are stored in
// [out, out + (f - l)). Returns a pair containing the iterator past the
// last compute out iterator and an overflow borrow digit.
//
// If x is the number in [f, l), this computes 0 - x borrow n.
template<Digit D>
std::pair<D*, D>
subtract_right_overflow_digit(D* out, D const* f, D const* l, D n)
{
  std::pair<D, D> result{D(0), n};
  while (f != l) {
    result = subtract_digits(D(0), *f++, result.second);
    *out++ = result.first;
  }
  return {out, result.second};
}

// Subtract the significant digits in [fi2, li2) from those in [fi1, li1 + n) 
// and store the results in [out, out max(li1 - fi1, li2 - fi2)). 
template<Digit D>
std::pair<D*, D>
subtract_significant_digits(D* out, 
                            D const* fi1, D const* li1, 
                            D const* fi2, D const* li2)
{
  std::pair<D, D> result{D(0), D(0)};
  while (fi1 != li1 && fi2 != li2) {
    result = subtract_digits(*fi1++, *fi2++, result.second);
    *out++ = result.first;
  }

  if (fi1 != li1)
    return subtract_left_overflow_digit(out, fi1, li1, result.second);
  else if (fi2 != li2)
    return subtract_right_overflow_digit(out, fi2, li2, result.second);
  else
    return {out, result.second};
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

} // namespace numeric
} // namespace origin


namespace std
{

// Specializations of numeric limits for digits.
//
// TODO: There are a bunch of other properties that need to be filled
// out, but they generally only apply to floating point types.
template<::origin::numeric::Radix B>
struct numeric_limits<::origin::numeric::digit<B>>
{
  using T = ::origin::numeric::digit<B>;

  static constexpr bool is_specialized = true;
  static constexpr bool is_signed = false;
  static constexpr bool is_integer = true;
  static constexpr bool is_exact = true;
  static constexpr bool is_bunded = true;
  static constexpr bool is_modular = true;

  // The radix of a bit data type is B.
  static constexpr std::uintmax_t radix = B;

  // It's a digit. There's only one.
  static constexpr int digits = 1;

  // TODO: This is not a n
  static constexpr int bits = ::origin::numeric::bits(radix - 1);
  
  static constexpr T min() { return 0; }
  static constexpr T max() { return B - 1; }
};


} // namespace std

#endif
