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

#include <origin/concepts.hpp>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <climits>
#include <algorithm>
#include <limits>
#include <type_traits>
#include <iosfwd>


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


// Returns a mask that will return the bits in the least significant digit of
// the promoted type.
template<Digit D, Integral T = promoted_type_t<D>>
constexpr T
lsd_mask()
{
  return (1ull << bits((T)max<D>())) - 1;
}


template<Digit D, Integral T>
  requires Same_as<T, promoted_type_t<D>>()
constexpr typename D::storage_type
lsd(T n)
{
  return n & lsd_mask<D>();
}


template<Digit D, Integral T = promoted_type_t<D>>
constexpr typename D::storage_type
msd_mask()
{
  return lsd_mask<D>() << bits((T)max<D>());
}

template<Digit D, Integral T>
  requires Same_as<T, promoted_type_t<D>>()
constexpr T
msd(T n)
{
  return (n & msd_mask<D>()) >> ilog2(D::radix);
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



// Computes a + b + c -> (s, c) where s is the sum and c is the carry.
template<Digit D>
std::pair<D, D>
add(D a, D b, D c)
{
  assert(c <= 1);
  auto x = promote(a) + promote(b) + c;
  return {
    lsd<D>(x),     // The sum is the bits in the lsd
    msd<D>(x) != 0 // The carry is 1 if any bits are in the msd
  };
}


// Add two digits, returning the sum and carry.
template<Digit D>
std::pair<D, D>
add(D a, D b)
{
  return add(a, b, D(0));
}


// Computes a - (b + c) -> (d, 0) if b + c < a. Otherwise, computes
// (r + a) - (b + c) -> (d, 1). In both cases d is the difference of
// subtraction.
//
// Note that this is implemented using the Austrian method for subtraction. 
template<Digit D>
std::pair<D, D>
sub(D a, D b, D c)
{
  assert(c <= D(1));
  auto p = promote(a);
  auto q = promote(b) + promote(c);
  if (p >= q) {
    auto x = p - q;
    return { 
      lsd<D>(x), // Result of subtraction
      D(0)       // Did not borrow.
    };
  }
  else {
    auto x = p + radix<D>() - q;
    return { 
      lsd<D>(x), // Result of subtraction
      D(1)       // Did borrow.
    };
  }
}


// Equivalent to sub(a, b, 0).
template<Digit D>
std::pair<D, D>
sub(D a, D b)
{
  return sub(a, b, D(0));
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
