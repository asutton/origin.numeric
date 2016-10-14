// Copyright (c) 2008-2016 Andrew Sutton
// All rights reserved

#ifndef ORIGIN_NUMERIC_AP_NATURAL_HPP
#define ORIGIN_NUMERIC_AP_NATURAL_HPP

#include <origin/numeric/number.hpp>


namespace origin {
namespace numeric {

// -------------------------------------------------------------------------- //
// Arbitrary precision natural numbers                           [natural.ap] //

// A generic arbitrary precision natural number. This is a variable-length
// sequence of digits in base R.
//
// FIXME: Actually finish implementing me.
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
  : digs_(digits(n, R))
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
