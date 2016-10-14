// Copyright (c) 2008-2016 Andrew Sutton
// All rights reserved

#ifndef ORIGIN_NUMERIC_CONCEPTS_HPP
#define ORIGIN_NUMERIC_CONCEPTS_HPP

// FIXME: This should be replaced by the origin.generic library.

#include <type_traits>
#include <limits>

namespace origin {


template<typename T, typename U>
concept bool Same_as()
{
  return std::is_same<T, U>::value;
}


// True for all types that are standard integral types.
template<typename T>
concept bool Integral()
{
  return std::is_integral<T>::value;
}


// True for all types that are signed integral types.
template<typename T>
concept bool Signed()
{
  return Integral<T>() && std::is_signed<T>::value;
}


// True for all types that are unsigned integral types.
template<typename T>
concept bool Unsigned()
{
  return Integral<T>() && std::is_unsigned<T>::value;
}


// -------------------------------------------------------------------------- //

template<typename T>
struct value_type;

template<typename T>
struct value_type<T*>
{
  using type = T;
};

template<typename T>
struct value_type<T const*>
{
  using type = T;
};

template<typename T>
  requires requires { T::value_type; }
struct value_type<T const*>
{
  using type = typename T::value_type;
};

template<typename T>
using value_type_t = typename value_type<T>::type;


} // namespace origin

#endif
