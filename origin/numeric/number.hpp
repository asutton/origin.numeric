// Copyright (c) 2008-2016 Andrew Sutton
// All rights reserved

#ifndef ORIGIN_NUMERIC_NUMBER_HPP
#define ORIGIN_NUMERIC_NUMBER_HPP

// Common types and algorithms for numeric representations.

#include <origin/numeric/digit.hpp>

#include <array>
#include <vector>


namespace origin {
namespace numeric {

// -------------------------------------------------------------------------- //
// Concepts                                                 [number.concepts] //

// True for all unsigned integers that are non-zero.
//
// TODO: Find a better module for this concept.
template<std::size_t N>
concept bool Nonzero()
{
  return N > 0;
}


} // namespace numeric
} // namespace origin

#endif
