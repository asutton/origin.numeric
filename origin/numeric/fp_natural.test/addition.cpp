// Copyright (c) 2008-2016 Andrew Sutton
// All rights reserved

#include <origin/numeric/fp_natural.hpp>

#include <bitset>
#include <vector>
#include <iostream>


using namespace origin;
using namespace origin::numeric;

template<typename T>
void
check(int a, int b, int c) 
{
  T x = a, y = b, z = c;
  std::cout << x << " + " << y << " == " << z << '\n';
  assert(T(a) + T(b) == T(c));
}


int 
main()
{
  using N = fp_natural<4, 10>;
  check<N>(0, 0, 0);
  check<N>(0, 1, 1);
  check<N>(1, 0, 1);
  
  check<N>(10, 20, 30);
  check<N>(20, 10, 30);

  for (N i = 0; i < N(100); i = i + N(1))
    std::cout << i << ' ';
  std::cout << '\n';

  // FIXME: Check overflow conditions.
  // check<N>(5000, 5000, 1000);
}
