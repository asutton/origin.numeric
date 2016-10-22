// Copyright (c) 2008-2016 Andrew Sutton
// All rights reserved

#include <origin/numeric/fp_natural.hpp>

#include <bitset>
#include <vector>
#include <iostream>


using namespace origin;
using namespace origin::numeric;

int main()
{
  fp_natural<4, 10> n1 = 10;
  fp_natural<4, 10> n2 = 20;
  // std::cout << n1 << ' ' << n2 << '\n';
  assert(n1 == n1);
  assert(n1 != n2);
  assert(n1 < n2);
  assert(n2 > n1);
  assert(n1 <= n1);
  assert(n1 >= n1);

  fp_natural<10, 10> n3 = 10;
  fp_natural<10, 10> n4 = 20;
  assert(n1 == n3);
  assert(n1 != n4);
  assert(n1 < n4);
  assert(n4 > n1);
  assert(n1 <= n3);
  assert(n1 >= n3);
}
