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
  fp_natural<4, 10> n1;
  assert(n1.sig_digits() == 0);
  std::cout << n1 << '\n';

  fp_natural<8, 2> n2;
  assert(n1.sig_digits() == 0);
  std::cout << n2 << '\n';
}
