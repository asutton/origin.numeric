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
  fp_natural<4, 10> n1 = 12;
  assert(n1.sig_digits() == 2);
  std::cout << n1 << '\n';

  fp_natural<10, 10> n2 = 12;
  assert(n2.sig_digits() == 2);
  std::cout << n2 << '\n';


  fp_natural<30, 10> n3 = 12345678901234567890ull;
  assert(n3.sig_digits() == 20);
  std::cout << n3 << '\n';
}
