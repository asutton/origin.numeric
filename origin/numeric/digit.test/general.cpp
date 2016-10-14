// Copyright (c) 2008-2016 Andrew Sutton
// All rights reserved

#include <origin/numeric/digit.hpp>

#include <bitset>
#include <vector>
#include <iostream>


using namespace algo;

using digit4 = digit<4>;
using digit10 = digit<10>;
using digit256 = digit<1 << 8>;
using digit16b = digit<radix16>;
using digit24b = digit<1ull << 24>;
using digit32b = digit<radix32>;


template<Digit D>
  requires sizeof(D) == 1
void print_info()
{
  std::cout << "radix: " << (unsigned)radix<D>() << '\n';
  std::cout << "min: " << min<D>() << '\n';
  std::cout << "max: " << max<D>() << '\n';

  using P = promoted_type_t<D>;
  constexpr int N = sizeof(P) * CHAR_BIT;
  std::cout << "lo mask: " << std::bitset<N>(lsd_mask<D>()) << '\n';
  std::cout << "hi mask: " << std::bitset<N>(msd_mask<D>()) << '\n';
  std::cout << "-------------------\n";
}

template<Digit D>
void print_info()
{
  std::cout << "radix: " << radix<D>() << '\n';
  std::cout << "min: " << min<D>() << '\n';
  std::cout << "max: " << max<D>() << '\n';

  using P = promoted_type_t<D>;
  constexpr int N = sizeof(P) * CHAR_BIT;
  std::cout << "lo mask: " << std::bitset<N>(lsd_mask<D>()) << '\n';
  std::cout << "hi mask: " << std::bitset<N>(msd_mask<D>()) << '\n';
  std::cout << "-------------------\n";
}


template<typename T, typename U>
std::ostream&
operator<<(std::ostream& os, std::pair<T, U> const& p) {
  return os << '(' << p.first << ',' << p.second << ')';
}


template<Digit D>
void test_add(D* first, D* limit)
{
  std::cout << "addition\n";
  for (D* i = first; i != limit; ++i) {
    for (D* j = first; j != limit; ++j)
      std::cout << *i << " + " << *j << " == " << add(*i, *j) << '\n';
  }
}


template<Digit D>
void test_sub(D* first, D* limit)
{
  std::cout << "subtraction\n";
  for (D* i = first; i != limit; ++i) {
    for (D* j = first; j != limit; ++j)
      std::cout << *i << " - " << *j << " == " << sub(*i, *j) << '\n';
  }
}


template<Digit D>
void test_mul(D* first, D* limit)
{
  std::cout << "multiplication\n";
  for (D* i = first; i != limit; ++i) {
    for (D* j = first; j != limit; ++j)
      std::cout << *i << " * " << *j << " == " << mul(*i, *j) << '\n';
  }
}


template<Digit D>
void test_div(D* first, D* limit)
{
  std::cout << "division\n";
  for (D* i = first; i != limit; ++i) {
    for (D* j = first; j != limit; ++j)
      if (*j != D(0))
        std::cout << *i << " / " << *j << " == " << div(*i, *j) << '\n';
  }
}


template<Digit D>
void test_twice(D* first, D* limit)
{
  std::cout << "twice\n";
  for (D* i = first; i != limit; ++i) {
    std::cout << *i << " * 2" << " == " << twice(*i, D(0)) << '\n';
  }
}


template<Digit D>
void test_half(D* first, D* limit)
{
  std::cout << "half\n";
  for (D* i = first; i != limit; ++i) {
    std::cout << *i << " / 2" << " == " << half(*i, D(0)) << '\n';
  }
}


template<Digit D, int N>
void 
test_digit(D const (&a)[N])
{
  print_info<D>();
  test_add(a, a + N);
  test_sub(a, a + N);
  test_mul(a, a + N);
  test_div(a, a + N);
  test_twice(a, a + N);
  test_half(a, a + N);  
}


void
test_bit()
{
  bit vals[] = {0, 1};
  test_digit(vals);
}


void
test_digit4()
{
  digit4 vals[] = {0, 1, 2, 3};
  test_digit(vals);
}


void
test_digit10()
{
  using D = digit<10>;
  print_info<D>();
}


int main()
{
  // print_info<digit10>();
  // print_info<digit16b>();
  // print_info<digit24b>();
  // print_info<digit32b>();

  test_bit();
  test_digit4();
  test_digit10();
}
