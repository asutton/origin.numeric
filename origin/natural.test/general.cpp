
#include "natural.hpp"

#include <bitset>
#include <vector>
#include <iostream>


using namespace algo;


template<typename T>
struct fmt_decimal
{
  T num;
};

template<typename T>
fmt_decimal<T> decimal(T const& n) { return fmt_decimal<T>{n}; }


template<typename T>
void print_small(T const& n)
{
  std::uintmax_t s = 0;
  std::uintmax_t e = 1;
  for (int i = 0; i < n.size(); ++i) {
    s += n[i] * e;
    e *= radix<T>();
  }
  std::cout << s << '\n';
}

// Convert a small value into a system integer. Here, "small" means 
// representable within the system int type.
template<typename T>
int
to_int(T const& n)
{
  int s = 0;
  int e = 1;
  for (int i = 0; i < n.size(); ++i) {
    s += n[i] * e;
    e *= radix<T>();
  }
  return s;
}


// Returns the largest r^k for which a / r^k != 0.
template<typename T>
inline T 
largest_scaling(T a, T r)
{
  T b = 1;
  while (a >= r) {
    a /= r;
    b *= r;
  }
  return b;
}


template<typename T>
std::ostream&
operator<<(std::ostream& os, fmt_decimal<T> const& fmt)
{
  T n = fmt.num;
  T x = largest_scaling(n, T(10));
  while (n != T(0)) {
    T q = n / x;
    T r = n % x;
    os << to_int(q);
    n = r;
    x = x / T(10);
  }
  return os;
}


void
test_binary()
{
  using Nat = fp_natural<8, 2>;
  Nat n1 = 25;
  Nat n2 = 3;
  std::cout << n1 << ' ' << n2 << '\n';
  // std::cout << decimal(n1) << ' ' << decimal(n2) << '\n';


  // std::cout << n1 + n2 << ' ' << decimal(n1 + n2) << '\n';

  // {
  //   Nat n = n1 - n2;
  //   std::cout << n << ' ' << decimal(n) << '\n';
  // }
  {
    Nat n = Nat(26) / n2;
    Nat r = Nat(26) % n2;
    std::cout << n << ' ' << decimal(n) << '\n';
    std::cout << n << ' ' << decimal(r) << '\n';
  }

  // std::cout << n1 * n2 << '\n';

  // std::cout << n1 / Nat(2) << '\n';

  // {
  //   Nat n = n2;
  //   std::cout << "init: " << n << '\n';
  //   for (int i = 0; i < 4; ++i) {
  //     twice(n, n);
  //     std::cout << "twice: " << n << '\n';
  //   }

  //   scale_up(n, n2, 4);
  //   std::cout << "scale: " << n << '\n';
  // }


  // {
  //   Nat n = 54;
  //   std::cout << "init: " << n << '\n';
  //   for (int i = 0; i < 4; ++i) {
  //     half(n, n);
  //     std::cout << "half:  " << n << '\n';
  //   }
  //   scale_down(n, Nat(54), 4);
  //   std::cout << "scale: " << n << '\n';
  // }


  // assert(n1 + n1 == Nat(50));
  // std::cout << n1 << " + " << n1 << " == " << (n1 + n1) << '\n';

  // assert(n1 + n2 == Nat(37));
  // std::cout << n1 << " + " << n2 << " == " << (n1 + n2) << '\n';

  // assert(n1 - n2 == Nat(13));
  // std::cout << n1 << " - " << n2 << " == " << (n1 - n2) << '\n';

  // std::cout << n1 << " up == " << (n1 << 1) << '\n';
  // std::cout << n1 << " down == " << (n1 >> 1) << '\n';

  // std::cout << n1 << " * " << n2 << " == " << (n1 * n2) << '\n';
}


void
test_quaternary()
{
  fp_natural<8, 8> n1 = 25;
  std::cout << n1 << '\n';
}


void
test_ap() 
{
  // using Nat = ap_natural<2>;
  // Nat n1 = 5;
  // Nat n2 = 17;
  // std::cout << n1 << " (" << decimal(n1) << ")\n";
  // std::cout << n2 << " (" << decimal(n2) << ")\n";
  // Nat n3 = n1 + n2;
  // std::cout << n3 << " (" << decimal(n3) << ")\n";
}


int main()
{
  using N = fp_natural<8, 2>;
  N big = max<N>();
  std::cout << big << '\n';
  std::cout << decimal(big) << '\n';

  // using D = fp_natural<4, 10>;
  // D n = 12;
  // std::cout << n << '\n';
  // std::cout << twice(n) << '\n';
  // std::cout << half(n) << '\n';
  // std::cout << half(half(n)) << '\n';


  using N2 = fp_natural<64, 32>;
  N2 big2 = max<N2>();
  std::cout << big2 << '\n';
  std::cout << decimal(big2) << '\n';
  

  // N2 zero;
  // std::cout << zero << '\n';
  // std::cout << (zero == zero) << '\n';

  // test_binary();
  // test_ap();
}
