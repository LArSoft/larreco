#include "QuadExpr.h"

#include <cstdlib>
#include <iostream>

QuadExpr QuadExpr::X()
{
  QuadExpr ret(0);
  ret.b = 1;
  return ret;
}

QuadExpr& QuadExpr::operator+=(const QuadExpr& e)
{
  a += e.a;
  b += e.b;
  c += e.c;
  return *this;
}

QuadExpr QuadExpr::operator+(const QuadExpr& e) const
{
  QuadExpr ret = *this;
  ret += e;
  return ret;
}

QuadExpr& QuadExpr::operator-=(const QuadExpr& e)
{
  a -= e.a;
  b -= e.b;
  c -= e.c;
  return *this;
}

QuadExpr QuadExpr::operator-(const QuadExpr& e) const
{
  QuadExpr ret = *this;
  ret -= e;
  return ret;
}

QuadExpr QuadExpr::operator*(const QuadExpr& e) const
{
  if((b != 0 && e.a != 0) ||
     (a != 0 && e.b != 0) ||
     (a != 0 && e.a != 0)){
    std::cout << "(" << *this << ") * (" << e << ")"
              << " does not result in a quadratic expression." << std::endl;
    abort();
  }

  QuadExpr ret(0);
  ret.c = c*e.c;
  ret.b = c*e.b + b*e.c;
  ret.a = c*e.a + a*e.c + b*e.b;

  return ret;
}

QuadExpr& QuadExpr::operator*=(const QuadExpr& e)
{
  *this = *this * e;
  return *this;
}

double QuadExpr::Eval(double x) const
{
  return a*x*x + b*x + c;
}

std::ostream& operator<<(std::ostream& os, const QuadExpr& e)
{
  os << e.Quadratic() << "*x^2 + " << e.Linear() << "*x + " << e.Constant();
  return os;
}
