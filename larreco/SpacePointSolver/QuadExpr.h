// Christopher Backhouse - bckhouse@fnal.gov

#ifndef QUADEXPR_H
#define QUADEXPR_H

#include <ostream>

class QuadExpr
{
public:
  QuadExpr(double v) : a(0), b(0), c(v) {}
  // Named constructor for independent variable
  static QuadExpr X();

  double Quadratic() const {return a;}
  double Linear()    const {return b;}
  double Constant()  const {return c;}

  double Eval(double x) const;

  QuadExpr& operator+=(const QuadExpr& e);
  QuadExpr& operator-=(const QuadExpr& e);
  QuadExpr& operator*=(const QuadExpr& e);

  QuadExpr operator+(const QuadExpr& e) const;
  QuadExpr operator-(const QuadExpr& e) const;
  QuadExpr operator*(const QuadExpr& e) const;

  QuadExpr operator+(double v) const {return (*this)+QuadExpr(v);}
  QuadExpr operator-(double v) const {return (*this)-QuadExpr(v);}
  QuadExpr operator*(double v) const {return (*this)*QuadExpr(v);}

protected:
  double a, b, c; // a*x^2 + b*x + c
};

inline QuadExpr operator+(double v, const QuadExpr& e){return e+v;}
inline QuadExpr operator-(double v, const QuadExpr& e){return QuadExpr(v)-e;}
inline QuadExpr operator*(double v, const QuadExpr& e){return e*v;}

std::ostream& operator<<(std::ostream&, const QuadExpr&);

#endif
