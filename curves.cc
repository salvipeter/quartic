#include "curves.hh"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <optional>

#include <Eigen/LU>

#include "nelder-mead.hh"

extern "C" {
#include "reduction.h"
}

double BezierCurve::bernstein(size_t i, size_t n, double u)
{
  DoubleVector tmp(n + 1, 0.0);
  tmp[n-i] = 1.0;
  double const u1 = 1.0 - u;
  for(size_t k = 1; k <= n; ++k)
    for(size_t j = n; j >= k; --j)
      tmp[j] = tmp[j] * u1 + tmp[j-1] * u;
  return tmp[n];
}

Point BezierCurve::evaluateOneByOne(double u) const
{
  Point p(0.0, 0.0, 0.0);
  for(size_t k = 0; k <= n; ++k)
    p += cp[k] * bernstein(k, n, u);
  return p;
}

void BezierCurve::bernsteinAll(size_t n, double u, DoubleVector &coeff)
{
  coeff.clear(); coeff.reserve(n + 1);
  coeff.push_back(1.0);
  double const u1 = 1.0 - u;
  for(size_t j = 1; j <= n; ++j) {
    double saved = 0.0;
    for(size_t k = 0; k < j; ++k) {
      double const tmp = coeff[k];
      coeff[k] = saved + tmp * u1;
      saved = tmp * u;
    }
    coeff.push_back(saved);
  }
}

Point BezierCurve::evaluate(double u) const
{
  DoubleVector coeff; bernsteinAll(n, u, coeff);
  Point p(0.0, 0.0, 0.0);
  for(size_t k = 0; k <= n; ++k)
    p += cp[k] * coeff[k];
  return p;
}

Point BezierCurve::evaluateWithCachedCofficients(DoubleVector const &coeff) const
{
  Point p(0.0, 0.0, 0.0);
  for(size_t k = 0; k <= n; ++k)
    p += cp[k] * coeff[k];
  return p;
}

Point BezierCurve::evaluateByDeCasteljau(double u) const
{
  PointVector tmp = cp;
  double const u1 = 1.0 - u;
  for(size_t k = 1; k <= n; ++k)
    for(size_t i = 0; i <= n - k; ++i)
      tmp[i] = tmp[i] * u1 + tmp[i+1] * u;
  return tmp[0];
}

void BezierCurve::derivativeControlPoints(size_t d, PointMatrix &dcp) const
{
  dcp.clear(); dcp.resize(d + 1);
  dcp[0] = cp;
  for(size_t k = 1; k <= d; ++k) {
    size_t tmp = n - k + 1;
    dcp[k].reserve(tmp);
    for(size_t i = 0; i <= n - k; ++i)
      dcp[k].push_back((dcp[k-1][i+1] - dcp[k-1][i]) * tmp);
  }
}

void BezierCurve::bernsteinAll(size_t n, double u, DoubleMatrix &coeff)
{
  coeff.clear(); coeff.resize(n + 1);
  coeff[0].push_back(1.0);
  double const u1 = 1.0 - u;
  for(size_t j = 1; j <= n; ++j) {
    coeff[j].reserve(j + 1);
    double saved = 0.0;
    for(size_t k = 0; k < j; ++k) {
      double const tmp = coeff[j-1][k];
      coeff[j].push_back(saved + tmp * u1);
      saved = tmp * u;
    }
    coeff[j].push_back(saved);
  }
}

Point BezierCurve::derivativesByControlPoints(double u, size_t d, VectorVector &der) const
{
  size_t const du = std::min(d, n);
  der.clear(); der.reserve(d + 1);
  DoubleMatrix coeff; bernsteinAll(n, u, coeff);
  PointMatrix dcp; derivativeControlPoints(du, dcp);
  for(size_t k = 0; k <= du; ++k) {
    der.push_back(Vector(0.0, 0.0, 0.0));
    for(size_t j = 0; j <= n - k; ++j)
      der[k] += dcp[k][j] * coeff[n-k][j];
  }
  for(size_t k = n + 1; k <= d; ++k)
    der.push_back(Vector(0.0, 0.0, 0.0));
  return der[0];
}

bool BezierCurve::hasInflections() const
{
  // only meaningful for 2D curves

  // the curve has an inflection, where det(f', f'') = 0, t in (0,1),
  // i.e., we need to find the roots of: f'x * f''y - f'y * f''x,
  // the derivative of which is f'x * f'''y - f'y * f'''x.

  // For now, we just use a simple check for cubic and quartic Bezier curves.
  if (n == 3) {
    Vector n = (cp[2] - cp[1]).unit();
    Vector d1 = cp[3] - cp[2];
    d1 = d1 - n * (d1 * n);
    Vector d2 = cp[0] - cp[1];
    d2 = d2 - n * (d2 * n);
    if (d1 * d2 < 0.0)
      return true;
    else
      return false;
  } else if (n == 4) {
    // ...
  }

  // ... and do a simple search for everything else (TODO)
  const size_t resolution = 1000;
  for (size_t i = 0; i <= resolution; ++i) {
    double u = (double)i / (double)resolution;
    VectorVector der;
    derivativesByControlPoints(u, 2, der);
    if (fabs(der[1].x * der[2].y - der[1].y * der[2].x) < 1.0e-2)
      return true;
  }
  return false;
}

BezierCurve BezierCurve::elevate() const {
  BezierCurve result;
  result.n = n + 1;
  result.cp.push_back(cp.front());
  for (size_t i = 1; i <= n; ++i) {
    double ratio = (double)i / (n + 1);
    result.cp.push_back(cp[i-1] * ratio + cp[i] * (1 - ratio));
  }
  result.cp.push_back(cp.back());
  return result;
}

// As described in
//   T. Várady, P. Salvi, I. Kovács:
//     Enhancement of a multi-sided Bézier surface representation
//       CAGD 55 (2017)
// and earlier in
//   M. Eck:
//     Degree reduction of Bézier curves
//       CAGD 10 (1993)
// and even eariler in
//   A.R. Forrest:
//     Interactive interpolation and approximation by Bézier polynomials
//       The Computer Journal 15 (1972)
BezierCurve BezierCurve::reduce() const {
  BezierCurve result;
  result.n = n - 1;
  size_t m = result.n / 2;
  result.cp.resize(n);
  result.cp[0] = cp.front();
  for (size_t j = 1; j <= m; ++j)
    result.cp[j] = (cp[j] * n - result.cp[j-1] * j) / (n - j);
  auto tmp = result.cp[m];
  result.cp[n-1] = cp.back();
  for (size_t j = n - 1; j >= n - m; --j)
    result.cp[j-1] = (cp[j] * n - result.cp[j] * (n - j)) / j;
  if (m == n - m - 1)           // two different middle control points
    result.cp[m] = (result.cp[m] + tmp) / 2;
  return result;
}

BezierCurve BezierCurve::reduce(size_t target) const {
  BezierCurve result;
  result.n = target;
  std::vector<double> Q((target + 1) * (n + 1));
  reduction_matrix(n, target, 1, 1, &Q[0]);
  for (size_t i = 0; i <= target; ++i) {
    Point p;
    for (size_t j = 0; j <= n; ++j)
      p += cp[j] * Q[i * (n + 1) + j];
    result.cp.push_back(p);
  }
  return result;
}

BezierCurve BezierCurve::interpolateUniform(const PointVector &points) {
  BezierCurve result;
  result.n = points.size() - 1;
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(points.size(), points.size());
  Eigen::MatrixXd b = Eigen::MatrixXd::Zero(points.size(), 3);
  for (size_t i = 0; i < points.size(); ++i) {
    double u = (double)i / result.n;
    DoubleVector coeff;
    bernsteinAll(result.n, u, coeff);
    for (size_t j = 0; j <= result.n; ++j)
      A(i, j) = coeff[j];
    b(i, 0) = points[i].x;
    b(i, 1) = points[i].y;
    b(i, 2) = points[i].z;
  }
  Eigen::MatrixXd x = A.fullPivLu().solve(b);
  for(size_t i = 0; i <= result.n; ++i)
    result.cp.push_back(Point(x(i, 0), x(i, 1), x(i, 2)));
  return result;
}

size_t BSplineCurve::findSpan(double u) const
{
  if(u == knots[n+1])
    return n;
  return (std::upper_bound(knots.begin() + p + 1, knots.end(), u) - knots.begin()) - 1;
}

void BSplineCurve::basisFunctions(size_t i, double u, DoubleVector &coeff) const
{
  coeff.clear(); coeff.reserve(p + 1);
  coeff.push_back(1.0);
  DoubleVector left(p + 1), right(p + 1);
  for(size_t j = 1; j <= p; ++j) {
    left[j]  = u - knots[i+1-j];
    right[j] = knots[i+j] - u;
    double saved = 0.0;
    for(size_t r = 0; r < j; ++r) {
      double tmp = coeff[r] / (right[r+1] + left[j-r]);
      coeff[r] = saved + tmp * right[r+1];
      saved = tmp * left[j-r];
    }
    coeff.push_back(saved);
  }
}

Point BSplineCurve::evaluate(double u) const
{
  size_t span = findSpan(u);
  DoubleVector coeff; basisFunctions(span, u, coeff);
  Point point(0.0, 0.0, 0.0);
  for(size_t i = 0; i <= p; ++i)
    point += cp[span - p + i] * coeff[i];
  return point;
}

void BSplineCurve::basisFunctionDerivatives(size_t i, double u, size_t d, DoubleMatrix &der) const
{
  der.clear(); der.resize(d + 1);
  DoubleVector left(p + 1), right(p + 1), a[2];
  a[0].resize(p + 1); a[1].resize(p + 1);
  DoubleMatrix ndu(p + 1);
  ndu[0].resize(p + 1); ndu[0][0] = 1.0;
  for(size_t j = 1; j <= p; ++j) {
    ndu[j].resize(p + 1);
    left[j] = u - knots[i+1-j];
    right[j] = knots[i+j] - u;
    double saved = 0.0;
    for(size_t r = 0; r < j; ++r) {
      // lower triangle
      ndu[j][r] = right[r+1] + left[j-r];
      double const tmp = ndu[r][j-1] / ndu[j][r];
      // upper triangle
      ndu[r][j] = saved + tmp * right[r+1];
      saved = tmp * left[j-r];
    }
    ndu[j][j] = saved;
  }
  for(size_t j = 0; j <= p; ++j)
    der[0].push_back(ndu[j][p]);
  for(size_t r = 0; r <= p; ++r) {
    size_t s1 = 0, s2 = 1;
    a[0][0] = 1.0;
    for(size_t k = 1; k <= d; ++k) {
      double dd = 0.0;
      int const rk = r - k;
      int const pk = p - k;
      if(r >= k) {
        a[s2][0] = a[s1][0] / ndu[pk+1][rk];
        dd = a[s2][0] * ndu[rk][pk];
      }
      size_t const j1 = rk >= -1 ? 1 : -rk;
      size_t const j2 = (int)r - 1 <= pk ? k - 1 : p - r;
      for(size_t j = j1; j <= j2; ++j) {
        a[s2][j] = (a[s1][j] - a[s1][j-1]) / ndu[pk+1][rk+j];
        dd += a[s2][j] * ndu[rk+j][pk];
      }
      if(r <= (size_t)pk) {
        a[s2][k] = -a[s1][k-1] / ndu[pk+1][r];
        dd += a[s2][k] * ndu[r][pk];
      }
      der[k].push_back(dd);
      std::swap(s1, s2);
    }
  }
  size_t r = p;
  for(size_t k = 1; k <= d; ++k) {
    for(size_t j = 0; j <= p; ++j)
      der[k][j] *= r;
    r *= p - k;
  }
}

Point BSplineCurve::derivatives(double u, size_t d, VectorVector &der) const
{
  size_t const du = std::min(d, p);
  der.clear();
  size_t span = findSpan(u);
  DoubleMatrix nder; basisFunctionDerivatives(span, u, du, nder);
  for(size_t k = 0; k <= du; ++k) {
    der.push_back(Vector(0.0, 0.0, 0.0));
    for(size_t j = 0; j <= p; ++j)
      der[k] += cp[span-p+j] * nder[k][j];
  }
  for(size_t k = p + 1; k <= d; ++k)
    der.push_back(Vector(0.0, 0.0, 0.0));
  return der[0];
}

void BSplineCurve::derivativeControlPoints(size_t d, size_t r1, size_t r2, PointMatrix &dcp) const
{
  dcp.clear(); dcp.resize(d + 1);
  size_t r = r2 - r1;
  dcp[0].reserve(r + 1);
  for(size_t i = 0; i <= r; ++i)
    dcp[0].push_back(cp[r1+i]);
  for(size_t k = 1; k <= d; ++k) {
    dcp[k].reserve(r + 1 - k);
    size_t tmp = p - k + 1;
    for(size_t i = 0; i <= r - k; ++i)
      dcp[k].push_back((dcp[k-1][i+1] - dcp[k-1][i]) * tmp / (knots[r1+i+p+1] - knots[r1+i+k]));
  }
}

void BSplineCurve::basisFunctionsAll(size_t i, double u, DoubleMatrix &coeff) const
{
  coeff.clear(); coeff.resize(p + 1);
  coeff[0].push_back(1.0);
  DoubleVector left(p + 1), right(p + 1);
  for(size_t j = 1; j <= p; ++j) {
    coeff[j].reserve(j + 1);
    left[j]  = u - knots[i+1-j];
    right[j] = knots[i+j] - u;
    double saved = 0.0;
    for(size_t r = 0; r < j; ++r) {
      double tmp = coeff[j-1][r] / (right[r+1] + left[j-r]);
      coeff[j].push_back(saved + tmp * right[r+1]);
      saved = tmp * left[j-r];
    }
    coeff[j].push_back(saved);
  }
}

Point BSplineCurve::derivativesByControlPoints(double u, size_t d, VectorVector &der) const
{
  size_t const du = std::min(d, p);
  der.clear();
  size_t span = findSpan(u);
  DoubleMatrix coeff; basisFunctionsAll(span, u, coeff);
  PointMatrix dcp; derivativeControlPoints(du, span - p, span, dcp);
  for(size_t k = 0; k <= du; ++k) {
    der.push_back(Vector(0.0, 0.0, 0.0));
    for(size_t j = 0; j <= p - k; ++j)
      der[k] += dcp[k][j] * coeff[p-k][j];
  }
  for(size_t k = p + 1; k <= d; ++k)
    der.push_back(Vector(0.0, 0.0, 0.0));
  return der[0];
}

size_t BSplineCurve::findSpanWithMultiplicity(double u, size_t &multi) const
{
  auto range = std::equal_range(knots.begin(), knots.end(), u);
  multi = range.second - range.first;

  if(u == knots[n+1])
    return n;
  return (range.second - knots.begin()) - 1;
}

Point BSplineCurve::evaluateByKnotInsertion(double u) const
{
  if(u == knots[0])
    return cp[0];
  if(u == knots[n+p+1])
    return cp[n];
  size_t s, k = findSpanWithMultiplicity(u, s), r = p - s;
  PointVector tmp; tmp.reserve(r + 1);
  std::copy_n(cp.begin() + k - p, r + 1, std::back_inserter(tmp));
  for(size_t j = 1; j <= r; ++j)
    for(size_t i = 0; i <= r - j; ++i) {
      double alpha = (u - knots[k-p+j+i]) / (knots[i+k+1] - knots[k-p+j+i]);
      tmp[i] = tmp[i+1] * alpha + tmp[i] * (1.0 - alpha);
    }
  return tmp[0];
}

Point BSplineCurve::evaluate2DRational(double u) const
{
  Point p = evaluate(u);
  return Point(p.x / p.z, p.y / p.z, 1.0);
}

size_t BSplineCurve::binomial(size_t n, size_t k)
{
  if(k > n)
    return 0;
  size_t result = 1;
  for(size_t d = 1; d <= k; ++d, --n)
    result = result * n / d;
  return result;
}

Point BSplineCurve::derivatives2DRational(double u, size_t d, VectorVector &der) const
{
  der.clear(); der.reserve(d + 1);
  VectorVector der3d; derivativesByControlPoints(u, d, der3d);
  for(size_t k = 0; k <= d; ++k) {
    Vector v = der3d[k];
    for(size_t i = 1; i <= k; ++i)
      v = v - der[k-i] * der3d[i].z * binomial(k, i);
    der.push_back(v / der3d[0].z);
  }
  return der[0];
}

BSplineCurve BSplineCurve::insertKnot(double u, size_t k, size_t s, size_t r) const
{
  PointVector tmp; tmp.reserve(p - s + 1);

  BSplineCurve result;
  result.p = p; result.n = n + r;

  result.knots.reserve(knots.size() + r);
  std::copy_n(knots.begin(), k + 1, std::back_inserter(result.knots));
  std::fill_n(std::back_inserter(result.knots), r, u);
  std::copy(knots.begin() + k + 1, knots.end(), std::back_inserter(result.knots));

  result.cp.resize(cp.size() + r);
  std::copy_n(cp.begin(), k - p + 1, result.cp.begin());
  std::copy(cp.begin() + k - s, cp.end(), result.cp.begin() + r + k - s);

  std::copy_n(cp.begin() + k - p, p - s + 1, std::back_inserter(tmp));
  size_t L;
  for(size_t j = 1; j <= r; ++j) {
    L = k - p + j;
    for(size_t i = 0; i <= p - j - s; ++i) {
      double alpha = (u - knots[L+i]) / (knots[i+k+1] - knots[L+i]);
      tmp[i] = tmp[i+1] * alpha + tmp[i] * (1.0 - alpha);
    }
    result.cp[L] = tmp[0];
    result.cp[k+r-j-s] = tmp[p-j-s];
  }
  if (p > s + r + 1)
    std::copy_n(tmp.begin() + 1, p - s - 1 - r, result.cp.begin() + L);

  return result;
}

BSplineCurve BSplineCurve::refineKnots(DoubleVector new_knots) const
{
  size_t const r = new_knots.size();
  size_t const a = findSpan(new_knots[0]);
  size_t const b = findSpan(new_knots[r-1]) + 1;

  BSplineCurve result;
  result.p = p; result.n = n + r;
  result.knots.resize(knots.size() + r);
  result.cp.resize(cp.size() + r);

  std::copy_n(knots.begin(), a + 1, result.knots.begin());
  std::copy(knots.begin() + b + p, knots.end(), result.knots.begin() + b + p + r);
  std::copy_n(cp.begin(), a - p + 1, result.cp.begin());
  std::copy(cp.begin() + b - 1, cp.end(), result.cp.begin() + b - 1 + r);

  size_t i = b + p - 1,         // next position in knots
         k = b + p + r,         // next position in result.knots
         j = r;                 // next position in new_knots
  do {
    --j; --k;
    for(; new_knots[j] <= knots[i] && i > a; --i, --k) {
      result.knots[k] = knots[i];
      result.cp[k-p-1] = cp[i-p-1];
    }
    result.cp[k-p-1] = result.cp[k-p];
    for(size_t l = 1; l <= p; ++l) {
      size_t const index = k - p + l;
      double alpha = result.knots[k+l] - new_knots[j];
      if(fabs(alpha) == 0.0)
        result.cp[index-1] = result.cp[index];
      else {
        alpha /= result.knots[k+l] - knots[i-p+l];
        result.cp[index-1] = result.cp[index-1] * alpha +
                             result.cp[index] * (1.0 - alpha);
      }
    }
    result.knots[k] = new_knots[j];
  } while(j > 0);
  return result;
}

Point BSplineCurve::projectPoint(Point const &point, double &u, double &distance,
                                 size_t resolution, double distance_tol, double cosine_tol) const
{
  // If we know that point is on the curve,
  // we could check only the spans of the convex hull point is in
  double const span_min = knots[p], span_max = knots[n+1];
  distance = std::numeric_limits<double>::max();
  for(size_t i = 0; i < resolution; ++i) {
    double const v = span_min + (span_max - span_min) * (double)i / (double)(resolution - 1);
    double const d = (evaluate(v) - point).norm();
    if(d < distance) {
      distance = d;
      u = v;
    }
  }

  VectorVector der;
  derivativesByControlPoints(u, 2, der);
  Vector deviation = der[0] - point;

  while(distance > distance_tol) {
    double const scaled_error = der[1] * deviation;
    double const cosine_err = fabs(scaled_error) / (der[1].norm() * distance);
    if(cosine_err < cosine_tol)
      break;
    
    double const old = u;
    u -= scaled_error / (der[2] * deviation + der[1] * der[1]);
    u = std::min(std::max(u, span_min), span_max);

    if((der[1] * (u - old)).norm() < distance_tol)
      break;

    derivativesByControlPoints(u, 2, der);
    deviation = der[0] - point;
    distance = deviation.norm();
  }

  return der[0];
}

std::list<BezierCurve> BSplineCurve::convertToBezierCurves() const
{
  BSplineCurve bsp = *this;
  std::list<BezierCurve> result;
  size_t s;
  for (DoubleVector::const_iterator i = knots.begin(), ie = knots.end(); i != ie; ++i) {
    size_t k = bsp.findSpanWithMultiplicity(*i, s);
    if (s < p)
      bsp = bsp.insertKnot(*i, k, s, p - s);
  }
  double lastKnot = bsp.knots.front();
  for (size_t i = 0, ie = bsp.knots.size(); i != ie; ++i) {
    if (bsp.knots[i] == lastKnot)
      continue;
    lastKnot = bsp.knots[i];
    BezierCurve bc;
    bc.n = p;
    for (size_t j = 0; j <= p; ++j)
      bc.cp.push_back(bsp.cp[i - p - 1 + j]);
    result.push_back(bc);
  }
  return result;
}

BSplineCurve BSplineCurve::interpolate(PointVector const &points, size_t p, bool centripetal)
{
  size_t const n = points.size() - 1;
  size_t const m = n + p + 1;
  DoubleVector u; u.reserve(n + 1);
  u.push_back(0.0);
  double total = 0.0;
  for(size_t i = 1; i <= n; ++i) {
    if (centripetal)
      u.push_back(std::sqrt((points[i] - points[i-1]).norm()));
    else
      u.push_back((points[i] - points[i-1]).norm());
    total += u.back();
  }
  for(size_t i = 1; i < n; ++i)
    u[i] = u[i-1] + u[i] / total;
  u.back() = 1.0;

  BSplineCurve result;
  result.p = p; result.n = n;
  result.knots.reserve(m + 1);
  result.cp.reserve(n + 1);

  std::fill_n(std::back_inserter(result.knots), p + 1, 0.0);
  for(size_t j = 1; j <= n - p; ++j) {
    double t = 0.0;
    for(size_t i = j; i <= j + p - 1; ++i)
      t += u[i];
    result.knots.push_back(t / p);
  }
  std::fill_n(std::back_inserter(result.knots), p + 1, 1.0);

  Eigen::MatrixXd A(n + 1, n + 1); A.setZero();
  for(size_t i = 0; i <= n; ++i) {
    size_t span = result.findSpan(u[i]);
    DoubleVector coeff; result.basisFunctions(span, u[i], coeff);
    for(size_t j = 0; j <= p; ++j)
      A(i, span - p + j) = coeff[j];
  }
  Eigen::MatrixXd b(n + 1, 3);
  for(size_t i = 0; i <= n; ++i) {
    b(i, 0) = points[i].x;
    b(i, 1) = points[i].y;
    b(i, 2) = points[i].z;
  }
  Eigen::MatrixXd x = A.fullPivLu().solve(b);
  for(size_t i = 0; i <= n; ++i)
    result.cp.push_back(Point(x(i, 0), x(i, 1), x(i, 2)));

  return result;
}

// Fits a parabola on the given points, and returns the tangent at u1
Vector BSplineCurve::parabolaRule(double u1, Point const &p1,
                                  double u2, Point const &p2,
                                  double u3, Point const &p3) {
  Eigen::MatrixXd A(3, 3);
  for (size_t i = 0; i < 3; ++i) {
    A(0, i) = BezierCurve::bernstein(i, 2, u1);
    A(1, i) = BezierCurve::bernstein(i, 2, u2);
    A(2, i) = BezierCurve::bernstein(i, 2, u3);
  }
  Eigen::MatrixXd b(3, 3);
  b(0, 0) = p1.x; b(0, 1) = p1.y; b(0, 2) = p1.z;
  b(1, 0) = p2.x; b(1, 1) = p2.y; b(1, 2) = p2.z;
  b(2, 0) = p3.x; b(2, 1) = p3.y; b(2, 2) = p3.z;
  Eigen::MatrixXd x = A.fullPivLu().solve(b);
  BezierCurve curve;
  curve.n = 2;
  for (size_t i = 0; i < 3; ++i)
    curve.cp.push_back(Point(x(i, 0), x(i, 1), x(i, 2)));
  VectorVector der;
  curve.derivativesByControlPoints(u1, 1, der);
  return der[1];
}

BSplineCurve BSplineCurve::interpolateWithParabolaRule(PointVector const &points, size_t p,
                                                       bool centripetal)
{
  size_t const n = points.size() - 1;
  size_t const m = n + p + 3;
  DoubleVector u; u.reserve(n + 1);
  u.push_back(0.0);
  double total = 0.0;
  for(size_t i = 1; i <= n; ++i) {
    if (centripetal)
      u.push_back(std::sqrt((points[i] - points[i-1]).norm()));
    else
      u.push_back((points[i] - points[i-1]).norm());
    total += u.back();
  }
  for(size_t i = 1; i < n; ++i)
    u[i] = u[i-1] + u[i] / total;
  u.back() = 1.0;

  BSplineCurve result;
  result.p = p; result.n = n + 2;
  result.knots.reserve(m + 1);
  result.cp.reserve(n + 1);

  std::fill_n(std::back_inserter(result.knots), p + 1, 0.0);
  for(size_t j = 0; j <= n - p + 1; ++j) {
    double t = 0.0;
    for(size_t i = j; i <= j + p - 1; ++i)
      t += u[i];
    result.knots.push_back(t / p);
  }
  std::fill_n(std::back_inserter(result.knots), p + 1, 1.0);

  Vector left_tangent = parabolaRule(u[0], points[0], u[1], points[1], u[2], points[2]);
  Vector right_tangent = parabolaRule(u[n], points[n], u[n-1], points[n-1], u[n-2], points[n-2]);
  Point left_control_point =
    points.front() + left_tangent * (result.knots[p+1] - result.knots[p]) / p;
  Point right_control_point =
    points.back() + right_tangent * (result.knots[m-p-1] - result.knots[m-p]) / p;

  Eigen::MatrixXd A(n + 3, n + 3); A.setZero();
  for(size_t i = 0; i <= n; ++i) {
    size_t span = result.findSpan(u[i]);
    DoubleVector coeff; result.basisFunctions(span, u[i], coeff);
    for(size_t j = 0; j <= p; ++j)
      A(i, span - p + j) = coeff[j];
  }
  A(n + 1, 1) = 1.0;
  A(n + 2, n + 1) = 1.0;
  Eigen::MatrixXd b(n + 3, 3);
  for(size_t i = 0; i <= n; ++i) {
    b(i, 0) = points[i].x;
    b(i, 1) = points[i].y;
    b(i, 2) = points[i].z;
  }
  b(n + 1, 0) = left_control_point.x;
  b(n + 1, 1) = left_control_point.y;
  b(n + 1, 2) = left_control_point.z;
  b(n + 2, 0) = right_control_point.x;
  b(n + 2, 1) = right_control_point.y;
  b(n + 2, 2) = right_control_point.z;
  Eigen::MatrixXd x = A.fullPivLu().solve(b);
  for(size_t i = 0; i <= n + 2; ++i)
    result.cp.push_back(Point(x(i, 0), x(i, 1), x(i, 2)));

  return result;
}

BSplineCurve BSplineCurve::interpolateAsInSketches(PointVector const &points, size_t p,
                                                   bool centripetal)
{
  size_t const n = points.size() - 1; // # of data points - 1
  size_t const m = n + p - 1;         // # of control points - 1
  DoubleVector u; u.reserve(n + 1);
  u.push_back(0.0);
  double total = 0.0;
  for(size_t i = 1; i <= n; ++i) {
    if (centripetal)
      u.push_back(std::sqrt((points[i] - points[i-1]).norm()));
    else
      u.push_back((points[i] - points[i-1]).norm());
    total += u.back();
  }
  for(size_t i = 1; i < n; ++i)
    u[i] = u[i-1] + u[i] / total;
  u.back() = 1.0;

  BSplineCurve result;
  result.p = p; result.n = m;
  result.knots.reserve(m + p + 2);
  result.cp.reserve(m + 1);

  std::fill_n(std::back_inserter(result.knots), p + 1, 0.0);
  for(size_t j = 1; j < n; ++j) {
    result.knots.push_back(u[j]);
  }
  std::fill_n(std::back_inserter(result.knots), p + 1, 1.0);

  Vector left_tangent = parabolaRule(u[0], points[0], u[1], points[1], u[2], points[2]);
  Vector right_tangent = parabolaRule(u[n], points[n], u[n-1], points[n-1], u[n-2], points[n-2]);
  Point left_control_point =
    points.front() + left_tangent * (result.knots[p+1] - result.knots[p]) / p;
  Point right_control_point =
    points.back() + right_tangent * (result.knots[m] - result.knots[m+1]) / p;

  Eigen::MatrixXd A(n + 3, m + 1); A.setZero();
  for(size_t i = 0; i <= n; ++i) {
    size_t span = result.findSpan(u[i]);
    DoubleVector coeff; result.basisFunctions(span, u[i], coeff);
    for(size_t j = 0; j <= p; ++j)
      A(i, span - p + j) = coeff[j];
  }
  A(n + 1, 1) = 1.0;
  A(n + 2, m - 1) = 1.0;

  Eigen::MatrixXd b(n + 3, 3);
  for(size_t i = 0; i <= n; ++i) {
    b(i, 0) = points[i].x;
    b(i, 1) = points[i].y;
    b(i, 2) = points[i].z;
  }
  b(n + 1, 0) = left_control_point.x;
  b(n + 1, 1) = left_control_point.y;
  b(n + 1, 2) = left_control_point.z;
  b(n + 2, 0) = right_control_point.x;
  b(n + 2, 1) = right_control_point.y;
  b(n + 2, 2) = right_control_point.z;

  Eigen::MatrixXd x;
  if (m + 1 == n + 3)
    x = A.fullPivLu().solve(b);
  else if (m + 1 > n + 3) {            // underdefined equation system - use least norm solution
    Eigen::MatrixXd At = A.transpose();
    x = At * (A * At).fullPivLu().solve(b);
  } else {                             // overdefined equation system - use LSQ solution
    Eigen::MatrixXd At = A.transpose();
    x = (At * A).fullPivLu().solve(At * b);
  }

  for(size_t i = 0; i <= m; ++i)
    result.cp.push_back(Point(x(i, 0), x(i, 1), x(i, 2)));

  return result;
}

// Simple version - no constraints, weights or derivatives as input.
// Approximates `points' with a B-spline of degree `p' with `n'+1 control points.
// This is essentially the same algorithm as `interpolate', just the knotvector generation differs.
BSplineCurve BSplineCurve::approximate(PointVector const &points, size_t p, size_t n,
                                       bool centripetal)
{
  size_t m = points.size() - 1;
  assert(n <= m);

  BSplineCurve result;
  result.p = p; result.n = n;
  result.knots.reserve(n + p + 2);
  result.cp.reserve(n + 1);

  // Generate parameters for the data points by the centripetal method
  DoubleVector u; u.reserve(m + 1);
  u.push_back(0.0);
  double total = 0.0;
  for(size_t i = 1; i <= m; ++i) {
    if (centripetal)
      u.push_back(std::sqrt((points[i] - points[i-1]).norm()));
    else
      u.push_back((points[i] - points[i-1]).norm());
    total += u.back();
  }
  for(size_t i = 1; i < m; ++i)
    u[i] = u[i-1] + u[i] / total;
  u.back() = 1.0;

  // Generate knotvector
  std::fill_n(std::back_inserter(result.knots), p + 1, 0.0);
  for(size_t j = 1; j <= n - p; ++j) {
    double d = (double)(m + 1) / (n - p + 1);
    size_t i = d * j;
    double alpha = d * j - i;
    double knot = (1.0 - alpha) * u[i-1] + alpha * u[i];
    result.knots.push_back(knot);
  }
  std::fill_n(std::back_inserter(result.knots), p + 1, 1.0);
  
  // Set up the equation
  Eigen::MatrixXd N(m + 1, n + 1); N.setZero();
  Eigen::MatrixXd S(m + 1, 3);
  for (size_t i = 0; i <= m; ++i) {
    size_t span = result.findSpan(u[i]);
    DoubleVector coeff; result.basisFunctions(span, u[i], coeff);
    for(size_t j = 0; j <= p; ++j)
      N(i, span - p + j) = coeff[j];
    S(i, 0) = points[i].x;
    S(i, 1) = points[i].y;
    S(i, 2) = points[i].z;
  }

  // Fill the control points
  Eigen::MatrixXd x = N.fullPivLu().solve(S);
  for(size_t i = 0; i <= n; ++i)
    result.cp.push_back(Point(x(i, 0), x(i, 1), x(i, 2)));

  return result;
}

static BSplineCurve bezierToBSpline(const BezierCurve &curve) {
  BSplineCurve result;
  result.p = curve.n;
  result.n = curve.n;
  std::fill_n(std::back_inserter(result.knots), result.p + 1, 0.0);
  std::fill_n(std::back_inserter(result.knots), result.p + 1, 1.0);
  result.cp = curve.cp;
  return result;
}

static BezierCurve subdivide(const BezierCurve &curve, size_t depth, double alpha) {
  BezierCurve result;
  size_t d = (curve.n + 2 * depth) / (2 * depth + 1); // original degree
  result.n = curve.n + (d - 1) * 2;
  if (!alpha)
    alpha = 2.0 / 3.0;
  result.cp.push_back(curve.cp.front());
  for (size_t j = 1; j < d; ++j) {
    size_t i = j + (2 * j - 1) * depth;
    for (size_t k = i - depth; k < i; ++k)
      result.cp.push_back(curve.cp[k]);
    result.cp.push_back(curve.cp[i-1] * (1 - alpha) + curve.cp[i] * alpha);
    result.cp.push_back(curve.cp[i]);
    result.cp.push_back(curve.cp[i+1] * (1 - alpha) + curve.cp[i] * alpha);
    for (size_t k = i + 1; k <= i + depth; ++k)
      result.cp.push_back(curve.cp[k]);
  }
  result.cp.push_back(curve.cp.back());
  return result;
}

BSplineCurve BSplineCurve::proximity(PointVector const &points, size_t depth, double alpha) {
  BezierCurve result;
  result.n = points.size() - 1;
  result.cp = points;
  if (result.n > 2)
    for (size_t i = 0; i < depth; ++i)
      result = subdivide(result, i, alpha);
  return bezierToBSpline(result);
}

BSplineCurve BSplineCurve::proximityMultiplicity(PointVector const &points, size_t depth) {
  size_t n = points.size();
  BezierCurve result;
  result.n = (n - 2) * (depth + 1) + 1;
  result.cp.push_back(points.front());
  for (size_t i = 1; i < n - 1; ++i)
    std::fill_n(std::back_inserter(result.cp), depth + 1, points[i]);
  result.cp.push_back(points.back());
  return bezierToBSpline(result);
}

[[maybe_unused]]
static Point controlFramePoint(const BSplineCurve &curve, double u) {
  if (u >= curve.knots[curve.knots.size() - curve.p - 1])
    return curve.cp.back();
  double last_greville = 0, greville = 0;
  size_t i = 0;
  for (; i <= curve.n; ++i) {
    last_greville = greville;
    greville = 0;
    for (size_t j = 1; j <= curve.p; ++j)
      greville += curve.knots[i+j];
    greville /= curve.p;
    if (greville > u)
      break;
  }
  double alpha = (u - last_greville) / (greville - last_greville);
  return curve.cp[i-1] * (1 - alpha) + curve.cp[i] * alpha;
}

static PointVector sampleCurve(const BSplineCurve &curve, size_t resolution) {
  PointVector result;
  for (size_t i = 0; i < resolution; ++i) {
    double const t = (double)i / ((double)resolution - 1.0);
    if (curve.knots[curve.p] <= t && t <= curve.knots[curve.knots.size() - curve.p - 1])
      result.push_back(curve.evaluate(t));
  }
  return result;
}

// Assumes uniform parameters
static BezierCurve fitBezier(const BezierCurve &curve, const PointVector &points, double alpha) {
  size_t resolution = points.size();
  BezierCurve result = curve;
  size_t degree = curve.n;

  Eigen::MatrixXd N = Eigen::MatrixXd::Zero(resolution + degree - 1, degree + 1);
  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(resolution + degree - 1, 3);
  for (size_t i = 0; i < resolution; ++i) {
    double const u = (double)i / ((double)resolution - 1.0);
    DoubleVector coeff; result.bernsteinAll(degree, u, coeff);
    for(size_t j = 0; j <= degree; ++j)
      N(i, j) = coeff[j];
    S(i, 0) = points[i].x;
    S(i, 1) = points[i].y;
    S(i, 2) = points[i].z;
  }
  double smoothing = alpha / 10.0;
  for (size_t i = 1; i < degree; ++i) {
    N(resolution + i - 1, i - 1) = smoothing;
    N(resolution + i - 1, i) = -2 * smoothing;
    N(resolution + i - 1, i + 1) = smoothing;
  }

  // Fill the control points
  Eigen::MatrixXd x = N.fullPivLu().solve(S);
  for (size_t i = 1; i < degree; ++i)
    result.cp[i] = { x(i, 0), x(i, 1), x(i, 2) };

#ifdef DEBUG
  double axis = 2.0 * sqrt(2);
  std::cerr << "Error: " << (N * x - S).norm() * 100 / axis << '%' << std::endl;
  double dmax = 0;
  for (size_t i = 0; i < result.n; ++i) {
    double d = (result.cp[i+1] - result.cp[i]).norm();
    if (d > dmax)
      dmax = d;
  }
  std::cerr << "Longest segment: " << dmax * 100 / axis << '%' << std::endl;
#endif

  return result;
}

BSplineCurve BSplineCurve::uniformCubic(const PointVector &points) {
  BSplineCurve result;
  result.p = 3;
  result.n = points.size() - 1;
  std::fill_n(std::back_inserter(result.knots), result.p + 1, 0.0);
  size_t segments = result.n - result.p + 1;
  for (size_t i = 1; i < segments; ++i)
    result.knots.push_back((double)i / segments);
  std::fill_n(std::back_inserter(result.knots), result.p + 1, 1.0);
  result.cp = points;
  return result;
}

BSplineCurve BSplineCurve::proximityFit(PointVector const &points, size_t depth, double alpha) {
  const size_t resolution = 100;
  auto base = uniformCubic(points);
  auto data = sampleCurve(base, resolution);
  BezierCurve result;
  result.n = points.size() - 1;
  result.cp = points;
  for (size_t i = 0; i < depth; ++i)
    result = result.elevate();
  result = fitBezier(result, data, alpha);
  return bezierToBSpline(result);
}

[[maybe_unused]]
static void fillPoints(BezierCurve &curve, size_t degree, size_t k) {
  for (size_t i = 1; i <= degree; ++i) {
    for (size_t j = 1; j < k; ++j) {
      auto p = curve.cp[(i-1)*k];
      auto q = curve.cp[i*k];
      double alpha = (double)j / k;
      if (i == 1) {
        if (j == k - 1)
          continue;
        p = curve.cp[1];
        alpha = (double)j / (k - 1);
      } else if (i == degree) {
        if (j == 1)
          continue;
        q = curve.cp[curve.n-1];
        alpha = (double)(j - 1) / (k - 1);
      }
      curve.cp[i*k-j] = q + (p - q) * alpha;
    }
  }
}

static bool rightTurn(const Point &a, const Point &b, const Point &c) {
  auto u = b - a, v = c - a;
  return u.x * v.y < u.y * v.x;
}

static PointVector convexHull(PointVector points) {
  size_t n = points.size();
  std::sort(points.begin(), points.end(),
            [](const Point &p, const Point &q) {
              return p.x < q.x || (p.x == q.x && p.y < q.y);
            });
  PointVector hull = { points[0], points[1] };
  for (size_t i = 2; i < n; ++i) {
    while (hull.size() >= 2 && !rightTurn(hull.rbegin()[1], hull.back(), points[i]))
      hull.pop_back();
    hull.push_back(points[i]);
  }
  PointVector lower = { points.back(), points.rbegin()[1] };
  for (size_t i = 3; i <= n; ++i) {
    while (lower.size() >= 2 && !rightTurn(lower.rbegin()[1], lower.back(), points[n-i]))
      lower.pop_back();
    lower.push_back(points[n-i]);
  }
  lower.pop_back();
  hull.insert(hull.end(), lower.begin() + 1, lower.end());
  return hull;
}

std::optional<Point> segmentIntersection(const Point &a, const Point &b,
                                         const Point &c, const Point &d) {
  // When the segments have a shared part,
  // this function returns the closest common point to `b`.
  constexpr double eps2 = 1e-10;
  auto A = b - a, C = d - c;
  auto D = A.x * C.y - C.x * A.y;
  if (std::pow(D, 2) < eps2 * (A * A) * (C * C)) {
    // Parallel lines
    auto E = c - a;
    if (std::pow(E.x * A.y - E.y * A.x, 2) < eps2 * (E * E) * (A * A)) {
      auto s0 = A * E / (A * A), s1 = s0 + A * C / (A * A);
      auto smax = std::max(s0, s1);
      if (smax >= 1) {
        if (std::min(s0, s1) <= 1)
          return b;
      } else if (smax >= 0)
        return a + A * smax;
    }
    return std::nullopt;
  }
  auto s = (a.x * (d.y - c.y) + c.x * (a.y - d.y) + d.x * (c.y - a.y)) / D;
  auto t = (a.x * (c.y - b.y) + b.x * (a.y - c.y) + c.x * (b.y - a.y)) / D;
  if (s < 0 || s > 1 || t < 0 || t > 1)
    return std::nullopt;
  return a + A * s;             // = c + C * t
}

static bool insideHull(const PointVector &hull, const Point &p) {
  size_t n = hull.size();
  for (size_t i = 0; i < n; ++i)
    if (!rightTurn(hull[i], hull[(i+1)%n], p))
      return false;
  return true;
}

static Point displacementInCHull(const PointVector &hull, const Point &p, const Vector &d) {
  auto q = p + d;
  return q; // comment to coerce the displacement curve to stay in the convex hull
  if (insideHull(hull, q))
    return q;
  size_t n = hull.size();
  for (size_t i = 0; i < n; ++i) {
    size_t j = (i + 1) % n;
    auto x = segmentIntersection(p, q, hull[i], hull[j]);
    if (x)
      return x.value();
  }
  return p;
}

BSplineCurve BSplineCurve::proximityDisplacement(PointVector const &points, size_t depth,
                                                 double alpha, size_t iterations) {
  // Compute convex hull
  auto hull = convexHull(points);
  // Original curve
  BezierCurve curve;
  curve.n = points.size() - 1;
  curve.cp = points;
  // auto base = uniformCubic(points);
  PBezierCurve base(points, std::max(0.0, std::min(1.0, 1.0 - alpha)));

  BezierCurve result = curve;
  for (size_t i = 0; i < depth; ++i) {
    result = result.elevate();

    for (size_t iter = 0; iter < iterations; ++iter) {
      // Create displacements
      VectorVector displacements;
      for (size_t i = 0; i <= result.n; ++i) {
        double u = (double)i / result.n;
        // displacements.push_back((base.evaluate(u) - result.evaluate(u)) * alpha);
        displacements.push_back(base.evaluate(u) - result.evaluate(u));
      }

      // Add displacements

      // Special case: first control points
      auto d1 = (result.cp[1] - result.cp[0]).unit();
      displacements[1] = d1 * (displacements[1] * d1);
      
      // Special case: last control points
      d1 = (result.cp.rbegin()[1] - result.cp.back()).unit();
      displacements.rbegin()[1] = d1 * (displacements.rbegin()[1] * d1);

      // General case
      for (size_t i = 1; i < result.n; ++i)
        result.cp[i] = displacementInCHull(hull, result.cp[i], displacements[i]);
    }
  }

#ifdef DEBUG
  size_t resolution = 100;
  double dmax = 0;
  for (size_t i = 0; i <= resolution; ++i) {
    double u = (double)i / resolution;
    double d = (base.evaluate(u) - result.evaluate(u)).norm();
    if (d > dmax)
      dmax = d;
  }
  double axis = 2.0 * sqrt(2);
  std::cerr << "Error: " << dmax * 100 / axis << '%' << std::endl;
  dmax = 0;
  for (size_t i = 0; i < result.n; ++i) {
    double d = (result.cp[i+1] - result.cp[i]).norm();
    if (d > dmax)
      dmax = d;
  }
  std::cerr << "Longest segment: " << dmax * 100 / axis << '%' << std::endl;
#endif

  return bezierToBSpline(result);
}

BSplineCurve BSplineCurve::proximitySlider(PointVector const &points, size_t depth) {
  BezierCurve result;
  result.n = points.size() + depth - 1;
  result.cp.resize(points.size() + depth);
  result.cp.front() = points.front();
  result.cp.back() = points.back();
  PBezierCurve base(points, 0.0);

  size_t n = points.size() + depth - 2;

  auto slideClubs = [&](const DoubleVector &x) {
    auto sorted = x;
    std::sort(sorted.begin(), sorted.end());
    for (size_t k = 1; k <= n; ++k)
      result.cp[k] = base.evaluate(sorted[k-1]);
  };

  auto f = [&](const DoubleVector &x) {
    slideClubs(x);
    double err = 0;
    for (size_t k = 1; k <= n; ++k) {
      double u = (double)k / (n + 1);
      err += (base.evaluate(u) - result.evaluate(u)).normSqr();
    }
    return err;
  };

  DoubleVector x;
  for (size_t k = 1; k <= n; ++k)
    x.push_back((double)k / (n + 1));

  double step = 0.1 / (n + 1);
  NelderMead::optimize(f, x, 100, 0.0, step);
  slideClubs(x);

  return bezierToBSpline(result);
}

BSplineCurve BSplineCurve::proximityRational(PointVector const &points, size_t depth) {
  BezierCurve result;
  result.n = points.size() - 1;
  result.cp = points;
  for (size_t i = 1; i < result.n; ++i) {
    result.cp[i] *= depth + 1;
    result.cp[i].z = depth + 1;
  }
  result.cp.front().z = 1;
  result.cp.back().z = 1;
  return bezierToBSpline(result);
}

BSplineCurve BSplineCurve::proximityReduced(PointVector const &points, size_t depth, double alpha) {
  size_t const resolution = 100;
  PBezierCurve base(points, std::max(0.0, std::min(1.0, 1.0 - alpha)));
  BezierCurve result;
  result.n = resolution;
  for (size_t i = 0; i <= resolution; ++i) {
    double u = (double)i / resolution;
    result.cp.push_back(base.evaluate(u));
  }
  // while (result.n > points.size() + depth)
  //   result = result.reduce();
  result = result.reduce(points.size() + depth);
  return bezierToBSpline(result);
}

PBezierCurve::PBezierCurve(const PointVector &cp, double gamma)
  : n(cp.size() - 1), cp(cp), gamma(gamma)
{
}

Point PBezierCurve::evaluate(double u) const {
  Point result(0, 0, 0);
  DoubleVector coeff;
  BezierCurve::bernsteinAll(n, u, coeff);

  // Compute r
  DoubleVector r;
  for (size_t i = 0; i <= n; ++i) {
    double ui = (double)i / n;
    double sum = 0;
    for (size_t j = 0; j < i; ++j)
      sum += (i - j) * coeff[j];
    double fi = std::pow(u - ui + 2.0 / n * sum, 2) - std::pow(u - ui, 2);
    r.push_back(std::sqrt(std::pow(ui - u, 2) + gamma * fi));
  }

  // Compute M
  DoubleVector M;
  M.push_back(0.5 + (r[1] - r[0]) * n / 2);
  for (size_t i = 1; i < n; ++i)
    M.push_back((r[i+1] - 2 * r[i] + r[i-1]) * n / 2);
  M.push_back(0.5 - (r[n] - r[n-1]) * n / 2);

  for (size_t i = 0; i <= n; ++i)
    result += cp[i] * M[i];
  return result;
}
