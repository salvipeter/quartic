#include "vector.hh"

#include <cstddef>
#include <list>

typedef Vector Point;
typedef std::vector<Point> PointVector;
typedef std::vector<Vector> VectorVector;
typedef std::vector<double> DoubleVector;
typedef std::vector<DoubleVector> DoubleMatrix;
typedef std::vector<PointVector> PointMatrix;

struct Curve
{
  virtual Point evaluate(double u) const = 0;
  virtual ~Curve() { }
};

struct BezierCurve : public Curve
{
  size_t n;                     // degree => n + 1 control points
  PointVector cp;               // control points

  static double bernstein(size_t i, size_t n, double u);
  Point evaluateOneByOne(double u) const;
  static void bernsteinAll(size_t n, double u, DoubleVector &coeff);
  Point evaluate(double u) const override;
  Point evaluateWithCachedCofficients(DoubleVector const &coeff) const;
  Point evaluateByDeCasteljau(double u) const;
  void derivativeControlPoints(size_t d, PointMatrix &dcp) const;
  static void bernsteinAll(size_t n, double u, DoubleMatrix &coeff);
  Point derivativesByControlPoints(double u, size_t d, VectorVector &der) const;
  bool hasInflections() const;
  BezierCurve elevate() const;
  static BezierCurve interpolateUniform(const PointVector &points);
};

struct BSplineCurve : public Curve
{
  size_t p;                     // degree
  size_t n;                     // n + 1 = cp.size()
  DoubleVector knots;           // first and last p+1 values are the same ("clamped")
  PointVector cp;               // knots.size() = cp.size() + p + 1

  size_t findSpan(double u) const;
  void basisFunctions(size_t i, double u, DoubleVector &coeff) const;
  Point evaluate(double u) const override;
  void basisFunctionDerivatives(size_t i, double u, size_t d, DoubleMatrix &der) const;
  Point derivatives(double u, size_t d, VectorVector &der) const;
  void derivativeControlPoints(size_t d, size_t r1, size_t r2, PointMatrix &dcp) const;
  void basisFunctionsAll(size_t i, double u, DoubleMatrix &coeff) const;
  Point derivativesByControlPoints(double u, size_t d, VectorVector &der) const;
  size_t findSpanWithMultiplicity(double u, size_t &multi) const;
  Point evaluateByKnotInsertion(double u) const;
  Point evaluate2DRational(double u) const;
  static size_t binomial(size_t n, size_t k);
  Point derivatives2DRational(double u, size_t d, VectorVector &der) const;
  BSplineCurve insertKnot(double u, size_t k, size_t s, size_t r) const;
  BSplineCurve refineKnots(DoubleVector new_knots) const;
  Point projectPoint(Point const &p, double &u, double &distance,
                     size_t resolution, double distance_tol, double cosine_tol) const;
  std::list<BezierCurve> convertToBezierCurves() const;
  static BSplineCurve interpolate(PointVector const &points, size_t degree, bool centripetal);
  static Vector parabolaRule(double u1, Point const &p1,
                             double u2, Point const &p2,
                             double u3, Point const &p3);
  static BSplineCurve interpolateWithParabolaRule(PointVector const &points, size_t degree,
                                                  bool centripetal);
  static BSplineCurve interpolateAsInSketches(PointVector const &points, size_t degree,
                                              bool centripetal);
  static BSplineCurve approximate(PointVector const &points, size_t degree, size_t n,
                                  bool centripetal);
  static BSplineCurve uniformCubic(const PointVector &points);
  static BSplineCurve proximity(PointVector const &points, size_t depth, double alpha);
  static BSplineCurve proximityFit(PointVector const &points, size_t depth, double alpha);
  static BSplineCurve proximityDisplacement(PointVector const &points, size_t depth, double alpha,
                                            size_t iterations);
};

struct PBezierCurve : public Curve
{
  size_t n;                     // degree => n + 1 control points
  PointVector cp;               // control points
  double gamma;                 // proximity - 1 for Bezier, 0 for control polygon

  PBezierCurve(const PointVector &cp, double gamma);
  Point evaluate(double u) const override;
};
