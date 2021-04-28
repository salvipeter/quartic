#include <cmath>
#include <vector>

struct Vector
{
  double x, y, z;
  Vector(double x = 0.0, double y = 0.0, double z = 0.0) : x(x), y(y), z(z) { }
  Vector &operator=(Vector const &p)      { x = p.x; y = p.y; z = p.z; return *this; }
  double  normSqr()                  const { return x * x + y * y + z * z; }
  double  norm()                     const { return std::sqrt(x * x + y * y + z * z); }
  Vector  unit()                     const { return (*this) / norm(); }
  void    normalize()                      { (*this) = unit(); }
  Vector  operator+(Vector const &p) const { return Vector(x + p.x, y + p.y, z + p.z); }
  Vector &operator+=(Vector const &p)      { x += p.x; y += p.y; z += p.z; return *this; }
  Vector  operator-(Vector const &p) const { return Vector(x - p.x, y - p.y, z - p.z); }
  Vector &operator-=(Vector const &p)      { x -= p.x; y -= p.y; z -= p.z; return *this; }
  Vector  operator*(double scale)    const { return Vector(x * scale, y * scale, z * scale); }
  Vector &operator*=(double scale)         { x *= scale; y *= scale; z *= scale; return *this; }
  Vector  operator/(double scale)    const { return Vector(x / scale, y / scale, z / scale); }
  Vector &operator/=(double scale)         { x /= scale; y /= scale; z /= scale; return *this; }
  double  operator*(Vector const &p) const { return x * p.x + y * p.y + z * p.z; }
  Vector  operator^(Vector const &p) const { return Vector(y * p.z - z * p.y,
                                                           z * p.x - x * p.z,
                                                           x * p.y - y * p.x); }
};
