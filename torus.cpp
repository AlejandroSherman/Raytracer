#include "torus.h"
#include "ray.h"
#include <algorithm>
#include <type_traits>
#include <array>
#include <cmath>

// class QuarticEquation
// {
// private:
//   const double d4, d3, d2, d1, d0;
//
// public:
//   QuarticEquation(double d4, double d3, double d2, double d1, double d0) :
//       d4(d4), d3(d3), d2(d2), d1(d1), d0(d0) { }
//
//   void Function(const double *x, double *y) const
//   {
//     const double x2 = x[0]*x[0];
//     const double x3 = x[0]*x2;
//     const double x4 = x[0]*x3;
//     y[0] = d4*x4 + d3*x3 + d2*x2 + d1*x[0] + d0;
//   }
//   void Jacobian(const double *x, double *J) const
//   {
//     const double x2 = x[0]*x[0];
//     const double x3 = x[0]*x2;
//     J[0] = 4*d4*x3 + 3*d3*x2 + 2*d2*x[0] + d1;
//   }
//   int Solve(double *const roots) const;
//
// private:
//   void solve_cubic_equation
//   (double c3, double c2, double c1, double c0,
//    double &x1, double &x2, double &x3, int &nr) const;
// };
//
// int QuarticEquation::Solve(double *const roots) const
// {
//   double r1, r2, r3, r4;
//   int nr12, nr34;
//
//   double a3 = d3/d4;
//   double a2 = d2/d4;
//   double a1 = d1/d4;
//   double a0 = d0/d4;
//
//   double au2 = -a2;
//   double au1 = (a1*a3 - 4.0*a0) ;
//   double au0 = 4.0*a0*a2 - a1*a1 - a0*a3*a3;
//
//   double x1, x2, x3;
//   int nr;
//   solve_cubic_equation(1.0, au2, au1, au0, x1, x2, x3, nr);
//
//   double u1;
//   if (nr==1) {
//     u1 = x1;
//   }
//   else {
//     u1 = (x1>x3) ? x1 : x3;
//   }
//
//   double R2 = 0.25*a3*a3 + u1 - a2;
//   double R = (R2>0.0) ? sqrt(R2) : 0.0;
//
//   double D2, E2;
//   if (R != 0.0) {
//     double foo1 = 0.75*a3*a3 - R2 - 2.0*a2;
//     double foo2 = 0.25*(4.0*a3*a2 - 8.0*a1 - a3*a3*a3) / R;
//     D2 = foo1 + foo2;
//     E2 = foo1 - foo2;
//   }
//   else {
//     double foo1 = 0.75*a3*a3 - 2.0*a2;
//     double foo2 = 2.0 * sqrt(u1*u1 - 4.0*a0);
//     D2 = foo1 + foo2;
//     E2 = foo1 - foo2;
//   }
//
//   if (D2 >= 0.0) {
//     double D = sqrt(D2);
//     r1 = -0.25*a3 + 0.5*R - 0.5*D;
//     r2 = -0.25*a3 + 0.5*R + 0.5*D;
//     nr12 = 2;
//   }
//   else {
//     r1 = r2 = -0.25*a3 + 0.5*R;
//     nr12 = 0;
//   }
//
//   if (E2 >= 0.0) {
//     double E = sqrt(E2);
//     r3 = -0.25*a3 - 0.5*R - 0.5*E;
//     r4 = -0.25*a3 - 0.5*R + 0.5*E;
//     nr34 = 2;
//   }
//   else {
//     r3 = r4 = -0.25*a3 - 0.5*R;
//     nr34 = 0;
//   }
//
//   int i=0;
//   if (nr12 != 0) {
//     roots[i++] = r1;
//     roots[i++] = r2;
//   }
//   if (nr34 != 0) {
//     roots[i++] = r3;
//     roots[i++] = r4;
//   }
//
//   return nr12 + nr34;
// }
//
// void QuarticEquation::solve_cubic_equation
// (double  c3, double  c2, double  c1, double c0,
//  double& x1, double& x2, double& x3, int& nr) const
// {
//   double a2 = c2/c3;
//   double a1 = c1/c3;
//   double a0 = c0/c3;
//
//   double q = a1/3.0 - a2*a2/9.0;
//   double r = (a1*a2 - 3.0*a0)/6.0 - a2*a2*a2 / 27.0;
//   double delta = q*q*q + r*r;
//
//   if (delta>0.0) {
//     double s1 = r + sqrt(delta);
//     s1 = (s1>=0.0) ? pow(s1,1./3.) : -pow(-s1,1./3.);
//
//     double s2 = r - sqrt(delta);
//     s2 = (s2>=0.0) ? pow(s2,1./3.) : -pow(-s2,1./3.);
//
//     x1 = (s1+s2) - a2/3.0;
//     x2 = x3 = -0.5 * (s1+s2) - a2/3.0;
//
//     nr = 1;
//   }
//   else if (delta < 0.0) {
//     double theta = acos(r/sqrt(-q*q*q)) / 3.0;
//     double costh = cos(theta);
//     double sinth = sin(theta);
//     double sq = sqrt(-q);
//
//     x1 = 2.0*sq*costh - a2/3.0;
//     x2 = -sq*costh - a2/3.0 - sqrt(3.) * sq * sinth;
//     x3 = -sq*costh - a2/3.0 + sqrt(3.) * sq * sinth;
//
//     nr = 3;
//   }
//   else {
//     double s = (r>=0.0) ? pow(r,1./3.) : -pow(-r,1./3.);
//     x1 = 2.0*s - a2/3.0;
//     x2 = x3 = -s - a2/3.0;
//     nr = 3;
//   }
// }




Torus::Torus(const vec3& center,const vec3& axis,
    double inner_radius,double outer_radius)
    :center(center),axis(axis),inner_radius(inner_radius),
    outer_radius(outer_radius)
{
    w = gsl_poly_complex_workspace_alloc(5);
}

Torus::~Torus()
{
    gsl_poly_complex_workspace_free(w);
}

// Determine if the ray intersects with the torus
void Torus::Intersection(const Ray& ray, std::vector<Hit>& hits) const
{
  TODO;
  // Vector rayOriginPosition = ray.getOriginPosition();
  // Vector rayDirection = ray.getDirection();
  //
  // Vector centerToRayOrigin = rayOriginPosition - mCenter;
  // const float centerToRayOriginDotDirection = rayDirection.dotProduct(centerToRayOrigin);
  // float	centerToRayOriginDotDirectionSquared = centerToRayOrigin.dotProduct(centerToRayOrigin);
  // float innerRadiusSquared = mInnerRadius * mInnerRadius;
  // float outerRadiusSquared = mOuterRadius * mOuterRadius;
  //
  // float	axisDotCenterToRayOrigin	= mAxis.dotProduct(centerToRayOrigin);
  // float	axisDotRayDirection = mAxis.dotProduct(rayDirection);
  // float	a = 1 - axisDotRayDirection * axisDotRayDirection;
  // float	b = 2 * (centerToRayOrigin.dotProduct(rayDirection) - axisDotCenterToRayOrigin * axisDotRayDirection);
  // float c = centerToRayOriginDotDirectionSquared - axisDotCenterToRayOrigin * axisDotCenterToRayOrigin;
  // float	d = centerToRayOriginDotDirectionSquared + outerRadiusSquared - innerRadiusSquared;
  //
  // // Solve quartic equation with coefficients A, B, C, D and E
  // float A = 1;
  // float B = 4 * centerToRayOriginDotDirection;
  // float C = 2 * d + B * B * 0.25f - 4 * outerRadiusSquared * a;
  // float D = B * d - 4 * outerRadiusSquared * b;
  // float E = d * d - 4 * outerRadiusSquared * c;
  //
  // // Maximum number of roots is 4
  // QuarticEquation equation(A, B, C, D, E);
  // const int maxRootsCount = 4;
  // double roots[maxRootsCount] = {-1.0, -1.0, -1.0, -1.0};
  // int rootsCount = equation.Solve(roots);
  //
  // if (rootsCount == 0) {
  //   return RayIntersection();
  // }
  //
  // // Find closest to zero positive solution
  // float closestRoot = MAX_DISTANCE_TO_INTERSECTON;
  // std::vector<float> intersectionDistances;
  // for (int idx = 0; idx < maxRootsCount; ++idx) {
  //   float root = roots[idx];
  //   if (root > FLOAT_ZERO && root < closestRoot) {
  //     closestRoot = root;
  //     intersectionDistances.push_back(root);
  //   }
  // }
  //
  // if (closestRoot != MAX_DISTANCE_TO_INTERSECTON) {
  //   std::sort(intersectionDistances.begin(), intersectionDistances.end());
  //   TorusPointer pointer = TorusPointer(new Torus(*this));
  //   return RayIntersection(true, pointer, closestRoot, getNormal(ray, closestRoot), intersectionDistances);
  // }
  //
  // return RayIntersection();
}

vec3 Torus::Normal(const vec3& point) const
{
    vec3 normal;
    //TODO; // compute the normal direction
    vec3 centerToPoint = point - center;
    double	centerToPointDotAxis = dot(centerToPoint, axis);
    vec3 direction = centerToPoint - axis * centerToPointDotAxis;
    direction.normalized();
    normal = point - center + direction * outer_radius;
    normal.normalized();
    return normal;
}

void Torus::Update_Bounding_Box()
{
    double x=axis[0]*axis[0],y=axis[1]*axis[1],z=axis[2]*axis[2];
    vec3 u(sqrt(y+z),sqrt(x+z),sqrt(x+y));
    vec3 v(u*outer_radius+inner_radius);
    bounding_box.lo=center-v;
    bounding_box.hi=center+v;
    infinite_box=false;
}
