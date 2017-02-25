#ifndef DSL_LIB_UTILS_H_
#define DSL_LIB_UTILS_H_

#include <sys/time.h>
#include <Eigen/Dense>
#include <vector>
#include "utilsimg.h"
#include <cmath>

namespace dsl {

void save_map(const char* map, int width, int height, const char* filename);

char* load_map(int &width, int &height, const char* filename);

void timer_start(struct timeval* time);

long timer_us(struct timeval* time);

double fangle(double a);

/**
 * Convert SE2 Matrix to its coordinate representation
 * @param m
 * @param q
 */
void se2_q2g(Eigen::Matrix3d& m, const Eigen::Vector3d& q);

/**
 * Obtain SE2 Matrix from its coordinate representation
 * @param q
 * @param m
 */
void se2_g2q(Eigen::Vector3d& q, const Eigen::Matrix3d& m);

/**
 * Invert a SE2 Matrix. Not same as matrix inverse.
 * @param mi
 * @param m
 */
void se2_inv(Eigen::Matrix3d& mi, const Eigen::Matrix3d& m);

/**
 * Exponentiate the se2 element(twist) to obtain SE2 Matrix
 * @param m
 * @param v
 * @param tol
 */
void se2_exp(Eigen::Matrix3d& m, const Eigen::Vector3d& v, double tol = 1e-16);

/**
 * Take the log of the SE2 Matrix to obtain se2 element(twist)
 * @param v
 * @param m
 * @param tol
 */
void se2_log(Eigen::Vector3d& v, const Eigen::Matrix3d& m, double tol = 1e-16);

/**
 * se2 twist(only w and vx) that takes you exactly only to (xf,yf) and not angle
 * @param xf
 * @param yf
 * @return twist(only w and vx). vy=0
 */
Eigen::Vector2d se2_get_wvx( double xf,double yf);

/**
 * Removes the specified dimension. [10,56,75] = removedDim( [10,14,56,75], 1 );
 * @param in input std::vector
 * @param dim dimension to remove
 * @return output std::vector
 */
template<typename T,int n>
Eigen::Matrix<T,n-1,1> RemoveDimension(const Eigen::Matrix<T,n,1>& in, int dim){
  assert(dim>=0 && dim <n);
  dim = dim<0 ? 0 : dim;
  dim = dim>=n ? n-1 : dim;

  Eigen::Matrix<T,n-1,1> out;
  out = in.template head<n-1>();
  out.tail(n-1-dim) = in.tail(n-1-dim);
  return out;
}

/**
 *Inserts an extra dim and sets the value at that dim.
 *e.g. [100, 10, 14, 56, 75] = InsertDimension( [10,14,56,75], 100, 0);
 *e.g. [ 10,100, 14, 56, 75] = InsertDimension( [10,14,56,75], 100, 1);
 *e.g. [ 10, 14, 56,100, 75] = InsertDimension( [10,14,56,75], 100, 3);
 *e.g. [ 10, 14, 56, 75,100] = InsertDimension( [10,14,56,75], 100, 4);
 * @param in input std::vector
 * @param dim dimension where to insert
 * @param val value to insert
 * @return output std::vector
 */
template<typename T,int n>
Eigen::Matrix<T,n+1,1> InsertDimension(const Eigen::Matrix<T,n,1>& in, int dim, T val){
  assert(dim>=0 && dim <=n);
  dim = dim<0 ? 0 : dim;
  dim = dim>n ? n : dim;

  Eigen::Matrix<T,n+1,1> out;
  out.template head<n>() = in;
  out(dim) = val;
  out.tail(n-dim) = in.tail(n-dim);
  return out;
}


/**
 * Find the barycentre coordinates for a triangle
 */
template<typename T>
class Barycentre3{
public:
  using Vector3T = Eigen::Matrix<T,3,1>;
  Barycentre3(const Vector3T& p0, const Vector3T& p1, const Vector3T& p2 )
    : p0_(p0), v01_(p1 - p0), v02_(p2 - p0) {
    d00_ = v01_.dot(v01_);
    d01_ = v01_.dot(v02_);
    d11_ = v02_.dot(v02_);

    denom_ = d00_ * d11_ - d01_ * d01_;
  }

  Vector3T Find(const Vector3T& p){
    Vector3T v = p - p0_;
    T d20 = v.dot(v01_);
    T d21 = v.dot(v02_);

    Vector3T uvw;
    uvw(1) = (d11_ * d20 - d01_ * d21) / denom_;
    uvw(2) = (d00_ * d21 - d01_ * d20) / denom_;
    uvw(0) = 1.0f - uvw(1) - uvw(2);

    return uvw;
  }

private:
  Vector3T p0_, v01_, v02_;
  T d00_, d01_, d11_, denom_;
};
using Barycentre3d = Barycentre3<double>;
using Barycentre3f = Barycentre3<float>;

/**
 * Function that returns a sign for object of any class as long as the operator
 * - and operator < are defined
 * for that class
 * @param val the object whose sign we need to check
 * @return sign of the object
 */
template < typename T >
int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

/**
 * Change the extension of a std::string.
 * @param filename The filename whose extension we want to change
 * @param new_extension new extension
 * @return filename with new extension
 */
std::string ReplaceExtension(const std::string& filename, const std::string& new_extension);

}

#endif
