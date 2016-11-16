#ifndef DSL_UTILS_H
#define DSL_UTILS_H

#include <sys/time.h>
#include <Eigen/Dense>
#include <vector>
#include "utilsimg.h"

namespace dsl {

void save_map(const char* map, int width, int height, const char* filename);

char* load_map(int &width, int &height, const char* filename);

void timer_start(struct timeval* time);

long timer_us(struct timeval* time);

double fangle(double a);

void se2_q2g(Eigen::Matrix3d& m, const Eigen::Vector3d& q);

void se2_g2q(Eigen::Vector3d& q, const Eigen::Matrix3d& m);

void se2_inv(Eigen::Matrix3d& mi, const Eigen::Matrix3d& m);

void se2_exp(Eigen::Matrix3d& m, const Eigen::Vector3d& v, double tol = 1e-16);

void se2_log(Eigen::Vector3d& v, const Eigen::Matrix3d& m, double tol = 1e-16);

void replaceExt(std::string& s, const std::string& newExt);

/**
 * Removes the specified dimension. (1 3 5) = removedDim( (1,2,3,5), 1 );
 * @param in input vector
 * @param dim dimension to remove
 * @return output vector
 */
template<typename T,int n>
Eigen::Matrix<T,n-1,1> removeDim(const Eigen::Matrix<T,n,1>& in, int dim){
  Eigen::Matrix<T,n-1,1> out;
  out = in.template head<n-1>();
  out.tail(n-1-dim) = in.tail(n-1-dim);
  return out;
}

/**
 *Inserts an extra dim and sets the value at that dim. (1 2 3 4 5) = insertDim( (1,2,3,5), 3, 4);
 * @param in input vector
 * @param dim dimension where to insert
 * @param val value to insert
 * @return output vector
 */
template<typename T,int n>
Eigen::Matrix<T,n+1,1> insertDim(const Eigen::Matrix<T,n,1>& in, int dim, T val){
  Eigen::Matrix<T,n+1,1> out;
  out.template head<n>() = in;
  out(dim) = val;
  out.tail(n-dim) = in.tail(n-dim);
  return out;
}

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
}

#endif
