#ifndef DSL_UTILS_H
#define DSL_UTILS_H

#include <sys/time.h>
#include <Eigen/Dense>

namespace dsl {
  
  void save_map(const char* map, int width, int height, const char* filename);
  
  char* load_map(int* width, int* height, const char* filename);
  
  void timer_start(struct timeval *time);
  
  long timer_us(struct timeval *time) ;

  double fangle(double a);

  void se2_q2g(Eigen::Matrix3d &m, const Eigen::Vector3d &q);

  void se2_g2q(Eigen::Vector3d &q, const Eigen::Matrix3d &m);

  void se2_inv(Eigen::Matrix3d &mi, const Eigen::Matrix3d &m);

  void se2_exp(Eigen::Matrix3d &m, const Eigen::Vector3d &v, double tol = 1e-16);
  
  
}

#endif
