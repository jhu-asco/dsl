#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "utils.h"

namespace dsl {

using namespace Eigen;

void save_map(const char* map, int width, int height, const char* filename) {
  int i, ind;
  char data[width * height * 3];
  FILE* file = fopen(filename, "w");
  assert(file);
  fprintf(file, "P6\n%d %d 255\n", width, height);
  ind = 0;
  for (i = 0; i < width * height; i++, ind += 3) {
    data[ind] = data[ind + 1] = data[ind + 2] = (char)(map[i] * 100);
  }
  assert(ind == 3 * width * height);
  assert((int)fwrite(data, sizeof(char), ind, file) == ind);
  fclose(file);
}

char* load_map(int* width, int* height, const char* filename) {
  int i, size;
  char* map, *data;
  FILE* file = fopen(filename, "r");
  assert(file);
  assert(fscanf(file, "P6\n%d %d 255\n", width, height));
  size = (*width * *height);
  map = (char*)malloc(size);
  data = (char*)malloc(size * 3);
  assert(fread(data, sizeof(char), size * 3, file));
  for (i = 0; i < size; i++)
    map[i] = (data[3 * i] ? 1 : 0);
  free(data);
  fclose(file);
  return map;
}

void timer_start(struct timeval* time) {
  gettimeofday(time, (struct timezone*)0);
}

long timer_us(struct timeval* time) {
  struct timeval now;
  gettimeofday(&now, (struct timezone*)0);
  return 1000000 * (now.tv_sec - time->tv_sec) + now.tv_usec - time->tv_usec;
}

double fangle(double a) {
  double d = fmod(a, 2 * M_PI);
  if (d > M_PI)
    d -= 2 * M_PI;
  else if (d < -M_PI)
    d += 2 * M_PI;
  // theta should now be in [-pi,pi]
  return d;
}

void se2_q2g(Matrix3d& m, const Vector3d& q) {
  double ct = cos(q[0]);
  double st = sin(q[0]);

  m(0, 0) = ct;
  m(0, 1) = -st;
  m(0, 2) = q[1];
  m(1, 0) = st;
  m(1, 1) = ct;
  m(1, 2) = q[2];
  m(2, 0) = 0;
  m(2, 1) = 0;
  m(2, 2) = 1;
}

void se2_g2q(Vector3d& q, const Matrix3d& m) {
  q[0] = atan2(m(1, 0), m(0, 0));
  q[1] = m(0, 2);
  q[2] = m(1, 2);
}

void se2_inv(Matrix3d& mi, const Matrix3d& m) {
  const double& ct = m(0, 0);
  const double& st = m(1, 0);
  const double& x = m(0, 2);
  const double& y = m(1, 2);

  mi(0, 0) = ct;
  mi(0, 1) = st;
  mi(0, 2) = -x * ct - y * st;
  mi(1, 0) = -st;
  mi(1, 1) = ct;
  mi(1, 2) = x * st - y * ct;
  mi(2, 0) = 0;
  mi(2, 1) = 0;
  mi(2, 2) = 1;
}

void se2_exp(Matrix3d& m, const Vector3d& v, double tol) {
  const double& w = v[0];

  if (fabs(w) < tol) {
    m(0, 0) = 1;
    m(0, 1) = 0;
    m(0, 2) = v[1];
    m(1, 0) = 0;
    m(1, 1) = 1;
    m(1, 2) = v[2];
    m(2, 0) = 0;
    m(2, 1) = 0;
    m(2, 2) = 1;
    return;
  }

  double c = cos(w);
  double s = sin(w);

  double ax = v[2] / w;
  double ay = -v[1] / w;

  m(0, 0) = c;
  m(0, 1) = -s;
  m(0, 2) = (c - 1) * ax - s * ay;
  m(1, 0) = s;
  m(1, 1) = c;
  m(1, 2) = s * ax + (c - 1) * ay;
  m(2, 0) = 0;
  m(2, 1) = 0;
  m(2, 2) = 1;
}
}
