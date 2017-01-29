#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "utils.h"

namespace dsl {

using namespace Eigen;
using namespace std;

void save_map(const char* map, int width, int height, const char* filename) {
  int i, ind;
  char data[width*height*3];
  //FILE* file = fopen(filename, "w");
  //assert(file);
  std::fstream fs(filename, std::fstream::out);
  assert(fs.is_open());
  fs << "P6" << std::endl << width << " " << height << std::endl << "255" << std::endl;
  //fprintf(file, "P6\n%d %d\n255\n", width, height);

  ind = 0;
  for (i = 0; i < width * height; i++, ind += 3) {
    data[ind] = data[ind + 1] = data[ind + 2] = (char)(map[i] * 100);
  }
  assert(ind == 3*width*height);
  fs.write(data, ind);
  //assert((int)fwrite(data, sizeof(char), ind, file) == ind);
  //fclose(file);
  fs.close();
}

char* load_map(int &width, int &height, const char* filename) {
  int i, size;

  char *map, *data;
  std::string header;

  int max_col = 0;

  //FILE* file = fopen(filename, "r");
  //assert(file);
  //assert(fscanf(file, "P6\n%d %d 255\n", width, height));
  std::fstream fs (filename, std::fstream::in);
  assert(fs.is_open());
  fs >> header;
  assert (header == std::string("P6"));

  fs >> width >> height >> max_col;
  assert(width > 0);
  assert(height > 0);

  size = width*height;
  int raster_size = (max_col > 255 ? size*6 : size*3);

  data = (char*)malloc(raster_size);
  map = (char*)malloc(size);

  fs.read(data, raster_size);

  int step = max_col > 255 ? 6 : 3;

  for (i = 0; i < size; i++) 
    map[i] = (data[step * i] ? 1 : 0);
  free(data);
  fs.close();
  //fclose(file);
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

  if (std::fabs(w) < tol) {
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

void se2_log(Vector3d &v, const Matrix3d &m, double tol){
  double th = atan2(m(1,0),m(0,0));
  v(0) = th;
  const double& x = m(0,2);
  const double& y = m(1,2);
  if (std::fabs(th) < tol) {
    v(1) = x;
    v(2) = y;
    return;
  }
  double th2 = th/2;
  double a = th2/tan(th/2);
  v(1) = a*x + th2*y;
  v(2) = -th2*x + a*y;
}

Vector2d se2_get_wvx( double xf,double yf){
  double w, t;
  if(std::fabs(yf)<1e-12){
    w=0;
    t=xf;
  }else{
    w = 2*yf/(xf*xf+yf*yf);
    t = atan2(w*xf, 1-w*yf)/w;
  }
  return Vector2d(w*t,t);
}

std::string ReplaceExtension(const std::string& filename, const std::string& new_extension) {
  string filename_new = filename;
  std::string::size_type i = filename.rfind('.', filename.length());
  if (i != std::string::npos)
    filename_new.replace(i+1, new_extension.length(), new_extension);
  return filename_new;
}

}
