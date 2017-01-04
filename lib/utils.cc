#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "utils.h"

namespace dsl {

using namespace Eigen;

void ImageRGB::changeBitDepth(BitDepth bitdepth_new){
  switch(bitdepth){
    case BitDepth::BD8:
      if(bitdepth_new == BitDepth::BD8){
        return;
      }else{
        for(int id = 0; id < w*h; id++){
          rdata[id] = rdata[id]<<8;
          gdata[id] = gdata[id]<<8;
          bdata[id] = bdata[id]<<8;
        }
      }
      break;
    case BitDepth::BD16:
      if(bitdepth_new == BitDepth::BD16){
        return;
      }else{
        for(int id = 0; id < w*h; id++){
          rdata[id] = rdata[id]>>8;
          gdata[id] = gdata[id]>>8;
          bdata[id] = bdata[id]>>8;
        }
      }
      break;
  }
  bitdepth = bitdepth_new;
}

void removeComment(ifstream &f){
  char linebuf[1024];
  char ppp;
  while (ppp = f.peek(), ppp == '\n' || ppp == '\r')
    f.get();
  if (ppp == '#')
    f.getline(linebuf, 1023);
}

bool loadPPM(ImageRGB &img, const string &filename){
  if(filename.compare(filename.size() - 4, 4, ".ppm")){
    cout<<"File doesn't have .ppm extension"<<endl;
    return false;
  }

  ifstream f(filename.c_str(), ios::binary);
  if (f.fail()){
    cout << "Could not open file: " << filename << endl;
    return false;
  }

  // get type of file
  removeComment(f);
  string header;
  f >> header;

  // get w
  removeComment(f);
  f >> img.w;

  // get h
  removeComment(f);
  f >> img.h;

  // get bitdepth
  removeComment(f);
  long bitdepth = 0;
  f >> bitdepth;
  img.bitdepth = bitdepth;

  // error checking
  if (header != "P6"){
    cout << "Unsupported magic number" << endl;
    f.close();
    return false;
  }

  if (img.w < 1){
    cout << "Unsupported width: " << img.w << endl;
    f.close();
    return false;
  }

  if (img.h < 1){
    cout << "Unsupported height: " << img.h << endl;
    f.close();
    return false;
  }

  if (bitdepth!=0xff && bitdepth != 0xffff){
    cout << "Unsupported bitdepth: " << bitdepth << endl;
    f.close();
    return false;
  }

  // load image data
  img.rdata.resize(img.w * img.h);
  img.gdata.resize(img.w * img.h);
  img.bdata.resize(img.w * img.h);

  f.get();

  if(bitdepth == 0xff){
    vector<uint8_t> data(img.w * img.h * 3);
    f.read((char*)&data[0], data.size());
    for(int i=0; i< img.w * img.h ; i++){
      int id = i*3;
      img.rdata[i] = data[id];
      img.gdata[i] = data[id+1];
      img.bdata[i] = data[id+2];
    }
  }

  if(bitdepth == 0xffff){
    vector<uint8_t> data(img.w * img.h * 3 *2); //3 channel and 2 byte per channel
    f.read((char*)&data[0], data.size());
    uint16_t temp;
    for(int i=0; i< img.w * img.h ; i++){
      int id = i*6;
      //Most Significant Byte fist according to Netpbm standard
      temp = data[id];
      img.rdata[i] = (temp<<8) | data[id+1];
      temp = data[id+2];
      img.gdata[i] = (temp<<8) | data[id+3];
      temp = data[id+4];
      img.bdata[i] = (temp<<8) | data[id+5];
    }

  }

  // close file
  f.close();

  return true;
}

bool savePPM(ImageRGB& img,  const string& filename) {

  if(filename.compare(filename.size() - 4, 4, ".ppm")){
    cout<<"File doesn't have .ppm extension"<<endl;
    return false;
  }

  std::fstream fs(filename, std::fstream::out);
  if(!fs.is_open()){
    cout<<"couldn't open the image file to write"<<endl;
    return false;
  }

  if(img.bitdepth == 0xff){
    fs << "P6" << std::endl << img.w << " " << img.h << std::endl << "255" << std::endl;
    uint8_t data[img.h*img.w*3];
    for(int i=0; i< img.w * img.h ; i++){
      int id = i*3;
      data[id] = (uint8_t)img.rdata[i];
      data[id + 1] = (uint8_t)img.gdata[i];
      data[id + 2] = (uint8_t)img.bdata[i];
    }
    fs.write((char*)data, img.h*img.w*3);
  }else if(img.bitdepth == 0xffff ){
    fs << "P6" << std::endl << img.w << " " << img.h << std::endl << "65535" << std::endl;
    uint8_t data[img.h*img.w*3*2];
    for(int i=0; i< img.w * img.h ; i++){
      int id = i*6;
      //Most Significant Byte fist according to Netpbm standard
      data[id]   = (uint8_t)(img.rdata[i]>>8);
      data[id+1] = (uint8_t)img.rdata[i];
      data[id+2] = (uint8_t)(img.gdata[i]>>8);
      data[id+3] = (uint8_t)img.gdata[i];
      data[id+4] = (uint8_t)(img.bdata[i]>>8);
      data[id+5] = (uint8_t)img.bdata[i];
    }
    fs.write((char*)data, img.h*img.w*6);
  }else{
    cout<<"The bitdepth of the image is not supported"<<endl;
    return false;
  }

  fs.close();
  return true;
}

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

void se2_log(Vector3d &v, const Matrix3d &m, double tol){
  double th = atan2(m(1,0),m(0,0));
  v(0) = th;
  const double& x = m(0,2);
  const double& y = m(1,2);
  if (fabs(th) < tol) {
    v(1) = x;
    v(2) = y;
    return;
  }
  double th2 = th/2;
  double a = th2/tan(th/2);
  v(1) = a*x + th2*y;
  v(2) = -th2*x + a*y;
}

Vector2d getWVx( double xf,double yf){
  double w, t;
  if(abs(yf)<1e-12){
    w=0;
    t=xf;
  }else{
    w = 2*yf/(xf*xf+yf*yf);
    t = atan2(w*xf, 1-w*yf)/w;
  }
  return Vector2d(w*t,t);
}

void replaceExt(std::string& s, const std::string& newExt) {
  std::string::size_type i = s.rfind('.', s.length());
  if (i != std::string::npos) {
    s.replace(i+1, newExt.length(), newExt);
  }
}

}
