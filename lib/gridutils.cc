#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "gridutils.h"

namespace dsl {

using namespace Eigen;

Map<bool, 2> load(const char* filename, const Vector2d &cs) {
  std::string header;
  int max_col = 0;
  
  std::fstream fs (filename, std::fstream::in);
  assert(fs.is_open());
  fs >> header;
  assert (header == std::string("P6"));

  int width, height;
  fs >> width >> height >> max_col;
  assert(width > 0);
  assert(height > 0);

  dsl::Map<bool, 2> map(Vector2d(0,0), Vector2d(cs[0]*width, cs[1]*height), cs);

  int size = width*height;
  int raster_size = (max_col > 255 ? size*6 : size*3);
  
  char *data = (char*)malloc(raster_size);    
  fs.read(data, raster_size);

  int step = max_col > 255 ? 6 : 3;  
  for (int i = 0; i < size; i++) 
    map.cells[i] = (data[step * i] ? 1 : 0);
  free(data);
  fs.close();
  return map;
}


void save(const dsl::Map<bool, 2> &map, const char* filename, const std::vector<Vector2d> *path) {
  const int &width = map.gs[0];
  const int &height = map.gs[1];
  
  char data[width*height*3];
  std::fstream fs(filename, std::fstream::out);
  assert(fs.is_open());
  fs << "P6" << std::endl << width << " " << height << std::endl << "255" << std::endl;

  int ind = 0;
  for (int i = 0; i < width * height; i++, ind += 3) {
    data[ind] = data[ind + 1] = data[ind + 2] = (char)(map.cells[i] * 100);
  }
  assert(ind == 3*width*height);

  if (path) {
    for (size_t i=0 ; i< path->size(); i++) {
      int id = map.Id(path->at(i));
      int id3 = 3*id;

      if(i==0){
        data[id3] = 255; data[id3+1] = 0; data[id3 + 2] = 0; // Start in green
      }else if(i==path->size()-1){
        data[id3] = 255; data[id3+1] = 0; data[id3 + 2] = 0; // Stop in red
      }else{
        data[id3] = 255; data[id3+1] = 0; data[id3 + 2] = 0; // All other path in blue
      }

    }
  }

  fs.write(data, ind);
  fs.close();
}

void saveMapWithPath(const dsl::Map<bool, 2>& omap, std::string filename,
                     const std::vector<Vector3d>& path, const CarGeom& geom, int scale){


  shared_ptr< dsl::Map<bool,2> > psmap = omap.ScaleUp(scale);
  dsl::Map<bool,2>& smap = *psmap; //scaled up for clarity

  //The image to be saved
  char data[smap.nc*3];
  int ind = 0;
  for (int i = 0; i < smap.nc; i++, ind += 3)
    data[ind] = data[ind + 1] = data[ind + 2] = (char)(smap.cells[i] * 100); //copy the map in grey


  Vector2d xy0(0,0); smap.CellCenter(xy0,Vector2i::Zero()); //find cell center of the cell with index (0,0,0)
  Vector3d axy0(0,xy0(0),xy0(1));

  for(size_t i=0; i<path.size(); i++){
    Vector3d axy = path[i] - axy0;
    Vector2d xy = axy.tail<2>();

    //draw the 4 lines
    Matrix<double,2,4> vs;
    getCarCorners(vs, geom, axy(0));  //get the tilted car corners at origin
    vs.colwise() += xy;               //move the car corners at origin to where it actually is
    vs.row(0) = vs.row(0)/smap.cs[0]; //convert the corners to grid coordinate
    vs.row(1) = vs.row(1)/smap.cs[1]; //convert the corners to grid coordinate
    int seq2[]={3,0,1,2};
    double lwm = 0.05; //meters
    int lw = ceil(lwm/smap.cs[0]);//line width in pixels
    smap.cells = vector<bool>(smap.nc,false); //reset
    for(int i=0;i<4;i++)
      addLine<bool>(smap.cells, smap.gs[0], smap.gs[1], vs.col(i), vs.col(seq2[i]), true,lw);

    if(i==0){
      int ind = 0;
      for (int i = 0; i < smap.nc; i++, ind += 3)
        if(smap.cells[i])
          data[ind + 1] = 255; //copy the path in green
    }else if(i==path.size()-1){
      int ind = 0;
      for (int i = 0; i < smap.nc; i++, ind += 3)
        if(smap.cells[i])
          data[ind] = 255; //copy the path in red
    }else{
      int ind = 0;
      for (int i = 0; i < smap.nc; i++, ind += 3)
        if(smap.cells[i])
          data[ind + 2] = 255; //copy the path in blue
    }
  }

  std::fstream fs(filename, std::fstream::out);
  assert(fs.is_open());
  fs << "P6" << std::endl << smap.gs[0] << " " << smap.gs[1] << std::endl << "255" << std::endl;

  assert(ind == smap.nc*3);

  fs.write(data, ind);
  fs.close();
}

void getCarCorners(Matrix2x4d& verts2d_rotd_pix, const CarGeom& geom, double theta){
  double l = geom.l;
  double b = geom.b;
  double ox = geom.ox;
  double oy = geom.oy;

  Vector2d org_m(-geom.ox, -geom.oy);
  Matrix2x4d verts_wrt_center;
  verts_wrt_center.col(0) << -l / 2, -b / 2; // rb(right back)
  verts_wrt_center.col(1) << l / 2, -b / 2;  // rf(right front)
  verts_wrt_center.col(2) << l / 2, b / 2;   // lf(left front)
  verts_wrt_center.col(3) << -l / 2, b / 2;  // lb(left back)

  // rotate vertices about origin and convert vertices and origin to pixel
  // coordinates
  Transform2d tfm2d = Rotation2Dd(theta) * Translation2d(org_m);
  verts2d_rotd_pix = tfm2d * verts_wrt_center;
}

}
