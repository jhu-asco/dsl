#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "gridutils.h"
#include <numeric>

namespace dsl {

using namespace Eigen;

Map<bool, 2>::Ptr load(const string& filename, const Vector2d &cs) {
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

  dsl::Map<bool, 2>::Ptr pmap(new Map<bool, 2>(Vector2d(0,0), Vector2d(cs[0]*width, cs[1]*height), cs));

  int size = width*height;
  int raster_size = (max_col > 255 ? size*6 : size*3);

  char *data = (char*)malloc(raster_size);    
  fs.read(data, raster_size);

  int step = max_col > 255 ? 6 : 3;  
  for (int i = 0; i < size; i++) 
    pmap->cells[i] = (data[step * i] ? 1 : 0);
  free(data);
  fs.close();
  return pmap;
}


void save(const dsl::Map<bool, 2> &map, const string& filename, const std::vector<Vector2d> *path) {
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

bool saveSlices(const dsl::Map<bool, 3> &cmap, string folder) {

  int slices = cmap.gs[0];
  int width  = cmap.gs[1];
  int height = cmap.gs[2];

  char data[width*height*3];//rgb channels

  if(folder.back() != '/')
    folder = folder+'/';

  string file = "map";
  string ext = ".ppm";
  string filename = folder+file+"0"+ext;

  std::fstream fs(filename, std::fstream::out);
  if(!fs.is_open()){
    cout<<"Couldn't write occupancy slices. Folder doesn't exist"<<endl;
    return false;
  }
  fs.close();

  dsl::Map<bool, 3>::SlicePtr pomap;

  for( int idx_a = 0; idx_a < slices; idx_a++){
    filename = folder+file+to_string(idx_a)+ext;
    fs.open(filename, std::fstream::out);
    assert(fs.is_open());
    if(!fs.is_open())
      continue;

    if(!pomap)
      pomap = cmap.GetSlice(idx_a,0);
    else
      cmap.GetSlice(*pomap,idx_a,0); //no reallocation of memory

    fs << "P6" << std::endl << width << " " << height << std::endl << "255" << std::endl;
    int ind = 0;
    for (int i = 0; i < width * height; i++, ind += 3) {
      data[ind] = data[ind + 1] = data[ind + 2] = (char)(pomap->cells[i] * 100);
    }
    assert(ind == 3*width*height);
    fs.write(data, ind);
    fs.close();
  }
  return true;
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

  for(size_t i=0; i<path.size(); i++){
    Vector3d axy = path[i];
    vector<Vector2d> vs,vs_grid; geom.GetTrueCorners(vs,axy(0));//Get corners of car relative to it's origin
    for_each(vs.begin(),vs.end(),[&](Vector2d& v){v +=axy.tail<2>();}); //corners relative to world origin

    smap.ToGridCoordinates(vs); //convert to grid coordinates
    int seq[]={3,0,1,2};
    double lwm = 0.05; //meters
    int lw = ceil(lwm/smap.cs[0]);//line width in pixels
    smap.cells = vector<bool>(smap.nc,false); //reset
    for(int i=0;i<4;i++)
      addLine<bool>(smap.cells, smap.gs[0], smap.gs[1], vs[i], vs[seq[i]], true,lw);

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

void saveMapWithPrims(const dsl::Map<bool, 2>& omap, std::string filename,
                      const std::vector<vector<Vector2d>>& prims,int scale){
  shared_ptr< dsl::Map<bool,2> > psmap = omap.ScaleUp(scale);
  dsl::Map<bool,2>& smap = *psmap; //scaled up for clarity

  //The image to be saved
  char data[smap.nc*3];
  int ind = 0;
  for (int i = 0; i < smap.nc; i++, ind += 3)
    data[ind] = data[ind + 1] = data[ind + 2] = (char)(smap.cells[i] * 100); //copy the map in grey

  for(auto& prim: prims){
    for(size_t i=0; i< prim.size(); i++){
      int id = smap.Id(prim[i]);
      if(i==0){
        data[3*id + 1] = 255; //start point in green
      }else if(i==prim.size()-1){
        data[3*id]     = 255; //goal point in red
      }else{
        data[3*id + 2] = 255; //point along the way in blue
      }
    }
  }

  std::fstream fs(filename, std::fstream::out);
  assert(fs.is_open());
  fs << "P6" << std::endl << smap.gs[0] << " " << smap.gs[1] << std::endl << "255" << std::endl;

  assert(ind == smap.nc*3);

  fs.write(data, ind);
  fs.close();
}

bool MakeSE2Map(const Map<bool, 2> &omap, Map<bool, 3> &cmap) {
  assert(omap.gs[0] == cmap.gs[1]);
  assert(omap.gs[1] == cmap.gs[2]);

  if(omap.gs[0] != cmap.gs[1] || omap.gs[1] != cmap.gs[2])
    return false;

  for (int i = 0; i < cmap.gs[0]; ++i) {
    for (int j = 0; j < cmap.gs[1]; ++j) {
      for (int k = 0; k < cmap.gs[2]; ++k) {
        int id2 = j + cmap.gs[1]*k; //2d index ito map
        assert(id2 < omap.nc);
        int id3 = i + cmap.gs[0]*j + cmap.gs[0]*cmap.gs[1]*k; // 3d index into cmap
        assert(id3 < cmap.nc);
        cmap.cells[id3] = omap.cells[id2];
      }
    }
  }
  return false;
}

//bool MakeSE2Map(const CarGeom& geom, const Map<bool, 2> &omap, Map<bool, 3> &cmap, int nthreads) {
//  assert(omap.gs[0] == cmap.gs[1]);
//  assert(omap.gs[1] == cmap.gs[2]);
//
//  vector<Vector2d> points;
//  geom.Raster(omap.cs, points);
//
//  Matrix2d R;
//
//  int dim_a(0); //dimension for angle
//  int dim_x(1); //dimension for angle
//  int dim_y(2); //dimension for angle
//
//  for (int idx_a = 0; idx_a < cmap.gs[dim_a]; ++idx_a) {
//    // dilate map for a particular angle
//    double theta = cmap.CellCenterIth(idx_a,dim_a);
//
//    // make a rotation matrix
//    double ct = cos(theta);
//    double st = sin(theta);
//    R(0,0) = ct; R(0,1) = -st;
//    R(1,0) = st; R(1,1) = ct;
//
//    for (int idx_x = 0; idx_x < cmap.gs[dim_x]; ++idx_x) {
//      double x = cmap.CellCenterIth(idx_x,dim_x);
//      for (int idx_y = 0; idx_y < cmap.gs[2]; ++idx_y) {
//        // index into workspace
//        int id_omap = omap.Id(Vector2i(idx_x,idx_y));
//        assert(id_omap < omap.nc);
//
//        if (!omap.cells[id_omap])
//          continue;        // if free continue
//
//        double y = cmap.CellCenterIth(idx_y,dim_y);
//
//        Vector2d p0(x,y); // position of car origin
//        for (auto&& dp : points) {
//          Vector2d p = p0 + R*dp; // point on the car
//          cmap.Set(Vector3d(theta, p[0], p[1]), true);
//        }
//      }
//    }
//  }
//  return true;
//}


bool MakeSE2Map(const CarGeom& geom, const Map<bool, 2> &omap, Map<bool, 3> &cmap, int nthreads) {
  assert(omap.gs[0] == cmap.gs[1]);
  assert(omap.gs[1] == cmap.gs[2]);
  assert(nthreads>0);

  if(omap.gs[0] != cmap.gs[1] || omap.gs[1] != cmap.gs[2])
    return false;

  nthreads = nthreads<1 ? 1:nthreads;
  vector< vector<bool> > dmaps(nthreads); //dilated map one for each thread
  for(auto& dmap:dmaps){dmap.resize(omap.nc);} //allocate memory for all threads

  int dim_a = 0; //The dimension corresponding to angle
  int n_a = cmap.gs[dim_a]; //number of different grid angles
  if(nthreads==1){
    for (int idx_a = 0; idx_a < n_a; ++idx_a) {
      // create a dilated map for a particular angle
      double theta = cmap.CellCenterIth(idx_a,dim_a);

      DilateMap(dmaps[0], omap, geom, theta); //Dilate map

      for (int idx_x = 0; idx_x < cmap.gs[1]; ++idx_x) {
        for (int idx_y = 0; idx_y < cmap.gs[2]; ++idx_y) {
          int id_omap = omap.Id( Vector2i(idx_x, idx_y) );
          int id_cmap = cmap.Id( Vector3i(idx_a, idx_x, idx_y) );
          cmap.cells[id_cmap] = dmaps[0][id_omap];
        }
      }
    }
  }else{
    std::vector<std::thread> threads(nthreads);

    for(size_t t = 0;t<nthreads;t++){
      threads[t] = std::thread(std::bind(
          [&](const int idx_a_start, const int idx_a_end, const int t)
          {            // loop over all items
        for(int idx_a = idx_a_start; idx_a< idx_a_end; idx_a++){
          // create a dilated map for a particular angle
          double theta = cmap.CellCenterIth(idx_a,dim_a);

          DilateMap(dmaps[t], omap, geom, theta); //Dilate map

          for (int idx_x = 0; idx_x < cmap.gs[1]; ++idx_x) {
            for (int idx_y = 0; idx_y < cmap.gs[2]; ++idx_y) {
              int id_omap = omap.Id( Vector2i(idx_x, idx_y) );
              int id_cmap = cmap.Id( Vector3i(idx_a, idx_x, idx_y) );
              cmap.cells[id_cmap] = dmaps[t][id_omap];
            }
          }
        }

          },t*n_a/nthreads,(t+1)==nthreads?n_a:(t+1)*n_a/nthreads,t));
    }
    std::for_each(threads.begin(),threads.end(),[](std::thread& x){x.join();});
  }
  return true;
}

void DilateMap(vector<bool>& dilated, const Map<bool,2>& omap, const CarGeom& geom, double theta){

  Matrix2x4d verts2d_rotd_pix;
  getRotdVertsInPixWrtOrg(verts2d_rotd_pix, geom.le(), geom.be(), geom.ox(), geom.oy(), omap.cs[0], omap.cs[1], theta);

  // round of the pixel values of the vertices above such that the rectange
  // formed by the rounded off
  //  vertices surrounds the rotated rectange
  Vector2i org2i_rotd_pix(0,0); // because it's wrt org itself and rounding doesn't matter
  Matrix2x4i verts2i_rotd_pix;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 4; j++)
      verts2i_rotd_pix(i, j) = verts2d_rotd_pix(i, j) > 0 ?
          ceil(verts2d_rotd_pix(i, j)) :
          floor(verts2d_rotd_pix(i, j));

  // Size of kernel is given by the horizontal rectange that bounds the rounded
  // off rectange above
  Vector2i verts2i_rotd_pix_min = verts2i_rotd_pix.rowwise().minCoeff();
  Vector2i verts2i_rotd_pix_max = verts2i_rotd_pix.rowwise().maxCoeff();
  Vector2i size2i_k = verts2i_rotd_pix_max - verts2i_rotd_pix_min;

  // Shift everything such that verts2i_rotd_pix_min is the [0,0] pixel of the
  // kernel
  Matrix2x4i verts2i_rotd_pospix =
      verts2i_rotd_pix.colwise() - verts2i_rotd_pix_min;
  Vector2i org2i_rotd_pospix = org2i_rotd_pix - verts2i_rotd_pix_min;

  // create dilation kernel by filling the inside of the rotated rectanges with
  // zero
  int w_k = size2i_k(0);
  int h_k = size2i_k(1);
  vector<bool> data_k(w_k * h_k);
  fillQuad<bool>(data_k,
                 size2i_k(0),
                 size2i_k(1),
                 verts2i_rotd_pospix.cast< double >(),
                 1.0);

  // Dilate
  dilate<bool>(dilated,
               omap.cells,
               omap.gs[0],
               omap.gs[1],
               data_k,
               size2i_k(0),
               size2i_k(1),
               org2i_rotd_pospix(0),
               org2i_rotd_pospix(1));
}

//void DilateMap(Map<bool,2>& dmap, const Map<bool,2> omap, const vector<Vector2d>& vertices){
//
//}


}
