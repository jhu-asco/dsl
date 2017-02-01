#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "gridutils.h"
#include <numeric>
#include "params.h"
#include "ppm_reader.h"
#include <thread>
#include "utilsimg.h"

namespace dsl {

using namespace Eigen;
using namespace std;
using Vector2b = Matrix<bool,2,1>;

Map<bool, 2>::Ptr LoadPpm(const string& filename, const Vector2d &cs) {
  if(filename.compare(filename.size() - 4, 4, ".ppm")){
    cout<<"File doesn't have .ppm extension"<<endl;
    return nullptr;
  }
  ImageRGB img;
  if(!LoadPpm(img,filename)){
    cout<<"Had problems reading the image file"<<endl;
    return nullptr;
  }

  vector<bool> cells(img.h*img.w);
  for (int id = 0; id < img.h*img.w; id++)
    cells[id] = img.rdata[id]? 1:0;

  dsl::Map<bool, 2>::Ptr omap(new Map<bool, 2>(Vector2d(0,0), Vector2d(cs[0]*img.w, cs[1]*img.h), Vector2i(img.w,img.h)));
  omap->set_cells(cells);

  return omap;
}

bool SavePpm(const dsl::Map<bool, 2> &map, const string& filename) {
  if(filename.compare(filename.size() - 4, 4, ".ppm")){
    cout<<"File doesn't have .ppm extension"<<endl;
    return false;
  }

  //The image to be saved
  ImageRGB img;
  img.w = map.gs()[0];
  img.h = map.gs()[1];
  img.bitdepth = ImageRGB::BD8; // only occupancy information. Higher bit depth not required.
  img.resize(img.w*img.h);

  for (int id = 0; id < map.nc(); id++){
    img.rdata[id] = map.Get(id)? img.bitdepth : 0;
    img.gdata[id] = map.Get(id)? img.bitdepth : 0;
    img.bdata[id] = map.Get(id)? img.bitdepth : 0;
  }

  if(!SavePpm(img,filename)){
    cout<<"Had problems saving the image file"<<endl;
    return false;
  }
  return true;
}

bool SavePpm(const dsl::Map<TerrainData, 2> &tmap, const string& filename) {
  if(filename.compare(filename.size() - 4, 4, ".ppm")){
    cout<<"File doesn't have .ppm extension"<<endl;
    return false;
  }

  //save the image file
  ImageRGB img;
  img.w = tmap.gs()[0];
  img.h = tmap.gs()[1];
  img.bitdepth = ImageRGB::BD16;
  img.resize(img.w*img.h);

  double maxh = numeric_limits<double_t>::lowest(); //max height
  double maxt = numeric_limits<double_t>::lowest(); //max traversibility
  for (int id = 0; id < tmap.nc(); id++){
    maxh = tmap.Get(id).height > maxh ? tmap.Get(id).height: maxh;
    maxt = tmap.Get(id).traversibility > maxt ? tmap.Get(id).traversibility: maxt;
  }
  double hscale = img.bitdepth/maxh;
  double tscale = img.bitdepth/maxt;

  for (int id = 0; id < tmap.nc(); id++){
    img.rdata[id] = 0.5*tmap.Get(id).height*hscale;
    img.gdata[id] = 0.5*tmap.Get(id).traversibility*tscale;

  }

  if(!SavePpm(img,filename)){
    cout<<"Had problems saving the image file"<<endl;
    return false;
  }

  return true;
}

bool SavePpm(const dsl::Map<bool, 3> &cmap, string folder) {

  int slices = cmap.gs()[0];
  int width  = cmap.gs()[1];
  int height = cmap.gs()[2];

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

  dsl::Map<bool, 3>::SlicePtr omap;

  for( int idx_a = 0; idx_a < slices; idx_a++){
    filename = folder+file+to_string(idx_a)+ext;
    fs.open(filename, std::fstream::out);
    assert(fs.is_open());
    if(!fs.is_open())
      continue;

    if(!omap)
      omap = cmap.GetSlice(idx_a, 0);
    else
      cmap.GetSlice(omap.get( ), idx_a, 0); //no reallocation of memory

    fs << "P6" << std::endl << width << " " << height << std::endl << "255" << std::endl;
    int ind = 0;
    for (int i = 0; i < width * height; i++, ind += 3) {
      data[ind] = data[ind + 1] = data[ind + 2] = (char)(omap->Get(i) * 100);
    }
    assert(ind == 3*width*height);
    fs.write(data, ind);
    fs.close();
  }
  return true;
}

Map<bool, 2>::Ptr LoadOmap(const string& omapfile){
  if(omapfile.compare(omapfile.size() - 5, 5, ".omap")){
    cout<<"File doesn't have .omap extension"<<endl;
    return false;
  }

  Params params(omapfile.c_str());

  Vector2d cs;
  if(!params.GetVector2d("cellsize",cs)){
    cout<<"cellsize parameter missing in .omap file"<<endl;
    return false;
  }

  string imgfilerel;
  if(!params.GetString("image",imgfilerel)){
    cout<<"image parameter missing in .omap file"<<endl;
    return false;
  }
  string imgfile = omapfile.substr(0, omapfile.find_last_of("\\/")) + "/" + imgfilerel;

  ImageRGB img;
  if(!LoadPpm(img,imgfile)){
    cout<<"Had problems reading the image file"<<endl;
    return false;
  }

  vector<bool> cells(img.h*img.w);
  for (int id = 0; id < img.h*img.w; id++)
    cells[id] = img.rdata[id]? 1:0;

  dsl::Map<bool, 2>::Ptr omap(new Map<bool, 2>(Vector2d(0,0), Vector2d(cs[0]*img.w, cs[1]*img.h), Vector2i(img.w,img.h)));
  omap->set_cells(cells);
  return omap;
}

Map<TerrainData, 2>::Ptr LoadTmap(const string& tmapfile){
  if(tmapfile.compare(tmapfile.size() - 5, 5, ".tmap")){
    cout<<"extention doesn't match .tmap. The file name is: "<<tmapfile<<endl;
    return false;
  }
  Params params(tmapfile.c_str());

   Vector2d cs;
   if(!params.GetVector2d("cellsize",cs)){
     cout<<"param cellsize not found in .tmap file"<<endl;
     return false;
   }
   if(!(cs(0)>0 && cs(1)>0)){
     cout<<"cellsize should be positive"<<endl;
     return false;
   }

   string imgfilerel;
   if(!params.GetString("image",imgfilerel)){
     cout<<"param image not found in .tmap file"<<endl;
     return false;
   }
   string imgfile = tmapfile.substr(0, tmapfile.find_last_of("\\/")) + "/" + imgfilerel;

   double hscale;
   if(!params.GetDouble("hscale",hscale)){
     cout<<"param hscale not found in .tmap file"<<endl;
     return false;
   }
   double tscale;
   if(!params.GetDouble("tscale",tscale)){
     cout<<"param tscale not found in .tmap file"<<endl;
     return false;
   }

   ImageRGB img;
   if(!LoadPpm(img,imgfile)){
     cout<<"Had problems reading the image file:"<< imgfile<<endl;
     return false;
   }

   vector<TerrainData> cells(img.h*img.w);
   for (int id = 0; id < img.h*img.w; id++){
     cells[id].height = img.rdata[id]*hscale;
     cells[id].traversibility = 1 + img.gdata[id]*tscale;
   }
   Map<TerrainData, 2>::Ptr tmap(new Map<TerrainData, 2>(Vector2d(0,0), Vector2d(cs[0]*img.w, cs[1]*img.h), Vector2i(img.w,img.h)));
   tmap->set_cells(cells);
   return tmap;
}

bool saveOmap(Map<bool, 2>& omap, const string& omapfile){
  if(omapfile.compare(omapfile.size() - 5, 5, ".omap")){
    cout<<"extention doesn't match .omap. The file name is: "<<omapfile<<endl;
    return false;
  }

  //save the image file
  std::string map_img_filename = ReplaceExtension(omapfile, ".ppm");
  if(!SavePpm(omap, map_img_filename))
    return false;

  //save the text file
  ofstream file(omapfile, std::ofstream::out);
  if(!file.is_open())
    return false;
  file << "image= "<<map_img_filename<<endl;
  file << "cellsize= "<<omap.cs()[0]<<", "<<omap.cs()[1]<<endl;
  file.close();

  return true;
}

bool saveTmap(Map<TerrainData, 2>& tmap, const string& tmapfile){
  if(tmapfile.compare(tmapfile.size() - 5, 5, ".tmap")){
    cout<<"extention doesn't match .tmap. The file name is: "<<tmapfile<<endl;
    return false;
  }

  //save the image file
  std::string map_img_filename = ReplaceExtension(tmapfile, ".ppm");
  ImageRGB img;
  img.w = tmap.gs()[0];
  img.h = tmap.gs()[1];
  img.bitdepth = ImageRGB::BD16;
  img.resize(img.w*img.h);

  double maxh = numeric_limits<double_t>::lowest(); //max height
  double maxt = numeric_limits<double_t>::lowest(); //max traversibility
  for (int id = 0; id < tmap.nc(); id++){
    maxh = tmap.Get(id).height > maxh ? tmap.Get(id).height: maxh;
    maxt = tmap.Get(id).traversibility > maxt ? tmap.Get(id).traversibility: maxt;
  }
  double hscale = img.bitdepth/maxh;
  double tscale = img.bitdepth/maxt;

  for (int id = 0; id < tmap.nc(); id++){
    img.rdata[id] = 0.5*tmap.Get(id).height*hscale;
    img.gdata[id] = 0.5*tmap.Get(id).traversibility*tscale;

  }

  if(!SavePpm(img,map_img_filename)){
    cout<<"Had problems saving the image file"<<endl;
    return false;
  }

  //save the text file
  ofstream file(tmapfile, std::ofstream::out);
  if(!file.is_open())
    return false;
  file << "image= "<<map_img_filename<<endl;
  file << "cellsize= "<<tmap.cs()[0]<<", "<<tmap.cs()[1]<<endl;
  file << "hscale= "<<hscale<<endl;
  file << "tscale= "<<tscale<<endl;
  file.close();

  return true;
}

bool SavePpmWithPath(const dsl::Map<bool, 2>& omap, std::string filename,
                      const std::vector<Vector3d>& path, int scale, const CarGeom* geom){
  //checks
  if(scale<1){
    cout<<"Scale should be >=1"<<endl;
    return false;
  }

  if(filename.compare(filename.size() - 4, 4, ".ppm")){
    cout<<"File doesn't have .ppm extension"<<endl;
    return false;
  }

  dsl::Map<bool,2>::Ptr smap = omap.ScaleUp(scale);

  //The image to be saved
  ImageRGB img;
  img.w = smap->gs()[0];
  img.h = smap->gs()[1];
  img.bitdepth = ImageRGB::BD8;
  img.resize(img.w*img.h);

  for (int id = 0; id < smap->nc(); id++)
    img.set_rgb(id, smap->Get(id)? 100:0);

  for(size_t i=0; i<path.size(); i++){
    if(geom){
      Vector3d axy = path[i];
      vector<Vector2d> vs; geom->GetTrueCorners(vs,axy(0));//Get corners of car relative to it's origin
      for_each(vs.begin(),vs.end(),[&](Vector2d& v){v +=axy.tail<2>();}); //corners relative to world origin

      smap->ToGridCoordinates(&vs); //convert to grid coordinates
      int seq[]={3,0,1,2};
      double lwm = 0.05; //line width meters
      int lw = ceil(lwm/smap->cs()[0]);//line width in pixels
      vector<bool> temp(smap->nc(), false);
      for(int i=0;i<4;i++){
        addLine<bool>(temp, smap->gs()[0], smap->gs()[1], vs[i], vs[seq[i]], true,lw);
      }

      if(i==0){
        for (int id = 0; id < smap->nc(); id++)
          if(temp[id])
             img.set_to_green(id,img.bitdepth);//Start of path in green
      }else if(i==path.size()-1){
        for (int id = 0; id < smap->nc(); id++)
          if(temp[id])
             img.set_to_red(id, img.bitdepth);//End of path in red
      }else{
        for (int id = 0; id < smap->nc(); id++)
          if(temp[id])
             img.set_to_blue(id, img.bitdepth);//Points in between in blue
      }
    }else{
      Vector2d pos = path.at(i).tail<2>();
      int id = smap->Id(pos);
      if(i==0)
         img.set_to_green(id,img.bitdepth);//Start of path in green
      else if(i==path.size()-1)
         img.set_to_red(id, img.bitdepth);//End of path in red
      else
         img.set_to_blue(id, img.bitdepth);//Points in between in blue

    }
  }

  if(!SavePpm(img,filename)){
    cout<<"Had problems saving the image file"<<endl;
    return false;
  }

  return true;
}

//TODO: requires testing
bool SavePpmWithPath(const dsl::Map<TerrainData, 2>& tmap, std::string filename,
                     const std::vector<Vector3d>& path, int scale, const CarGeom* geom){
  //checks
  if(scale<1){
    cout<<"Scale should be >=1"<<endl;
    return false;
  }

  if(filename.compare(filename.size() - 4, 4, ".ppm")){
    cout<<"File doesn't have .ppm extension"<<endl;
    return false;
  }

  dsl::Map<TerrainData,2>::Ptr smap = tmap.ScaleUp(scale);

  //The image to be saved
  ImageRGB img;
  img.w = smap->gs()[0];
  img.h = smap->gs()[1];
  img.bitdepth = ImageRGB::BD16;
  img.resize(img.w*img.h);

  double maxh = numeric_limits<double_t>::lowest(); //max height
  double maxt = numeric_limits<double_t>::lowest(); //max traversibility
  for (int id = 0; id < smap->nc(); id++){
    maxh = smap->Get(id).height > maxh ? smap->Get(id).height: maxh;
    maxt = smap->Get(id).traversibility > maxt ? smap->Get(id).traversibility: maxt;
  }
  double hscale = img.bitdepth/maxh;
  double tscale = img.bitdepth/maxt;

  for (int id = 0; id < smap->nc(); id++){
    img.rdata[id] = 0.5*smap->Get(id).height*hscale;
    img.gdata[id] = 0.5*smap->Get(id).traversibility*tscale;
  }

  for(size_t i=0; i<path.size(); i++){
    if(geom){
      Vector3d axy = path[i];
      vector<Vector2d> vs; geom->GetTrueCorners(vs,axy(0));//Get corners of car relative to it's origin
      for_each(vs.begin(),vs.end(),[&](Vector2d& v){v +=axy.tail<2>();}); //corners relative to world origin

      smap->ToGridCoordinates(&vs); //convert to grid coordinates
      int seq[]={3,0,1,2};
      double lwm = 0.05; //meters
      int lw = ceil(lwm/smap->cs()[0]);//line width in pixels
      vector<bool> temp(smap->nc(), false);
      for(int i=0;i<4;i++){
        addLine<bool>(temp, smap->gs()[0], smap->gs()[1], vs[i], vs[seq[i]], true,lw);
      }

      if(i==0){
        for (int id = 0; id < smap->nc(); id++)
          if(temp[id])
             img.set_to_green(id,img.bitdepth);//Start of path in green

      }else if(i==path.size()-1){
        for (int id = 0; id < smap->nc(); id++)
          if(temp[id])
             img.set_to_red(id, img.bitdepth);//End of path in red
      }else{
        for (int id = 0; id < smap->nc(); id++)
          if(temp[id])
             img.set_to_blue(id, img.bitdepth);//Middle points of path in blue
      }
    }else{
      Vector2d pos = path.at(i).tail<2>();
      int id = smap->Id(pos);
      if(i==0)
         img.set_to_green(id,img.bitdepth);//Start of path in green
      else if(i==path.size()-1)
         img.set_to_red(id, img.bitdepth);//End of path in red
      else
         img.set_to_blue(id, img.bitdepth);//Points in between in blue

    }
  }

  if(!SavePpm(img,filename)){
    cout<<"Had problems saving the image file"<<endl;
    return false;
  }
  return true;
}

//TODO: requires testing
bool SavePpmWithPrimitives(const dsl::Map<bool, 2>& omap, std::string filename,
                      const std::vector<vector<Vector2d>>& prims, int scale){
  //checks
  if(scale<1){
    cout<<"Scale should be >=1"<<endl;
    return false;
  }

  if(filename.compare(filename.size() - 4, 4, ".ppm")){
    cout<<"File doesn't have .ppm extension"<<endl;
    return false;
  }

  dsl::Map<bool,2>::Ptr smap = omap.ScaleUp(scale);

  //The image to be saved
  ImageRGB img;
  img.w = smap->gs()[0];
  img.h = smap->gs()[1];
  img.bitdepth = ImageRGB::BD8;
  img.resize(img.w*img.h);

  for (int id = 0; id < smap->nc(); id++)
    img.set_rgb(id, smap->Get(id)? 100:0);

  for(auto& prim: prims){
    for(size_t i=0; i< prim.size(); i++){
      int id = smap->Id(prim[i]);
      if(i==0)
         img.set_to_green(id,img.bitdepth); //start point in green
      else if(i==prim.size()-1)
         img.set_to_red(id, img.bitdepth); //goal point in red
      else
         img.set_to_blue(id, img.bitdepth); //point along the way in blue
    }
  }

  if(!SavePpm(img,filename)){
    cout<<"Had problems saving the image file"<<endl;
    return false;
  }
  return true;
}

bool SavePpmWithPrimitives(const dsl::Map<TerrainData, 2>& tmap, std::string filename,
                      const std::vector<vector<Vector2d>>& prims, int scale){
  //checks
  if(scale<1){
    cout<<"Scale should be >=1"<<endl;
    return false;
  }

  if(filename.compare(filename.size() - 4, 4, ".ppm")){
    cout<<"File doesn't have .ppm extension"<<endl;
    return false;
  }

  shared_ptr< dsl::Map<TerrainData,2> > smap = tmap.ScaleUp(scale);

  //The image to be saved
  ImageRGB img;
  img.w = smap->gs()[0];
  img.h = smap->gs()[1];
  img.bitdepth = ImageRGB::BD16;
  img.resize(img.w*img.h);

  double maxh = numeric_limits<double_t>::lowest(); //max height
  double maxt = numeric_limits<double_t>::lowest(); //max traversibility
  for (int id = 0; id < smap->nc(); id++){
    maxh = smap->Get(id).height > maxh ? smap->Get(id).height: maxh;
    maxt = smap->Get(id).traversibility > maxt ? smap->Get(id).traversibility: maxt;
  }
  double hscale = img.bitdepth/maxh;
  double tscale = img.bitdepth/maxt;
  for (int id = 0; id < smap->nc(); id++){
    img.rdata[id] = 0.5*smap->Get(id).height*hscale;
    img.gdata[id] = 0.5*smap->Get(id).traversibility*tscale;
  }

  char data[smap->nc()*3];

  for(auto& prim: prims){
    for(size_t i=0; i< prim.size(); i++){
      int id = smap->Id(prim[i]);
      if(i==0)
         img.set_to_green(id,img.bitdepth); //start point in green
      else if(i==prim.size()-1)
         img.set_to_red(id, img.bitdepth); //goal point in red
      else
         img.set_to_blue(id, img.bitdepth); //point along the way in blue
    }
  }

  if(!SavePpm(img,filename)){
    cout<<"Had problems saving the image file"<<endl;
    return false;
  }
  return true;
}

Map<bool, 3>::Ptr MakeCmap(const Map<bool, 2>& omap, double csa) {

  Map<bool, 3>::Ptr cmap = omap.GetStack(0,-M_PI, M_PI, csa, true);

  auto fun = [&](int id, const Vector3i& gidx){
    Vector2i gidx_omap = gidx.tail<2>();
    cmap->Set(id, omap.Get(gidx_omap) );
  };
  cmap->LoopOver(fun);

  return cmap;
}

// //rasterization based method
//Map<bool, 3>::Ptr MakeCmap(const Map<bool, 2> &omap, double csa, const CarGeom& geom, int nthreads) {
//  Map<bool, 3>::Ptr cmap = omap.GetStack(0,-M_PI, M_PI, csa, true);
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
//  for (int idx_a = 0; idx_a < cmap->gs()[dim_a]; ++idx_a) {
//    // dilate map for a particular angle
//    double theta = cmap->CellCenterIth(idx_a,dim_a);
//
//    // make a rotation matrix
//    double ct = cos(theta);
//    double st = sin(theta);
//    R(0,0) = ct; R(0,1) = -st;
//    R(1,0) = st; R(1,1) = ct;
//
//    for (int idx_x = 0; idx_x < cmap->gs()[dim_x]; ++idx_x) {
//      double x = cmap->CellCenterIth(idx_x,dim_x);
//      for (int idx_y = 0; idx_y < cmap->gs()[2]; ++idx_y) {
//        // index into workspace
//        int id_omap = omap.Id(Vector2i(idx_x,idx_y));
//        assert(id_omap < omap.nc());
//
//        if (!omap.cells_[id_omap])
//          continue;        // if free continue
//
//        double y = cmap->CellCenterIth(idx_y,dim_y);
//
//        Vector2d p0(x,y); // position of car origin
//        for (auto&& dp : points) {
//          Vector2d p = p0 + R*dp; // point on the car
//          cmap->Set(Vector3d(theta, p[0], p[1]), true);
//        }
//      }
//    }
//  }
//  return cmap;
//}


Map<bool, 3>::Ptr MakeCmap(const Map<bool, 2>& omap, double csa, const CarGeom& geom, int nthreads) {

  Map<bool, 3>::Ptr cmap = omap.GetStack(0,-M_PI, M_PI, csa, true);

  nthreads = nthreads<1 ? 1:nthreads;
  vector< vector<bool> > dmaps(nthreads); //dilated map one for each thread
  for(auto& dmap:dmaps){dmap.resize(omap.nc());} //allocate memory for all threads

  int dim_a = 0; //The dimension corresponding to angle
  int n_a = cmap->gs()[dim_a]; //number of different grid angles
  if(nthreads==1){
    for (int idx_a = 0; idx_a < n_a; ++idx_a) {
      // create a dilated map for a particular angle
      double theta = cmap->CellCenterIth(idx_a,dim_a);

      DilateMap(dmaps[0], omap, geom, theta); //Dilate map

      for (int idx_x = 0; idx_x < cmap->gs()[1]; ++idx_x) {
        for (int idx_y = 0; idx_y < cmap->gs()[2]; ++idx_y) {
          int id_omap = omap.Id( Vector2i(idx_x, idx_y) );
          int id_cmap = cmap->Id( Vector3i(idx_a, idx_x, idx_y) );
          cmap->Set(id_cmap,dmaps[0][id_omap]);
        }
      }
    }
  }else{
    std::vector<std::thread> threads(nthreads);

    for(int t = 0;t<nthreads;t++){
      threads[t] = std::thread(std::bind(
          [&](const int idx_a_start, const int idx_a_end, const int t)
          {            // loop over all items
        for(int idx_a = idx_a_start; idx_a< idx_a_end; idx_a++){
          // create a dilated map for a particular angle
          double theta = cmap->CellCenterIth(idx_a,dim_a);

          DilateMap(dmaps[t], omap, geom, theta); //Dilate map

          for (int idx_x = 0; idx_x < cmap->gs()[1]; ++idx_x) {
            for (int idx_y = 0; idx_y < cmap->gs()[2]; ++idx_y) {
              int id_omap = omap.Id( Vector2i(idx_x, idx_y) );
              int id_cmap = cmap->Id( Vector3i(idx_a, idx_x, idx_y) );
              cmap->Set(id_cmap,dmaps[t][id_omap]);
            }
          }
        }

          },t*n_a/nthreads,(t+1)==nthreads?n_a:(t+1)*n_a/nthreads,t));
    }
    std::for_each(threads.begin(),threads.end(),[](std::thread& x){x.join();});
  }
  return cmap;
}

Map<bool, 3>::Ptr MakeCmap(const Map<TerrainData, 2> & tmap, double csa){
  Map<bool, 2> omap(tmap.xlb(), tmap.xub(), tmap.gs());

  //Iterate over all cells
  for(int id = 0; id < omap.nc(); id++)
  omap.Set(id, std::isnan(tmap.Get(id).traversibility));//occupied

  return MakeCmap(omap,csa);
}

Map<bool, 3>::Ptr MakeCmap(const Map<TerrainData, 2>& tmap, double csa, const CarGeom& geom, int nthreads){
  Map<bool, 2> omap(tmap.xlb(), tmap.xub(), tmap.gs());

    //Iterate over all cells
    for(int id = 0; id < omap.nc(); id++)
      omap.Set(id, std::isnan(tmap.Get(id).traversibility));//occupied

    return MakeCmap(omap,csa,geom,nthreads);
}

void DilateMap(vector<bool>& dilated, const Map<bool,2>& omap, const CarGeom& geom, double theta){

  Matrix2x4d verts2d_rotd_pix;
  getRotdVertsInPixWrtOrg(verts2d_rotd_pix, geom.le(), geom.be(), geom.ox(), geom.oy(), omap.cs()[0], omap.cs()[1], theta);

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
               omap.cells(),
               omap.gs()[0],
               omap.gs()[1],
               data_k,
               size2i_k(0),
               size2i_k(1),
               org2i_rotd_pospix(0),
               org2i_rotd_pospix(1));
}

}
