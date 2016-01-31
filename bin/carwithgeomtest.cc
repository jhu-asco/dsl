// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Subhransu Mishra subhransu.kumar.mishra@gmail.com
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <string.h>
#include "gridsearch.h"
#include "cargrid.h"
#include "carcost.h"
#include "carconnectivity.h"
#include "utils.h"
#include "utilsimg.h"

#include <algorithm>
#include <cctype>
#include <iostream>

using namespace dsl;
using namespace std;
using namespace Eigen;

bool to_bool(std::string str) {
  std::transform(str.begin(), str.end(), str.begin(), ::tolower);
  std::istringstream is(str);
  bool b;
  is >> std::boolalpha >> b;
  return b;
}

int main(int argc, char** argv)
{
  if (argc!=4) {
    cout << "Usage: $./cartestdil map4.ppm true true" << endl;
    cout << "\t\t where map4.ppm is a map graphics file" << endl;
    cout << "\t\t The first bool is true/false depending on whether or not to use car geometry for planning" << endl;
    cout << "\t\t If the second bool is true, the car only moves forward(never reverses)" << endl;
    cout << "\t\t You can use only map4.ppm and map5.ppmfor this example" << endl;
    cout << "\t\t output will be written to graphics files path1.ppm path2.ppm and path3.ppm where the " << endl;
    cout << "\t\t\t path1.ppm is normal planning,  "<<endl;
    cout << "\t\t\t path2.ppm is replanning with start position moved to halfway point"<<endl;
    cout << "\t\t\t path3.ppm is replanning with start position moved to halfway point with a narrow passage along prev path closed off"<<endl;
    cout << "\t\t For this example the width and heigh of each pixel is 0.1m and the geometry of car is follows" << endl;
    cout << "\t\t Rectange of size 0.75(along x) and 0.43(along y) with the car origin being at 0.15m(along x) behind rect center" << endl;
    return 0;
  }

  //car geometry
  double l=0.75;
  double b=0.43;
  double ox = -0.15;
  double oy = 0.0;
  cout<<"The car geometry is given as a rectange of size 0.75(along x) and 0.43(along y) \n \t with the car origin being at 0.15m(along x) behind rect center"<<endl;

  //use geometry or not
  bool use_car_geom = to_bool(string(argv[2]));
  cout<<"use_geom:"<<(use_car_geom?"true":"false")<<endl;

  //use geometry or not
  bool only_fwd = to_bool(string(argv[3]));
  cout<<"car moves forward:"<<(only_fwd?"true":"false")<<endl;

  // what map is being used? only allow map4, map5 , map6 and map7
  // decide start goal based on the map type
  int len = strlen(argv[1]);
  int num_map = atoi(argv[1]+len-5);
  //for map4 and map5
  Vector3d goal, start;
  double sx=.1, sy=.1;// the size of pixel along x and y direction;
  if(num_map==4 || num_map==5)
  {
    start <<   0,  1, 11;
    goal  <<   0, 24, 12;
  }
  else
  {
    cout<<"map file is not map4 map5. Quitting"<<endl;
    return 1;
  }
  cout << "start=" << start.transpose() << endl;
  cout << "goal=" << goal.transpose() << endl;
  
  // load a map from ppm file and convert it into double*
  int width, height;
  char* chmap = load_map(width, height, argv[1]);
  double map_data[width*height];
  for (int i = 0; i < width*height; ++i)
    map_data[i] = 1000*(double)chmap[i];
  cout<<"Read the input map. Input image width and height"<<width<<","<<height<<endl;

  Map2d map(width, height, map_data);
  
  // create planner
  cout << "Creating a graph..." << endl;
  bool expand_at_start = true;
  struct timeval timer;
  timer_start(&timer);
  CarGrid* pgrid;
  CarGeom geom(l, b, ox, oy);
  if(use_car_geom)
    pgrid = new CarGrid(geom, map, sx, sy, M_PI/16, 1, 0.5);
  else 
    pgrid = new CarGrid(map, sx, sy, M_PI/16, 1,0.5);
  CarCost cost;
  double bp=1; //backward_penalty
  CarConnectivity connectivity(*pgrid,bp, only_fwd, 1);
  GridSearch<3, Matrix3d> search(*pgrid, connectivity, cost, expand_at_start);
  cout << "Created a graph with " << search.Vertices() << " vertices and " << search.Edges() << " edges. " << endl;
  printf("And the graph construction time= %ld  us\n", (long int)time);
  cout <<"Grid info:"<<endl;
  cout <<"\tcs:"<<pgrid->cs.transpose()<<endl;
  cout <<"\tgs:"<<pgrid->gs.transpose()<<endl;
  cout <<"\txlb:"<<pgrid->xlb.transpose()<<endl;
  cout <<"\txub:"<<pgrid->xub.transpose()<<endl;

  //Plan the path
  cout << "Planning a path..." << endl;
  SE2Path path, path2, path3;
  long time = timer_us(&timer);
  search.SetStart(start);
  search.SetGoal(goal);
  timer_start(&timer);
  search.Plan(path);
  time = timer_us(&timer);
  printf("plan path time= %ld  us\n", time);
  printf("path: edges=%lu len=%f\n", path.cells.size(), path.cost);
  cout << "After planning the graph now has " << search.Vertices() << " vertices and " << search.Edges() << " edges. " << endl;

  //   Print and save result for case 1
  int scale=3; //scale original image for clarity
  char mapPath[width*height*scale*scale];
  scaleMap<char>(mapPath, chmap,width,height,scale);
  vector<SE2Cell>::iterator it; 
  for (it = path.cells.begin(); it != path.cells.end(); ++it)
  {
    double theta = (it->c)(0);
    //    printf("(%d,%d) ",  path.cells[i].p[0], path.cells[i].p[1]);    
    int x = scale*pgrid->Index(it->c, 1);
    int y = scale*pgrid->Index(it->c, 2);
    Vector2d c(x,y);

    //draw the 4 lines
    Matrix2x4d vs;
    getRotdVertsInPixWrtOrg(vs, l, b, ox, oy, sx, sy, theta);
    vs*=scale;
    int seq2[]={3,0,1,2};
    int lw=1;//linewidth
    for(int i=0;i<4;i++)
      addLine<char>(mapPath,width*scale,height*scale, vs.col(i)+c, vs.col(seq2[i])+c, 2,lw);
  }
  save_map(mapPath, width*scale, height*scale, "path1.ppm");
  cout << "Map and path saved to path1.ppm" << endl;

  //Re-planning the planned path when the initial position is moved to a halfway point
  cout << "Re-planning the planned path when the initial position is moved to a halfway point" << endl;
  Vector3d cell_mid = path.cells[path.cells.size()/2].c;  // follow path until middle
  search.SetStart(cell_mid);

  timer_start(&timer);
  search.Plan(path2);
  time = timer_us(&timer);
  printf("re-planning path time= %ld  us\n", time);
  printf("path: edges=%lu len=%f\n", path.cells.size(), path2.cost);
  cout << "After re-planning the graph now has " << search.Vertices() << " vertices and " << search.Edges() << " edges. " << endl;

  //   Print and save result for case 2
  scaleMap<char>(mapPath, chmap,width,height,scale);
  for (it = path2.cells.begin(); it != path2.cells.end(); ++it)
  {
    double theta = (it->c)(0);
    //    printf("(%d,%d) ",  path.cells[i].p[0], path.cells[i].p[1]);
    int x = scale*pgrid->Index(it->c, 1);
    int y = scale*pgrid->Index(it->c, 2);
    Vector2d c(x,y);

    //draw the 4 lines
    Matrix2x4d vs;
    getRotdVertsInPixWrtOrg(vs, l, b, ox, oy, sx, sy, theta);
    vs*=scale;
    int seq2[]={3,0,1,2};
    int lw=1;//linewidth
    for(int i=0;i<4;i++)
      addLine<char>(mapPath,width*scale,height*scale, vs.col(i)+c, vs.col(seq2[i])+c, 2,lw);
  }
  save_map(mapPath, width*scale, height*scale, "path2.ppm");
  cout << "Map and path saved to path2.ppm" << endl;

  //Re-planning the planned path when the initial position is moved to a halfway point
  //  and a narrow passage along the prev path closed
  cout << "Re-planning the planned path when the initial position is moved to a halfway point and a narrow path along prev path closed off" << endl;

  //  cell_mid = path.cells[path.cells.size()/2].c;  // follow path until middle
  search.SetStart(cell_mid);

  for(int r=150; r<160;r++){
    for(int c=176;c<195;c++){
      int idx=c+r*width;
      chmap[idx]=2; //for displaying where we block the path
      for (double a = pgrid->xlb[0]+pgrid->cs[0]/2; a < pgrid->xub[0]; a += pgrid->cs[0]){
        Vector3d x(a, (c+0.5)*pgrid->cs[1], (r+0.5)*pgrid->cs[2]);
        // search.RemoveCell(x);
        search.SetCost(x, 1000);
      }
    }
  }

  for(int r=125; r<135;r++){
    for(int c=210;c<220;c++){
      int idx=c+r*width;
      chmap[idx]=2;//for displaying where we block the path
      for (double a = pgrid->xlb[0]+pgrid->cs[0]/2; a < pgrid->xub[0]; a += pgrid->cs[0]){
        Vector3d x(a,(c+0.5)*pgrid->cs[1],(r+0.5)*pgrid->cs[2]);
        // search.RemoveCell(x);
        search.SetCost(x, 1000);
      }
    }
  }
  
  timer_start(&timer);
  search.Plan(path3);
  time = timer_us(&timer);
  printf("re-planning path time= %ld  us\n", time);
  printf("path: edges=%lu len=%f\n", path3.cells.size(), path3.cost);
  cout << "After re-planning with narrow path closing the graph now has " << search.Vertices() << " vertices and " << search.Edges() << " edges. " << endl;

  //   Print and save result for case 3
  scaleMap<char>(mapPath, chmap,width,height,scale);
  for (it = path3.cells.begin(); it != path3.cells.end(); ++it)
  {
    double theta = (it->c)(0);
    //    printf("(%d,%d) ",  path.cells[i].p[0], path.cells[i].p[1]);
    int x = scale*pgrid->Index(it->c, 1);
    int y = scale*pgrid->Index(it->c, 2);
    Vector2d c(x,y);

    //draw the 4 lines
    Matrix2x4d vs;
    getRotdVertsInPixWrtOrg(vs, l, b, ox, oy, sx, sy, theta);
    vs*=scale;
    int seq2[]={3,0,1,2};
    int lw=1;//linewidth
    for(int i=0;i<4;i++)
      addLine<char>(mapPath,width*scale,height*scale, vs.col(i)+c, vs.col(seq2[i])+c, 2,lw);
  }
  save_map(mapPath, width*scale, height*scale, "path3.ppm");
  cout << "Map and path saved to path3.ppm" << endl;


  return 0;
}
