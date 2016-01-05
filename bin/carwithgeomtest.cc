#include <string.h>
#include "gridsearch.h"
#include "cargrid.h"
#include "carcost.h"
#include "carconnectivity.h"
#include "utils.h"

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
  if (argc!=3) {
    cout << "Usage: $./cartestdil map4.ppm true" << endl;
    cout << "\t\t where map4.ppm is a map graphics file" << endl;
    cout << "\t\t last param is true/false depending on whether or not to use car geometry for planning" << endl;
    cout << "\t\t You can use only map4.ppm, map5.ppm, map6.ppm and map7.ppm for this example" << endl;
    cout << "\t\t output will be written to graphics files path1.ppm and path2.ppm" << endl;
    cout << "\t\t For this example the width and heigh of each pixel is 0.1m and the geometry of car is follows" << endl;
    cout << "\t\t Rectange of size 0.75(along x) and 0.43(along y) with the car origin being at 0.15m(along x) behind rect center" << endl;
    return 0;
  }

  //car geometry
  double l=0.75;
  double b=0.43;
  double ox = -0.15;
  double oy = 0.0;

  //use geometry or not
  bool use_car_geom = to_bool(string(argv[2]));
  cout<<"use_geom:"<<(use_car_geom?"true":"false")<<endl;
  //what map is being used?
  int len = strlen(argv[1]);
  cout<<"strlen:"<<len<<endl;
  int num_map = atoi(argv[1]+len-5);
  cout<<"number:"<<argv[1]<<endl;
  cout<<"num_map:"<<num_map<<endl;

  //for map4 and map5
  Vector3d goal, start;
  if(num_map==4 || num_map==5)
  {
    goal <<      0,  1, 11;
    start<< M_PI/2, 24, 12;
  }
  else if(num_map==6 || num_map==7)
  {
    goal << M_PI/2, 8, 5;
    start<< M_PI/2, 2, 0.5;
  }
  else
  {
    cout<<"map file is not map4 map5 or map6. Quitting"<<endl;
    return 1;
  }

  int width, height;
  // load a map from ppm file
  char* chmap = load_map(&width, &height, argv[1]);


  cout<<"width and height"<<width<<","<<height<<endl;


  double map[width*height];
  for (int i = 0; i < width*height; ++i)
    map[i] = 1000*(double)chmap[i];


  cout << "Creating a graph..." << endl;

  // create planner
  struct timeval timer;
  timer_start(&timer);

  double sx=0.1, sy=0.1;// the size of pixel along x and y direction;

  CarGrid* pgrid;
  if(use_car_geom)
    pgrid = new CarGrid(l,b,ox,oy,width, height, map, sx, sy, M_PI/17, 1, 0.5);
  else
    pgrid = new CarGrid(width, height, map, sx, sy, M_PI/17, 1,0.5);


  CarCost cost;
  CarConnectivity connectivity(*pgrid);
  GridSearch<3, Matrix3d> search(*pgrid, connectivity, cost, false);
  SE2Path path, optPath;

  long time = timer_us(&timer);
  printf("graph construction time= %ld  us\n", time);

  search.SetStart(start);
  search.SetGoal(goal);
  //  search.SetGoal(Vector3d(grid.xub[0]*.5, grid.xub[1]*.58, 15.0/16*M_PI));

  cout << "Created a graph with " << search.Vertices() << " vertices and " << search.Edges() << " edges. " << endl;

  cout << "Planning a path..." << endl;
  // plan
  timer_start(&timer);
  search.Plan(path);
  time = timer_us(&timer);
  printf("plan path time= %ld  us\n", time);
  printf("path: edges=%lu len=%f\n", path.cells.size(), path.cost);

  cout << "Graph has " << search.Vertices() << " vertices and " << search.Edges() << " edges. " << endl;

  // Print and save result
  int scale=5; //scale original image for clarity
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
    //mapPath[y*width*scale + x] = 2;// path.cells[i].p[1]*width +  path.cells[i].p[0]] = 2;
    //draw the 4 lines
    Matrix2x4d vs;

    getRotdVertsInPixWrtOrg(vs, l, b, ox, oy, sx, sy, theta);
    vs*=scale;
    int seq2[]={3,0,1,2};
    int lw=1;//linewidth
    for(int i=0;i<4;i++)
      addLine<char>(mapPath,width*scale,height*scale, vs.col(i)+c, vs.col(seq2[i])+c, 2,lw);
  }

  // save it to image for viewing
  save_map(mapPath, width*scale, height*scale, "path1.ppm");
  cout << "Map and path saved to path1.ppm" << endl;

  return 1;

  // follow path until middle
  Vector3d c = path.cells[path.cells.size()/2].c;
  search.SetStart(c);

  // simulate closing the narrow passage
  if (0) {
    for (double a = pgrid->xlb[0]+pgrid->cs[0]/2; a < pgrid->xub[0] + pgrid->cs[0]/2; a += pgrid->cs[0]) {
      // by increasing the cost drastically
      //      search.SetCost(Vector3d(29,18,a), 1000);
      //      search.SetCost(Vector3d(30,18,a), 1000);
      //      search.SetCost(Vector3d(31,18,a), 1000);
    }
  } else {
    // a better way: by simply removing the passage
    for (double y = c[1] - 2; y < c[1] +2; ++y)  {
      for (double a = pgrid->xlb[0]+pgrid->cs[0]/2; a < pgrid->xub[0] + pgrid->cs[0]/2; a += pgrid->cs[0]) {
        search.RemoveCell(Vector3d(a, c[0],y));
      }
    }
  }

  // this is just for display
  memcpy(mapPath, chmap, width*height);
  mapPath[18*width + 29] = 1;  mapPath[18*width + 30] = 1;  mapPath[18*width + 31] = 1;

  // replan
  timer_start(&timer);
  search.Plan(path);
  time = timer_us(&timer);
  printf("replan path time= %ld us\n", time);
  printf("path: count=%lu len=%f\n", path.cells.size(), path.cost);
  fflush(stdout);


  // bypass the old vertex
  //  gdsl.AddEdge(29,18,31,18);  
  // // replan
  // timer_start(&timer);
  // gdsl.Plan(path);
  // time = timer_ns(&timer);
  // printf("replan path time= %ld\n", time);
  // printf("path: count=%d len=%f\n", path.count, path.cost);
  // fflush(stdout);


  // optimize path (experimental)

  /*
  timer_start(&timer);
  search.OptPath(path, optPath);
  time = timer_us(&timer);
  printf("opt path time= %ld us\n", time);
  printf("optPath: count=%lu len=%f\n", optPath.cells.size(), optPath.cost);
   */

  for (it = path.cells.begin(); it != path.cells.end(); ++it) {
    int x = pgrid->Index(it->c, 0);
    int y = pgrid->Index(it->c, 1);
    mapPath[y*width + x] = 2;// path.cells[i].p[1]*width +  path.cells[i].p[0]] = 2;
  }
  // save it to image for viewing
  save_map(mapPath, width, height, "path2.ppm");


  getchar();

  return 0;
}
