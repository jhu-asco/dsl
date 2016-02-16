#include <string.h>
#include "gridsearch.h"
#include "cargrid.h"
#include "carcost.h"
#include "carconnectivity.h"
#include "utils.h"
#include "params.h"
#include <fstream>

using namespace dsl;
using namespace std;
using namespace Eigen;

vector<Vector2d> ToVector2dPath(const SE2Path &path) {
  vector<Vector2d> path2d;
  for (auto&& cell : path.cells)
    path2d.push_back(Vector2d(cell.c.tail<2>()));
  return path2d;
}


int main(int argc, char** argv)
{
  if (argc!=2) {
    cout << "Usage: $./cartest params.cfg" << endl;
    cout << "\t\t where params.cfg is a parameter file" << endl;
    cout << "\t\t output will be written to graphics files path1.ppm and path2.pppm" << endl;
    return 0;
  }
  assert(argc == 2);

  Params params(argv[1]);

  // get start and goal
  Vector3d goal, start;
  params.GetVector3d("start", start);
  params.GetVector3d("goal", goal);

  // get map
  string mapName;
  params.GetString("map", mapName);

  // occupancy cell size
  Vector3d ocs;
  params.GetVector3d("ocs", ocs);  

  // grid cell size (normally larger than ocs)
  Vector3d gcs;
  params.GetVector3d("gcs", gcs);  

  // load an occupancy map from ppm file
  dsl::Map<bool, 2> omap = load(mapName.c_str(), ocs.tail<2>());

  // a map that we'll use for display
  dsl::Map<bool, 2> dmap = omap;


  // dimensions are determined from occupancy map
  Vector3d xlb(-M_PI + gcs[0]/2, omap.xlb[0], omap.xlb[1]);
  Vector3d xub(M_PI + gcs[0]/2, omap.xub[0], omap.xub[1]);

  // configuration-space map
  dsl::Map<bool, 3> cmap(xlb, xub, ocs);
  
  CarGrid::MakeMap(omap, cmap);
  
  CarGrid grid(cmap, gcs);
  CarCost cost;
  CarConnectivity connectivity(grid);

  double vx, w, dt;
  params.GetDouble("vx", vx);
  params.GetDouble("w", w);
  params.GetDouble("dt", dt);
  connectivity.SetPrimitives(vx, w, dt);

  cout << "Creating a graph..." << endl;
  // create planner
  struct timeval timer;
  timer_start(&timer);
  
  GridSearch<Vector3d, Matrix3d> search(grid, connectivity, cost, false);
  SE2Path path;

  long time = timer_us(&timer);
  printf("graph construction time= %ld  us\n", time);

  search.SetStart(start);
  search.SetGoal(goal);

  cout << "Created a graph with " << search.Vertices() << " vertices and " << search.Edges() << " edges." << endl;

  cout << "Planning a path..." << endl;
  // plan
  timer_start(&timer);
  search.Plan(path);
  time = timer_us(&timer);
  printf("plan path time= %ld  us\n", time);
  printf("path: edges=%lu len=%f\n", path.cells.size(), path.cost);

  cout << "Graph has " << search.Vertices() << " vertices and " << search.Edges() << " edges. " << endl;
  /*
  // draw lattice
  for (int i =0; i < grid.nc; ++i) {
    if (grid.cells[i]) {
      int id = dmap.Id(grid.cells[i]->c.tail<2>());
      dmap.cells[id] = 1;// path.cells[i].p[1]*width +  path.cells[i].p[0]] = 2;
    }
  }
  */

  /*
  printf("\n");
  for (y = 0; y < height; ++y) {
    for (x = 0; x < width; ++x) {
      printf("%d ", mapPath[y*width + x]);
    }
    printf("\n");
  }
  fflush(stdout);
  */
  
  // save it to image for viewing
  vector<Vector2d> path2d = ToVector2dPath(path);
  save(dmap, "path1.ppm", &path2d);

  cout << "Map and path saved to path1.ppm" << endl;

  return 0;
  /*
  // follow path until middle
  Vector3d c = path.cells[path.cells.size()/2].c;
  search.SetStart(c);

  // simulate closing the narrow passage
  if (0) {
    for (double a = grid.xlb[0]+grid.cs[0]/2; a < grid.xub[0] + grid.cs[0]/2; a += grid.cs[0]) {
      // by increasing the cost drastically
      //      search.SetCost(Vector3d(29,18,a), 1000);
      //      search.SetCost(Vector3d(30,18,a), 1000);
      //      search.SetCost(Vector3d(31,18,a), 1000);
    }
  } else {
  // a better way: by simply removing the passage
    for (double y = c[1] - 2; y < c[1] +2; ++y)  {
      for (double a = grid.xlb[0]+grid.cs[0]/2; a < grid.xub[0] + grid.cs[0]/2; a += grid.cs[0]) {     
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
  
  */
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

  /*
  for (it = path.cells.begin(); it != path.cells.end(); ++it) {
    int x = grid.Index(it->c, 0);
    int y = grid.Index(it->c, 1);
    mapPath[y*width + x] = 2;// path.cells[i].p[1]*width +  path.cells[i].p[0]] = 2;
  }
  // save it to image for viewing
  save_map(mapPath, width, height, "path2.ppm");
  */

  return 0;
}
