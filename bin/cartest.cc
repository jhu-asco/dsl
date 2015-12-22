#include <string.h>
#include "gridsearch.h"
#include "cargrid.h"
#include "carcost.h"
#include "carconnectivity.h"
#include "utils.h"

using namespace dsl;
using namespace std;
using namespace Eigen;


int main(int argc, char** argv)
{
  if (argc!=2) {
    cout << "Usage: $./test map.ppm" << endl;
    cout << "\t\t where map.ppm is a map graphics file" << endl;
    cout << "\t\t output will be written to graphics files path1.ppm and path2.pppm" << endl;
    return 0;
  }
  assert(argc == 2);
  int width, height; 
  char* chmap = load_map(&width, &height, argv[1]);
  GridPath<3> path, optPath;
  char mapPath[width*height];
  struct timeval timer;
  long time;

  // create a map
  double map[width*height];
  for (int i = 0; i < width*height; ++i)
    map[i] = 1000*(double)chmap[i];

  // this is just for display
  memcpy(mapPath, chmap, width*height);
  cout << "Creating a graph..." << endl;

  // create planner
  timer_start(&timer);

  CarGrid grid(width, height, map, .1, .1, M_PI/16, 1, 0.5);
  CarCost cost;
  CarConnectivity connectivity(grid);
  connectivity.dt = 1; // length in seconds of primitives
  GridSearch<3> search(grid, connectivity, cost, false);

  time = timer_us(&timer);
  printf("graph construction time= %ld  us\n", time);

  search.SetStart(Vector3d(.1, grid.xub[1]/2, 0));
  search.SetGoal(Vector3d(grid.xub[0] - .1, grid.xub[1]/2, 0));
  //  search.SetGoal(Vector3d(grid.xub[0]*.5, grid.xub[1]*.58, 15.0/16*M_PI));

  cout << "Created a graph with " << search.Vertices() << " vertices and " << search.Edges() << " edges. " << endl;

  cout << "Planning a path..." << endl;
  // plan
  timer_start(&timer);
  search.Plan(path);
  time = timer_us(&timer);
  printf("plan path time= %ld  us\n", time);
  printf("path: edges=%lu len=%f\n", path.cells.size(), path.len);

  cout << "Graph has " << search.Vertices() << " vertices and " << search.Edges() << " edges. " << endl;


  // print results
  vector<Cell<3> >::iterator it; 
  for (it = path.cells.begin(); it != path.cells.end(); ++it) {
    //    printf("(%d,%d) ",  path.cells[i].p[0], path.cells[i].p[1]);    
    int x = grid.Index(it->c, 0);
    int y = grid.Index(it->c, 1);
    mapPath[y*width + x] = 2;// path.cells[i].p[1]*width +  path.cells[i].p[0]] = 2;
  }
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
  save_map(mapPath, width, height, "path1.ppm");
  cout << "Map and path saved to path1.ppm" << endl;
  
  return 1;

  // follow path until middle
  Vector3d c = path.cells[path.cells.size()/2].c;
  search.SetStart(c);

  // simulate closing the narrow passage
  if (0) {
    for (double a = grid.xlb[2]+grid.cs[2]/2; a < grid.xub[2] + grid.cs[2]/2; a += grid.cs[2]) {
      // by increasing the cost drastically
      //      search.SetCost(Vector3d(29,18,a), 1000);
      //      search.SetCost(Vector3d(30,18,a), 1000);
      //      search.SetCost(Vector3d(31,18,a), 1000);
    }
  } else {
  // a better way: by simply removing the passage
    for (double y = c[1] - 2; y < c[1] +2; ++y)  {
      for (double a = grid.xlb[2]+grid.cs[2]/2; a < grid.xub[2] + grid.cs[2]/2; a += grid.cs[2]) {     
        search.RemoveCell(Vector3d(c[0],y,a));
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
  printf("path: count=%lu len=%f\n", path.cells.size(), path.len);
  fflush(stdout);
  

  // bypass the old vertex
  //  gdsl.AddEdge(29,18,31,18);  
  // // replan
 // timer_start(&timer);
 // gdsl.Plan(path);
 // time = timer_ns(&timer);
 // printf("replan path time= %ld\n", time);
 // printf("path: count=%d len=%f\n", path.count, path.len);
 // fflush(stdout);
  
  
  // optimize path (experimental)
  
  /*
  timer_start(&timer);
  search.OptPath(path, optPath);
  time = timer_us(&timer);
  printf("opt path time= %ld us\n", time);
  printf("optPath: count=%lu len=%f\n", optPath.cells.size(), optPath.len);
  */

  for (it = path.cells.begin(); it != path.cells.end(); ++it) {
    int x = grid.Index(it->c, 0);
    int y = grid.Index(it->c, 1);
    mapPath[y*width + x] = 2;// path.cells[i].p[1]*width +  path.cells[i].p[0]] = 2;
  }
  // save it to image for viewing
  save_map(mapPath, width, height, "path2.ppm");
  

  getchar();

  return 0;
}
