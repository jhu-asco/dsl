#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include "gridsearch.h"
#include "grid2d.h"
#include "gridcost.h"
#include "traversabilitycost.h"
#include "grid2dconnectivity.h"
#include "utils.h"

using namespace dsl;
using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{
  if (argc!=2) {
    cout << "Usage: $./test2d map.ppm" << endl;
    cout << "\t\t where map.ppm is a map graphics file" << endl;
    cout << "\t\t output will be written to graphics files path1.ppm and path2.pppm" << endl;
    return 0;
  }

  // load a map from ppm file
  assert(argc == 2);
  int width, height;
  char* chmap = load_map(width, height, argv[1]);
  char mapPath[width*height];
  double map[width*height];
  for (int i = 0; i < width*height; ++i)
    map[i] = 1000*(double)chmap[i];

  // this is just for display
  memcpy(mapPath, chmap, width*height);

  // create planner
  Grid2d grid(width, height, map, 1, 1, 1e16);
  TraversabilityCost<Vector2d, double> cost;
  // GridCost<Vector2d, double> cost;
  Grid2dConnectivity connectivity(grid);
  GridSearch<Vector2d, double, Vector2d> search(grid, connectivity, cost);
  GridPath<Vector2d, double, Vector2d> path, optPath;

  search.SetStart(Vector2d(1, height/2));
  search.AddGoal(Vector2d(width - 2, height/2));

  std::cout  << "Added start and goal" << std::endl;

  // plan
  struct timeval timer;
  timer_start(&timer);
  search.Plan(path);
  long time = timer_us(&timer);
  printf("plan path time= %ld\n", time);
  printf("path: count=%lu len=%f\n", path.cells.size(), path.cost);

  // print results
  vector<Cell<Vector2d, double> >::iterator it;
  for (it = path.cells.begin(); it != path.cells.end(); ++it) {
    int id = grid.Id(it->c);
    mapPath[id] = 2;
  }
  printf("\n");
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      printf("%d ", mapPath[y*width + x]);
    }
    printf("\n");
  }
  fflush(stdout);

  // save it to image for viewing
  save_map(mapPath, width, height, "path1.ppm");

  cout << "Map and path saved to path1.ppm... Press Enter to simulated replanning after closing a passage." << endl;
  getchar();

  // follow path until (28,18)
  search.SetStart(Vector2d(28,18));

  // simulate closing the narrow passage
  if (1) {
    // by increasing the cost drastically
    search.SetCost(Vector2d(29,18), 1000);
    search.SetCost(Vector2d(30,18), 1000);
    search.SetCost(Vector2d(31,18), 1000);
  } else {
    // another way: by simply removing the passage
    search.RemoveCell(Vector2d(30,18));
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

  timer_start(&timer);
  connectivity.OptPath(path, optPath);
  time = timer_us(&timer);
  printf("opt path time= %ld us\n", time);
  printf("optPath: count=%lu len=%f\n", optPath.cells.size(), optPath.cost);


  for (it = path.cells.begin(); it != path.cells.end(); ++it) {
    int id = grid.Id(it->c);
    mapPath[id] = 2;
  }
  printf("\n");

  for (it = optPath.cells.begin(); it != optPath.cells.end(); ++it) {
    int id = grid.Id(it->c);
    mapPath[id] = 3;
  }

  printf("\n");
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      printf("%d ", mapPath[y*width + x]);
    }
    printf("\n");
  }

  // save it to image for viewing
  save_map(mapPath, width, height, "path2.ppm");
  cout << "Map and path saved to path2.ppm... " << endl;
  free(chmap);

  return 0;
}
