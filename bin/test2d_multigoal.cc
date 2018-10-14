#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include "gridsearch.h"
//#include "gridlpastar.h"

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
  GridSearch<Vector2d, double, Vector2d> search(grid, connectivity, cost, Method::kDstar);
  // to try LpAStar forward search uncomment below and comment out the preceding line
  // GridSearch<Vector2d, double, Vector2d> search(grid, connectivity, cost, Method::kLpAstar);
  GridPath<Vector2d, double, Vector2d> path, opt_path;

  search.setStart(Vector2d(1, height/2));

  // create a goal region with size goalWidth x goalHeight on the far right middle
  double goalWidth = 10;
  double goalHeight = 10;
  double x_start = std::max(0.0,(double)width - goalWidth);
//  double y_start = std::max(height/2.0 - goalHeight/2.0, 0.0);
  double y_start = 1;
  for (double x = x_start; x < std::min(x_start + goalWidth, (double)width); x += 1)
    for (double y = y_start; y < std::min(y_start + goalHeight, (double)height); y += 1)
      search.addGoal(Vector2d(x, y));

  // plan
  struct timeval timer;
  timer_start(&timer);
  search.plan(path);
  long time = timer_us(&timer);
  printf("plan path time= %ld\n", time);
  printf("path: count=%lu len=%f\n", path.cells.size(), path.cost);

  // print results
  vector<Cell<Vector2d, double> >::iterator it;
  for (it = path.cells.begin(); it != path.cells.end(); ++it) {
    int id = grid.computeId(it->centr);
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

  return 0;
}
