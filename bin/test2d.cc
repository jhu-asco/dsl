#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include "gridsearch.h"
#include "grid2d.h"
#include "gridcost.h"
#include "grid2dconnectivity.h"

using namespace dsl;
using namespace std;
using namespace Eigen;


void save_map(const char* map, int width, int height, const char* filename)
{
  int i, ind;
  char data[width*height*3];
  FILE* file = fopen(filename, "w");
  assert(file);
  fprintf(file, "P6\n%d %d 255\n", width, height);
  ind = 0;
  for (i = 0; i < width*height; i++, ind+=3) {
    data[ind] = data[ind+1] = data[ind+2] = (char)(map[i]*100);
  }
  assert(ind == 3*width*height);
  assert((int)fwrite(data, sizeof(char), ind, file) == ind);
  fclose(file);
}


char* load_map(int* width, int* height, const char* filename)
{
  int i, size;
  char *map, *data;
  FILE* file = fopen(filename, "r");
  assert(file);
  assert(fscanf(file, "P6\n%d %d 255\n", width, height));
  size = (*width**height);
  map = (char*)malloc(size);
  data = (char*)malloc(size*3);
  assert(fread(data, sizeof(char), size*3, file));
  for (i = 0; i < size; i++)
    map[i] = (data[3*i] ? 1 : 0);
  free(data);
  fclose(file);
  return map;
}

void timer_start(struct timeval *time) 
{
  gettimeofday(time,(struct timezone*)0);
}

long timer_us(struct timeval *time) 
{
  struct timeval now;
  gettimeofday(&now,(struct timezone*)0);
  return 1000000*(now.tv_sec - time->tv_sec) + now.tv_usec - time->tv_usec;
}


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
  GridPath<2> path, optPath;
  int x, y;
  char mapPath[width*height];
  struct timeval timer;
  long time;

  // create a map
  double map[width*height];
  for (int i = 0; i < width*height; ++i)
    map[i] = 1000*(double)chmap[i];

  // this is just for display
  memcpy(mapPath, chmap, width*height);

  // create planner
  Grid2d grid(width, height, map, 1, 1, 1, 1e16);
  GridCost<2> cost;
  Grid2dConnectivity connectivity(grid);
  GridSearch<2> search(grid, connectivity, cost, true);

  search.SetStart(Vector2d(1, height/2));
  search.SetGoal(Vector2d(width - 2, height/2));

  // plan
  timer_start(&timer);
  search.Plan(path);
  time = timer_us(&timer);
  printf("plan path time= %ld\n", time);
  printf("path: count=%lu len=%f\n", path.cells.size(), path.len);
  // print results
  vector<Cell<2> >::iterator it;
  for (it = path.cells.begin(); it != path.cells.end(); ++it) {
    int id = grid.Id(it->c);
    mapPath[id] = 2;// path.cells[i].p[1]*width +  path.cells[i].p[0]] = 2;
  }
  printf("\n");
  for (y = 0; y < height; ++y) {
    for (x = 0; x < width; ++x) {
      printf("%d ", mapPath[y*width + x]);
    }
    printf("\n");
  }
  fflush(stdout);

  // save it to image for viewing
  save_map(mapPath, width, height, "path1.ppm");

  // follow path until (28,18)
  search.SetStart(Vector2d(28,18));

  // simulate closing the narrow passage
  if (1) {
    // by increasing the cost drastically
    search.SetCost(Vector2d(29,18), 1000);
    search.SetCost(Vector2d(30,18), 1000);
    search.SetCost(Vector2d(31,18), 1000);
  } else {
    // a better way: by simply removing the passage
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
  
  timer_start(&timer);
  search.OptPath(path, optPath);
  time = timer_us(&timer);
  printf("opt path time= %ld us\n", time);
  printf("optPath: count=%lu len=%f\n", optPath.cells.size(), optPath.len);
 

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
  for (y = 0; y < height; ++y) {
    for (x = 0; x < width; ++x) {
      printf("%d ", mapPath[y*width + x]);
    }
    printf("\n");
  }

  // save it to image for viewing
  save_map(mapPath, width, height, "path2.ppm");
 

  getchar();
  free(chmap);

  return 0;
}
