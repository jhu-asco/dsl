#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include "gridsearch.h"
#include "grid3d.h"
#include "gridcost.h"
#include "grid3dconnectivity.h"


#define STDOUT_DEBUG

using namespace dsl;
using namespace std;

void save_ply(double* map, int length, int width, int height, const char* filename)
{
  int i,j,k;
  int numPts = 0;
  FILE* file = fopen(filename, "w");
  
  for(k = 0; k < height; k++){
    for(j = 0; j < width; j++){
      for(i = 0; i < length; i++){
        if(map[k*length*width + j*length + i] > 0)
          numPts++;
      }
    }
  }

  assert(file);
  fprintf(file,"ply\nformat ascii 1.0\nelement vertex %d\nproperty float x\nproperty float y\nproperty float z\nend_header\n", numPts);
  for(k = 0; k < height; k++){
    for(j = 0; j < width; j++){
      for(i = 0; i < length; i++){
        if(map[k*length*width + j*length + i] > 0)
          fprintf(file, "%d %d %d\n", i, j, k);
      }
    }
  }
  fclose(file);
}


void save_map_with_path(unsigned const char* map, unsigned int width, unsigned int height,
  std::vector<Vector3d>& path, const char* filename)
{
  unsigned int i, ind;
  unsigned char data[width*height*3];
  FILE* file = fopen(filename, "w");
  assert(file);
  fprintf(file, "P6\n%d %d 255\n", width, height);
  ind = 0;
  for (i = 0; i < width*height; i++, ind+=3) {
    data[ind] = data[ind+1] = data[ind+2] = (map[i]);
  }
  for(i = 0; i < path.size(); i++){
    int pidx = int(floor(path[i](1))*width + floor(path[i](0)));
    data[3*pidx] = path[i](2)*255/50;
    data[3*pidx+1] = 0;
    data[3*pidx+2] = 0;
  }
  assert(ind == 3*width*height);
  assert(fwrite(data, sizeof(unsigned char), ind, file) == ind);
  fclose(file);
}

void save_map_with_cells(unsigned const char* map, unsigned int width, unsigned int height, 
  GridPath<3>& gp, const char* filename)
{
  std::vector<Vector3d> path;
  for(unsigned int i = 0; i < gp.cells.size(); i++)
  {
    path.push_back(gp.cells[i].c);
  }
  save_map_with_path(map, width, height, path, filename);
 
}

void save_map(const char* map, int width, int height, const char* filename)
{
  int i, ind;
  char data[width*height*3];
  FILE* file = fopen(filename, "w");
  assert(file);
  fprintf(file, "P6\n%d %d 255\n", width, height);
  ind = 0;
  for (i = 0; i < width*height; i++, ind+=3) {
    data[ind] = data[ind+1] = data[ind+2] = (char)(map[i]);
  }
  assert(ind == 3*width*height);
  assert((int)fwrite(data, sizeof(char), ind, file) == ind);
  fclose(file);
}


unsigned char* load_map(int* width, int* height, const char* filename)
{
  int i, size;
  unsigned char *map, *data;
  FILE* file = fopen(filename, "r");
  assert(file);
  assert(fscanf(file, "P6\n%d %d 255\n", width, height));
  size = (*width**height);
  map = (unsigned char*)malloc(size);
  data = (unsigned char*)malloc(size*3);
  assert(fread(data, sizeof(unsigned char), size*3, file));
  for (i = 0; i < size; i++)
    map[i] = (data[3*i] ? 1 : 0);
  free(data);
  fclose(file);
  return map;
}

unsigned char* load_heightmap(unsigned int* length, unsigned int* width, const char* filename)
{
  int i, size;
  unsigned char *map, *data;
  FILE* file = fopen(filename, "r");
  assert(file);
  assert(fscanf(file, "P6\n%d %d 255\n", length, width));
  size = (*length* (*width));
  map = (unsigned char*)malloc(size);
  data = (unsigned char*)malloc(size*3);
  assert(fread(data, sizeof(unsigned char), size*3, file));
  for (i = 0; i < size; i++)
  {
    map[i] = data[3*i];
  }
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
  if (argc<2) {
    cout << "Usage: $./test3D map.ppm" << endl;
    cout << "\t\t where map.ppm is a map graphics file" << endl;
    cout << "\t\t output will be written to graphics files path1.ppm" << endl;
    return 0;
  }
  assert(argc >= 2);
  unsigned int length, width, height = 50; 
  unsigned char* chmap; // = load_map(&width, &height, argv[1]);
  //if(argc == 3 && string(argv[2]) == "--height")
  chmap = load_heightmap(&length, &width, argv[1]);
  
  cout << "Height map loaded" << endl;
  cout << "length=" << length << endl;
  cout << "width=" << width << endl;

  GridPath<3> path, splineCells;
  std::vector<Vector3d> splinePath;
  unsigned char mapPath[length*width];
  struct timeval timer;
  long time;
  unsigned int i,j;

  // create a map
  double *map = (double*)malloc(length*width*height*sizeof(double));
  if(!map)
  {
    cout << "Failed to malloc map" << endl;
    return 0;
  }

  for(j = 0; j < height; ++j){
    for (i = 0; i < length*width; ++i){ 
      //cout << j*length*width + i << " " << (height > int(chmap[i]*50/255)) << " " << height << " " << int(chmap[i]*50/255) << endl;
      map[j*length*width + i] = (j > int(chmap[i])*50/255 ? 0 : 1000);
    }
  }
  save_ply(map, length, width, height, "map.ply"); 
  cout << "Internal map created" << endl;

  // this is just for display
  memcpy(mapPath, chmap, length*width);
  cout << "Copied display map" << endl;
  
  // create planner
  Grid3d grid(length, width, height, map, 1, 1, 1, 1, 1e16);
  GridCost<3> cost;
  Grid3dConnectivity connectivity(grid);
  GridSearch<3> gdsl(grid, connectivity, cost, true);

  cout << "Initialized planner" << endl;
  int startx = length/4.;
  int starty = width/2.;
  int startz = 15;
  int goalx = 5*length/8.;
  int goaly = width/2.+2;//width/4 - 10;
  int goalz = 15;
  gdsl.SetStart(Vector3d(startx, starty, startz));
  gdsl.SetGoal(Vector3d(goalx, goaly, goalz));
  printf("start (%d,%d,%d) \n", startx, starty, startz);
  printf("goal (%d,%d,%d) \n", goalx, goaly, goalz);
  

  cout << "Setup planner" << endl;

  // plan
  timer_start(&timer);
  gdsl.Plan(path);
  time = timer_us(&timer);
  printf("plan path time= %ld\n", time);
  printf("path: count=%lu len=%f\n", path.cells.size(), path.len);

  // print results
  for (i = 0; i < path.cells.size(); ++i) {
    printf("(%f,%f,%f) \n", path.cells[i].c[0], path.cells[i].c[1], path.cells[i].c[2]);
  }

  //gdsl.SplinePath(path, splinePath, /*splineCells,*/ .1);
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
  cout << "Saving path..." << endl;
  save_map_with_cells(mapPath, length, width, path, "path1.ppm");
  //cout << "Saving spline path..." << endl;
  //save_map_with_path(mapPath, length, width, splinePath, "spline_path1.ppm");

  /*
  // follow path until (28,18)
  gdsl.SetStart(28,18);

  // simulate closing the narrow passage
  if (0) {
    // one way
    gdsl.SetCost(29,18,1000);
    gdsl.SetCost(30,18,1000);
    gdsl.SetCost(31,18,1000);
  } else {
    // a better way
    //    gdsl.RemoveVertex(29,18);
    gdsl.RemoveVertex(30,18);
    //    gdsl.RemoveVertex(31,18);
  }
  */
  // this is just for display
  /*
  memcpy(mapPath, chmap, width*height);
  gdsl.OptPath(path,optPath);

  for (i = 0; i < path.count; ++i) {
    //printf("(%d,%d) ", path.pos[2*i], path.pos[2*i+1]);
    mapPath[optPath.pos[3*i+1]*length + optPath.pos[3*i]] = 2;
  }
  save_map(mapPath, length, width, "optpath1.ppm");
  */
  //mapPath[18*width + 29] = 1;  mapPath[18*width + 30] = 1;  mapPath[18*width + 31] = 1;
  /*
  // replan
  timer_start(&timer);
  gdsl.Plan(path);
  time = timer_us(&timer);
  printf("replan path time= %ld us\n", time);
  printf("path: count=%d len=%f\n", path.count, path.len);
  fflush(stdout);
  */

  /*
  // bypass the old vertex
  gdsl.AddEdge(29,18,31,18);
  
  // replan
  timer_start(&timer);
  gdsl.Plan(path);
  time = timer_ns(&timer);
  printf("replan path time= %ld\n", time);
  printf("path: count=%d len=%f\n", path.count, path.len);
  fflush(stdout);
  */
  /*
  // optimize path (experimental)
  timer_start(&timer);
  gdsl.OptPath(path, optPath);
  time = timer_us(&timer);
  printf("opt path time= %ld us\n", time);
  printf("optPath: count=%d len=%f\n", optPath.count, optPath.len);
  
  for (i = 0; i < path.count; ++i) {
    printf("(%d,%d) ", path.pos[2*i], path.pos[2*i+1]);
    mapPath[path.pos[2*i+1]*width + path.pos[2*i]] = 2;
  }
  printf("\n");
  for (i = 0; i < optPath.count; ++i) {
    printf("(%d,%d) ", optPath.pos[2*i], optPath.pos[2*i+1]);
    fflush(stdout);
    mapPath[optPath.pos[2*i+1]*width + optPath.pos[2*i]] = 3;
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
  */

  //getchar();

  return 0;
}
