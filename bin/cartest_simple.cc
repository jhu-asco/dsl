#include <carconnectivity.h>
#include <string.h>
#include "gridsearch.h"
#include "cargrid.h"
#include "carcost.h"
#include "utils.h"
#include "gridutils.h"
#include "params.h"
#include <fstream>

using namespace dsl;
using namespace std;
using namespace Eigen;

using CarPath = GridPath<SE2Cell::PointType, SE2Cell::DataType, SE2Path>;

/**
 * Converts CarPath to a vector of axy points representing the path. axy: yaw, x and y
 * @param path
 * @return vector of axy points
 */
vector<Vector3d> ToVector3dPath(const CarPath &path) {
  vector<Vector3d> path3d;
  for (auto&& connection : path.connections)
    for (auto&& g : connection){
      Vector3d axy; axy.tail<2>() = g.topRightCorner<2,1>();
      axy(0) = atan2(-g(0,1),g(0,0));
      path3d.push_back(axy);
    }
  return path3d;
}

int main(int argc, char** argv)
{
  if (argc!=2) {
    cout << "Usage: $./cartest params.cfg" << endl;
    cout << "\t\t where params.cfg is a parameter file" << endl;
    cout << "\t\t output will be written to graphics files path1.ppm and path2.ppm" << endl;
    return 0;
  }
  assert(argc == 2);

  //Initialize the parse for .cfg file. The .cfg file defines the planning problem.
  Params params(argv[1]);

  // The problem setup for planning the path of a car on a plane is done as follows:
  // 1. Load 2D(x and y) occupancy map, or omap, from image file(.ppm).
  // 2. Create or load configuration space( yaw, x and y) occupancy map, or cmap.
  //     a.Loading cmap is done from a .cmap file, if already available.
  //     b.Creating cmap is done by taking an omap and doing a Minkowski sum with the shape of the car.
  //        CarGeom object deals with the details of shape of the car. If car geometry is not available
  //        then the omap is simply stack(for all angles) to create cmap.
  // 3. Create a lower resolution(than cmap) grid, or lattice, that for every unoccupied cell. This grid
  //     is what the search algorithm will work on.
  // 4. Create cost object whose job is to provide, given a lattice, the cost of motion from one cell of
  //     the lattice to another. It also provides heuristics cost for the same. Calculation of heuristics
  //     cost is faster but the cost is an underestimator of the real cost.
  // 5. Create connectivity object whose job is to provide the search algorithm with motion primitives
  //     (along with cost), for a given cell, to its neighboring cells. In case of car these motion primitives
  //     are se2 twists. It uses the cost object above to get the primitive costs.
  // 6. Create the main search graph which uses the lattice, cost, connectivity create above.
  // 7. Plan the path from a start and goal position using the search graph created above
  // 8. Plot the results on the occupancy map we read before.

  //**********************************************************************************/
  //************************************omap******************************************/
  //**********************************************************************************/

  // get filename+ path for the image file(.ppm) representing our occupancy map
  string map_img_filename;
  if(!params.GetString("map", map_img_filename)){
    cout<<"There is no string parameter named map. Exiting."<<endl;
    return -1;
  }
  // occupancy cell size
  Vector3d ocs;
  params.GetVector3d("ocs", ocs);

  // load an occupancy map from ppm file
  dsl::Map<bool, 2>::Ptr omap = LoadPpm(map_img_filename, ocs.tail<2>());

  //**********************************************************************************/
  //************************************cmap******************************************/
  //**********************************************************************************/

  // get filename+ path for .cmap file that stores representing our occupancy map
  string cmap_filename;
  bool load_cmap = params.GetString("cmap", cmap_filename);

  // get geometry information
  Vector5d l_b_ox_oy_sb;
  bool use_geom = params.GetVector5d("geom", l_b_ox_oy_sb);
  CarGeom geom;
  if(use_geom)
    geom.set(l_b_ox_oy_sb);

  // dimensions are determined from occupancy map
  Vector3d xlb(-M_PI + ocs[0]/2, omap->xlb()[0], omap->xlb()[1]);
  Vector3d xub(M_PI + ocs[0]/2, omap->xub()[0], omap->xub()[1]);

  // configuration-space map
  dsl::Map<bool, 3>::Ptr cmap;

  struct timeval timer;
  if(!load_cmap){ //load cmap
    cmap.reset(new dsl::Map<bool, 3>(xlb, xub, ocs));
    std::cout << "Making cmap... " << std::endl;
    timer_start(&timer);
    if (use_geom) {
      int nthreads;
      params.GetInt("nthreads", nthreads);
      cmap = MakeCmap(*omap, ocs(0), geom, nthreads);
    } else {
      cmap = MakeCmap(*omap, ocs(0));
    }
    long time = timer_us(&timer);
    printf("cmap construction time= %ld  us\n", time);

    cmap_filename = ReplaceExtension(map_img_filename, string("cmap"));
    cmap->Save(cmap_filename);
    std::cout << "Saved cmap " << cmap_filename << " with xlb=" << cmap->xlb().transpose() << " xub=" << cmap->xub().transpose() << " gs=" << cmap->gs().transpose() << std::endl;

  }else{ //create cmap
    cmap = dsl::Map<bool,3>::Load(cmap_filename);
    std::cout << "Loaded map with xlb=" << cmap->xlb().transpose() << " xub=" << cmap->xub().transpose() << " gs=" << cmap->gs().transpose() << std::endl;
  }

  //save all slices of cmap in "cmap_slices" folder
  SavePpm(*cmap, "cmap_slices");

  //**********************************************************************************/
  //*************************************lattice**************************************/
  //**********************************************************************************/

  Vector3d gcs; // grid cell size (normally larger than ocs)
  params.GetVector3d("gcs", gcs);
  CarGrid grid(*cmap, gcs);

  //**********************************************************************************/
  //*************************************cost**************************************/
  //**********************************************************************************/

  Vector3d wt; //weights for twist element
  CarCost::Ptr cost;
  if( params.GetVector3d("wt", wt))
    cost.reset(new CarCost(grid, wt));
  else
    cost.reset(new CarCost(grid)); //use default wt

  //**********************************************************************************/
  //*************************************connectivity*********************************/
  //**********************************************************************************/

  // load car connectivity and set custom parameters
  double dt = .25;
  double vx = 4;
  double kmax = 0.57;
  int kseg = 4;
  bool onlyfwd = false;
  bool allow_slip = true;
  params.GetDouble("dt", dt);
  params.GetDouble("vx", vx);
  params.GetDouble("kmax", kmax);
  params.GetInt("kseg", kseg);
  params.GetBool("onlyfwd", onlyfwd);
  params.GetBool("allow_slip", allow_slip);
  CarPathConnectivity connectivity(grid,*cost,allow_slip);
  connectivity.SetPrimitives(dt, vx, kmax, kseg, onlyfwd);

  //**********************************************************************************/
  //***************************************search graph*******************************/
  //**********************************************************************************/

  bool initExpand = false;
  params.GetBool("initExpand", initExpand);

  cout << "Creating a graph..." << endl;
  // create planner
  timer_start(&timer);
  GridSearch<Vector3d, Matrix3d, SE2Path> search(grid, connectivity, *cost, initExpand);
  CarPath path;
  long time = timer_us(&timer);
  printf("graph construction time= %ld  us\n", time);

  //**********************************************************************************/
  //***************************************plan***************************************/
  //**********************************************************************************/

  // get start and goal
  Vector3d goal, start;
  params.GetVector3d("start", start);
  params.GetVector3d("goal", goal);

  bool start_set = search.SetStart(start);
  bool goal_set = search.SetGoal(goal);

  if(start_set && goal_set){
    cout << "Created a graph with " << search.Vertices() << " vertices and " << search.Edges() << " edges." << endl;
    // plan
    cout << "Planning a path..." << endl;
    timer_start(&timer);
    search.Plan(path);
    time = timer_us(&timer);
    printf("plan path time= %ld  us\n", time);
    printf("path: edges=%lu len=%f\n", path.cells.size(), path.cost);

    cout << "Graph has " << search.Vertices() << " vertices and " << search.Edges() << " edges. " << endl;
  }

  //**********************************************************************************/
  //***************************************plot***************************************/
  //**********************************************************************************/

  vector<Vector3d> path3d;
  if(start_set && goal_set){
    path3d = ToVector3dPath(path);
    cout << "Map and path saved to path.ppm" << endl;
  }else{
    path3d.push_back(grid.CellCenter(start));
    path3d.push_back(grid.CellCenter(goal));
    cout << "Map, start and goal (no path available) saved to path.ppm" << endl;
  }

  bool plot_car; //To plot the path as rectanges or just points
  params.GetBool("plot_car", plot_car);
  if(plot_car)
    SavePpmWithPath(*omap, "path.ppm", path3d, 3, &geom);
  else
    SavePpmWithPath(*omap, "path.ppm", path3d, 3);


  //Display the primitive at start
  vector<vector<Vector2d>> prims(0);
  bool gotprims = connectivity.GetPrims(start,prims);
  if(gotprims){
    SavePpmWithPrimitives(*omap, "prim.ppm", prims, 3);
    cout << "Map, with primitives from start saved to prim.ppm" << endl;
  }

  return 0;
}
