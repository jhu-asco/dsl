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
using CarTwistPath = dsl::GridPath<SE2Cell::PointType, SE2Cell::DataType, SE2Twist>;

vector<Vector3d> ToVector3dPath(const CarTwistPath &path, double gridcs) {
  vector<Vector3d> path3d;
  for (size_t i=0; i<path.connections.size(); i++){
    Vector3d v = path.connections[i];
    Vector3d axy = path.cells[i].c;
    Matrix3d g_from;
    se2_q2g(g_from, axy);
    double d = std::abs(v(1));
    int n_seg = ceil(d/ (2 * gridcs));
    double s = d/n_seg;
    for (int i_seg=0; i_seg<=n_seg; i_seg++){
      Matrix3d g, dg;
      se2_exp(dg, (s*i_seg / d) * v);
      g = g_from * dg;
      se2_g2q(axy, g);
      path3d.push_back(axy);
    }
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
  Params params(argv[1]);

  struct timeval timer;
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
  string mapName;
  if(!params.GetString("map", mapName)){
    cout<<"There is no string parameter named map. Exiting."<<endl;
    return -1;
  }
  // occupancy cell size
  Vector3d ocs;
  params.GetVector3d("ocs", ocs);

  // load an occupancy map from ppm file
  dsl::Map<bool, 2>::Ptr omap = LoadPpm(mapName, ocs.tail<2>());

  //**********************************************************************************/
  //************************************cmap******************************************/
  //**********************************************************************************/

  // get filename+ path for .cmap file that stores representing our occupancy map
  string cmapName;
  bool cmapValid = params.GetString("cmap", cmapName);

  // get geometry information
  Vector5d lboxoysb;
  bool use_geom = params.GetVector5d("geom", lboxoysb);
  CarGeom geom;
  if(use_geom)
    geom.set(lboxoysb);

  // dimensions are determined from occupancy map
  Vector3d xlb(-M_PI + ocs[0]/2, omap->xlb()[0], omap->xlb()[1]);
  Vector3d xub(M_PI + ocs[0]/2, omap->xub()[0], omap->xub()[1]);

  // configuration-space map
  shared_ptr<dsl::Map<bool, 3> > cmap;

  if(!cmapValid){
    cmap.reset(new dsl::Map<bool, 3>(xlb, xub, ocs));
    std::cout << "Making cmap... " << std::endl;
    timer_start(&timer);
    if (use_geom) {
      int nthreads; params.GetInt("nthreads", nthreads);
      cmap = MakeCmap(*omap, ocs(0), geom, nthreads);
    } else {
      cmap = MakeCmap(*omap, ocs(0));
    }
    long time = timer_us(&timer);
    printf("cmap construction time= %ld  us\n", time);

    cmapName = mapName;
    ReplaceExtension(cmapName, string("cmap"));
    cmap->Save(cmapName);
    std::cout << "Saved cmap " << cmapName << " with xlb=" << cmap->xlb().transpose() << " xub=" << cmap->xub().transpose() << " gs=" << cmap->gs().transpose() << std::endl;

  }else{
    cmap = dsl::Map<bool,3>::Load(cmapName);
    std::cout << "Loaded map with xlb=" << cmap->xlb().transpose() << " xub=" << cmap->xub().transpose() << " gs=" << cmap->gs().transpose() << std::endl;
  }

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
  shared_ptr<CarCost> cost;
  if( params.GetVector3d("wt", wt))
    cost.reset(new CarCost(grid, wt));
  else
    cost.reset(new CarCost(grid)); //use default wt

  //**********************************************************************************/
  //*************************************connectivity*********************************/
  //**********************************************************************************/

  // load car connectivity and set custom parameters
  CarPrimitiveConfig primcfg;
  CarTwistConnectivity connectivity(grid, *cost);
  int nl, na;
  params.GetBool("prim_fwdonly",primcfg.fwdonly);
  params.GetDouble("prim_tphioverlmax",primcfg.tphioverlmax);
  params.GetDouble("prim_lmin",primcfg.lmin);
  params.GetDouble("prim_lmax",primcfg.lmax);
  params.GetInt("prim_nl",nl); primcfg.nl = nl;
  params.GetDouble("prim_amax",primcfg.amax);
  params.GetInt("prim_na",na);primcfg.na = na;
  params.GetBool("prim_pert",primcfg.pert);
  connectivity.SetPrimitives(1,primcfg);

  //**********************************************************************************/
  //***************************************search graph*******************************/
  //**********************************************************************************/

  bool initExpand = false;
  params.GetBool("initExpand", initExpand);

  cout << "Creating a graph..." << endl;
  // create planner
  timer_start(&timer);
  GridSearch<Vector3d, Matrix3d, SE2Twist> search(grid, connectivity, *cost, initExpand);
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

  CarTwistPath path;
  if(start_set && goal_set){
    cout << "Created a graph with " << search.Vertices() << " vertices and " << search.Edges() << " edges." << endl;

    cout << "Planning a path..." << endl;
    // plan
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
    path3d = ToVector3dPath(path,grid.cs()[1]);
    cout << "Map and path saved to path.ppm" << endl;
  }else{
//    path3d.push_back(start);
//    path3d.push_back(goal);
    path3d.push_back(grid.CellCenter(start));
    path3d.push_back(grid.CellCenter(goal));
    cout << "Map, start and goal (no path available ) saved to path.ppm" << endl;
  }

  bool plot_car; //To plot the path as rectanges or just points
  params.GetBool("plot_car", plot_car);
  // save it to image for viewing
  if(plot_car){
    SavePpmWithPath(*omap, "path.ppm", 3, path3d, &geom);
  }else{
    SavePpmWithPath(*omap, "path.ppm", 3, path3d);
  }

  //Display the primitive at start
  vector<vector<Vector2d>> prims(0);
  bool gotprims = connectivity.GetPrims(start,prims);
  if(gotprims){
    SavePpmWithPrimitives(*omap, "prim.ppm",3, prims);
    cout << "Map, with primitives from start saved to prim.ppm" << endl;
  }

  return 0;
}
