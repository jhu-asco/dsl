#include <carconnectivity.h>
#include <string.h>
#include "gridsearch.h"
#include "terrainse2grid.h"
#include "terrainse2gridcost.h"
#include "utils.h"
#include "gridutils.h"
#include "params.h"
#include <fstream>

using namespace dsl;
using namespace std;
using namespace Eigen;

using CarTwistPath = dsl::GridPath<TerrainCell::PointType, TerrainCell::DataType, SE2Twist>;

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
    cout << "Usage: $./terraincartest params.cfg" << endl;
    cout << "\t\t where params.cfg is a parameter file containing planner "
        <<"configurations, terrain maps and car geometry" << endl;
    cout << "\t\t output will be written to graphics files path1.ppm and path2.ppm" << endl;
    return 0;
  }
  assert(argc == 2);
  Params params(argv[1]);

  struct timeval timer;
  IOFormat eigformat(StreamPrecision, DontAlignCols, ", ", ", ", "", "", "(", ")");

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
  cout<<"\nLoading terrain map from .tmap file"<<endl;
  // get filename+ path for the terrain file(.tmap) containing the height and
  //  traversibility of terrain
  string tmapfile; ;
  if(!params.GetString("tmap", tmapfile)){
    cout<<"There is no string parameter named tmap. Exiting."<<endl;
    return -1;
  }

  dsl::Map<TerrainData,2>::Ptr tmap =LoadTmap(tmapfile);
  if(!tmap){
    cout<<"Unable to load the .tmap file:"<<tmapfile<<endl;
    return -1;
  }
  cout<<"  Loaded terrain map with"<<endl;
  cout<<"    xlb: "<<tmap->xlb().transpose().format(eigformat)<<endl;
  cout<<"    xub: "<<tmap->xub().transpose().format(eigformat)<<endl;
  cout<<"    gs: "<<tmap->gs().transpose().format(eigformat)<<endl;

  //**********************************************************************************/
  //************************************cmap******************************************/
  //**********************************************************************************/
  cout<<"\nconfiguration space occupancy map(cmap)"<<endl;
  // get filename+ path for .cmap file that stores representing our occupancy map
  string cmapfile;
  bool cmapValid = params.GetString("cmap", cmapfile);

  // get geometry information
  Vector5d lboxoysb;
  bool use_geom = params.GetVector5d("geom", lboxoysb);
  CarGeom geom;
  if(use_geom)
    geom.set(lboxoysb);

  // configuration-space map
  shared_ptr<dsl::Map<bool, 3> > cmap;

  if(params.GetString("cmap", cmapfile)){ //load cmap from file
    cmap = dsl::Map<bool,3>::Load(cmapfile);
    if(!cmap)
      return -1;
    cout << "  Loaded cmap from "<<cmapfile<<" with"<<endl;
    cout << "    xlb=" << cmap->xlb().transpose().format(eigformat)<<endl;
    cout << "    xub=" << cmap->xub().transpose().format(eigformat)<<endl;
    cout << "    gs="  << cmap->gs().transpose().format(eigformat)<<endl;
  }else { // make cmap
    double ocsy; params.GetDouble("ocsy", ocsy);
    timer_start(&timer);
    if (use_geom) {
      std::cout << "  Making cmap with car geometry " << std::endl;
      int nthreads; params.GetInt("nthreads", nthreads);
      cmap = MakeCmap(*tmap, ocsy, geom, nthreads);
    } else {
      std::cout << "  Making cmap with point geometry " << std::endl;
      cmap = MakeCmap(*tmap, ocsy);
    }
    if(!cmap)
      return -1;
    long time = timer_us(&timer);
    printf("  cmap construction time= %ld  us\n", time);

    cmapfile = tmapfile;
    ReplaceExtension(cmapfile, string("cmap"));
    cmap->Save(cmapfile);
    cout << "  Saved cmap in "<<cmapfile <<" with"<<endl;
    cout << "    xlb=" << cmap->xlb().transpose().format(eigformat)<<endl;
    cout << "    xub=" << cmap->xub().transpose().format(eigformat)<<endl;
    cout << "    gs="  << cmap->gs().transpose().format(eigformat)<<endl;
  }
  SavePpm(*cmap,"cmap_slices");

  //**********************************************************************************/
  //*************************************lattice**************************************/
  //**********************************************************************************/
  cout<<"\nMain lattice"<<endl;

  Vector3d gcs;  // grid cell size (normally larger than ocs)
  params.GetVector3d("gcs", gcs);

  TerrainSE2Grid grid(*cmap, *tmap, gcs);
  cout<<"  Created the main lattice grid with"<<endl;
  cout<<"    xlb: "<<grid.xlb().transpose().format(eigformat)<<endl;
  cout<<"    xub: "<<grid.xub().transpose().format(eigformat)<<endl;
  cout<<"    gs: "<<grid.gs().transpose().format(eigformat)<<endl;

  //**********************************************************************************/
  //*************************************cost**************************************/
  //**********************************************************************************/

  //create cost
  SE2GridCostConfig cost_config;
  cost_config.eps = 1e-6;
  if( !params.GetVector3d("wt", cost_config.twist_weight))
    cout<<"param wt missing. Using default values "<<endl;
  shared_ptr<TerrainSE2GridCost> cost(new TerrainSE2GridCost(grid, cost_config));

  //**********************************************************************************/
  //*************************************connectivity*********************************/
  //**********************************************************************************/

  // load car connectivity and set custom parameters
  CarPrimitiveConfig primcfg;
  TerrainTwistConnectivity connectivity(grid, *cost);
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
  cout << "\nPlanner and graph creation" << endl;

  bool initExpand = false;
  params.GetBool("initExpand", initExpand);

  timer_start(&timer);
  GridSearch<TerrainCell::PointType, TerrainCell::DataType, SE2Twist> search(grid, connectivity, *cost, initExpand);
  long time = timer_us(&timer);
  printf("  graph construction time= %ld  us\n", time);
  cout <<"  Created a graph with " << search.Vertices() << " vertices and " << search.Edges() << " edges." << endl;

  //**********************************************************************************/
  //***************************************plan***************************************/
  //**********************************************************************************/

  // get start and goal
  cout << "\nSet start and goal" << endl;
  Vector3d goal, start;
  params.GetVector3d("start", start);
  params.GetVector3d("goal", goal);
  bool start_set = search.SetStart(start);
  bool goal_set = search.SetGoal(goal);
  if(start_set)
    cout<<"  feasible start:"<<start.transpose().format(eigformat)<<endl;
  else
    cout<<"  infeasible start:"<<start.transpose().format(eigformat)<<endl;
  if(goal_set)
    cout<<"  feasible goal:"<<goal.transpose().format(eigformat)<<endl;
  else
    cout<<"  infeasible goal:"<<goal.transpose().format(eigformat)<<endl;


  //Plan the path
  cout << "\nPlanning a path" << endl;
  CarTwistPath path;
  if(start_set && goal_set){
    cout << "  Start and goal feasible" << endl;
    timer_start(&timer);
    search.Plan(path);
    time = timer_us(&timer);
    printf("  plan path time= %ld  us\n", time);
    printf("  path: edges=%lu len=%f\n", path.cells.size(), path.cost);
    cout << "  Graph has " << search.Vertices() << " vertices and " << search.Edges() << " edges. " << endl;
  }else{
    cout<< "  Start or goal or both not feasible. Couldn't run planner"<<endl;
  }

  //**********************************************************************************/
  //***************************************plot***************************************/
  //**********************************************************************************/
  cout << "\nPlotting start goal path and primitives"<<endl;
  vector<Vector3d> path3d;
  if(start_set && goal_set){
    path3d = ToVector3dPath(path,grid.cs()[1]);
    cout << "  Map and path saved to terrain_car_path.ppm" << endl;
  }else{
    path3d.push_back(start); path3d.push_back(goal);
    cout<<"  No planning done so plotting just the start and goal cars to terrain_car_path.ppm"<<endl;
  }
  SavePpmWithPath(*tmap, "terrain_car_path.ppm", 3, path3d, &geom);

  //Plot the primitive at start position
  vector<vector<Vector2d>> prims(0);
  bool gotprims = connectivity.GetPrims(start,prims);
  SavePpmWithPrimitives(*tmap, "terrain_car_prim.ppm", 3, prims);
  cout << "  Map, with primitives from start saved to terrain_car_prim.ppm" << endl;

  return 0;
}
