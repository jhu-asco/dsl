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
    double d = abs(v(1));
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

vector<Vector2d> ToVector2dPath(const CarTwistPath &path, double cs) {
  vector<Vector3d> path3d = ToVector3dPath(path,cs);
  vector<Vector2d> path2d;
  for(auto& p:path3d)
    path2d.push_back(Vector2d(p(1),p(2)));
  return path2d;
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

  //Use Geometry of car or not
  Vector5d lboxoysb;
  bool use_geom = params.GetVector5d("geom", lboxoysb);
  CarGeom geom;
  if(use_geom)
    geom.set(lboxoysb);

  //load the terrain map
  cout<<"\nLoading terrain map from .tmap file"<<endl;
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


  //Create or load the configuration space occupancy map
  cout<<"\nconfiguration space occupancy map(cmap)"<<endl;
  string cmapfile;
  dsl::Map<bool, 3>::Ptr cmap;
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
    replaceExt(cmapfile, string("cmap"));
    dsl::Map<bool,3>::Save(*cmap, cmapfile);
    cout << "  Saved cmap in "<<cmapfile <<" with"<<endl;
    cout << "    xlb=" << cmap->xlb().transpose().format(eigformat)<<endl;
    cout << "    xub=" << cmap->xub().transpose().format(eigformat)<<endl;
    cout << "    gs="  << cmap->gs().transpose().format(eigformat)<<endl;
  }
  SavePpm(*cmap,"cmap_slices");

  //To plot the path as rectanges or just points
  bool plot_car; params.GetBool("plot_car", plot_car);

  // grid cell size (normally larger than ocs)
  Vector3d gcs;
  params.GetVector3d("gcs", gcs);

  //The main lattice grid

  cout<<"\nMain lattice grid"<<endl;
  TerrainSE2Grid grid(*cmap, *tmap, gcs);
  cout<<"  Created the main lattice grid with"<<endl;
  cout<<"    xlb: "<<grid.xlb().transpose().format(eigformat)<<endl;
  cout<<"    xub: "<<grid.xub().transpose().format(eigformat)<<endl;
  cout<<"    gs: "<<grid.gs().transpose().format(eigformat)<<endl;

  //create cost
  SE2GridCostConfig cost_config;
  cost_config.eps = 1e-6;
  if( !params.GetVector3d("wt", cost_config.twist_weight))
    cout<<"param wt missing. Using default values "<<endl;
  shared_ptr<TerrainSE2GridCost> cost;
  cost.reset(new TerrainSE2GridCost(grid, cost_config));

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


  // create planner
  cout << "\nPlanner and graph creation" << endl;
  bool initExpand = false;
  params.GetBool("initExpand", initExpand);
  timer_start(&timer);
  GridSearch<TerrainCell::PointType, TerrainCell::DataType, SE2Twist> search(grid, connectivity, *cost, initExpand);
  long time = timer_us(&timer);
  printf("  graph construction time= %ld  us\n", time);
  cout <<"  Created a graph with " << search.Vertices() << " vertices and " << search.Edges() << " edges." << endl;

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

  //Plot the start goal and the path
  cout << "\nPlotting start goal path and primitives"<<endl;
  if(start_set && goal_set){
    vector<Vector3d> path3d = ToVector3dPath(path,grid.cs()[1]);
    SavePpmWithPath(*tmap, "terrain_car_path.ppm", 3, path3d, &geom);
    cout << "  Map and path saved to terrain_car_path.ppm" << endl;
  }else{
    cout<<"  No planning done so plotting just the start and goal cars"<<endl;
    vector<Vector3d> path3d; path3d.push_back(start); path3d.push_back(goal);
    SavePpmWithPath(*tmap, "terrain_car_path.ppm", 3, path3d, &geom);
    cout << "  Map, start and goal (no path available ) saved to terrain_car_path.ppm" << endl;
  }

  //Plot the primitive at start position
  vector<vector<Vector2d>> prims(0);
  bool gotprims = connectivity.GetPrims(start,prims);
  SavePpmWithPrimitives(*tmap, "terrain_car_prim.ppm", 3, prims);
  cout << "  Map, with primitives from start saved to terrain_car_prim.ppm" << endl;

  return 0;
}
