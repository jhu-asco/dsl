#include <string.h>
#include "gridsearch.h"
#include "cargrid.h"
#include "carcost.h"
#include "carconntwistprim.h"
#include "utils.h"
#include "gridutils.h"
#include "params.h"
#include <fstream>

using namespace dsl;
using namespace std;
using namespace Eigen;

using CarPathTwist = dsl::GridPath<Vector3d, Matrix3d, dsl::SE2Prim>;

vector<Vector3d> ToVector3dPath(const CarPathTwist &path, double gridcs) {
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

vector<Vector2d> ToVector2dPath(const CarPathTwist &path, double cs) {
  vector<Vector3d> path3d = ToVector3dPath(path,cs);
  vector<Vector2d> path2d;
  for(auto& p:path3d)
    path2d.push_back(Vector2d(p(1),p(2)));
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

  //To plot the path as rectanges or just points
  bool plot_car; params.GetBool("plot_car", plot_car);

  // get start and goal
  Vector3d goal, start;
  params.GetVector3d("start", start);
  params.GetVector3d("goal", goal);

  // get name of 2d(in x and y) occupancy map
  string mapName;
  if(!params.GetString("map", mapName)){
    cout<<"There is no string parameter named map. Exiting."<<endl;
    return -1;
  }

  // get name of 3d(in angle, x and y) occupancy map
  string cmapName;
  bool cmapValid = params.GetString("cmap", cmapName);
  
  // occupancy cell size
  Vector3d ocs;
  params.GetVector3d("ocs", ocs);

  // grid cell size (normally larger than ocs)
  Vector3d gcs;
  params.GetVector3d("gcs", gcs);

  Vector5d whxy;
  bool useGeom = params.GetVector5d("geom", whxy);
  CarGeom geom;
  if(useGeom)
    geom.set(whxy(0),whxy(1),whxy(2),whxy(3),whxy(4));
  
  // load an occupancy map from ppm file
  dsl::Map<bool, 2> omap = load(mapName.c_str(), ocs.tail<2>());

  // a map that we'll use for display
  dsl::Map<bool, 2> dmap = omap;

  // dimensions are determined from occupancy map
  Vector3d xlb(-M_PI + gcs[0]/2, omap.xlb[0], omap.xlb[1]);
  Vector3d xub(M_PI + gcs[0]/2, omap.xub[0], omap.xub[1]);

  // configuration-space map
  shared_ptr<dsl::Map<bool, 3>> cmap = 0;

  //  CarGrid::MakeMap(omap, cmap);
  // for non-point geometry comment this out
  struct timeval timer;
  if (!cmapValid) {
    cmap.reset(new dsl::Map<bool, 3>(xlb, xub, ocs));
    std::cout << "Making cmap... " << std::endl;
    timer_start(&timer);
    if (useGeom) {
      CarGrid::MakeMap(geom, omap, *cmap);
    } else {
      CarGrid::MakeMap(omap, *cmap);
    }
    long time = timer_us(&timer);
    printf("cmap construction time= %ld  us\n", time);

    cmapName = mapName;
    replaceExt(cmapName, string("cmap"));
    dsl::Map<bool,3>::Save(*cmap, cmapName);
    std::cout << "Saved cmap " << cmapName << " with xlb=" << cmap->xlb.transpose() << " xub=" << cmap->xub.transpose() << " gs=" << cmap->gs.transpose() << std::endl;

  } else {
    cmap = dsl::Map<bool,3>::Load(cmapName);
    std::cout << "Loaded map with xlb=" << cmap->xlb.transpose() << " xub=" << cmap->xub.transpose() << " gs=" << cmap->gs.transpose() << std::endl;
  }

  CarGrid grid(*cmap, gcs);

  // create cost and set custom angular cost mixing factor ac
  CarCost cost;
  params.GetDouble("ac", cost.ac);

  // load car connectivity and set custom parameters
  CarPrimitiveCfg primcfg;
  int nl, na;
  params.GetBool("prim_fwdonly",primcfg.fwdonly);
  params.GetDouble("prim_tphioverlmax",primcfg.tphioverlmax);
  params.GetDouble("prim_lmin",primcfg.lmin);
  params.GetDouble("prim_lmax",primcfg.lmax);
  params.GetInt("prim_nl",nl); primcfg.nl = nl;
  params.GetDouble("prim_amax",primcfg.amax);
  params.GetInt("prim_na",na);primcfg.na = na;
  params.GetBool("prim_pert",primcfg.pert);
  params.GetBool("prim_tocenter",primcfg.tocenter);

  CarConnTwistPrim connectivity(grid, primcfg);

  cout << "Creating a graph..." << endl;
  // create planner

  timer_start(&timer);

  bool initExpand = false;
  params.GetBool("initExpand", initExpand);


  GridSearch<Vector3d, Matrix3d, SE2Prim> search(grid, connectivity, cost, initExpand);
  CarPathTwist path;

  long time = timer_us(&timer);
  printf("graph construction time= %ld  us\n", time);

  bool start_set = search.SetStart(start);
  bool goal_set = search.SetGoal(goal);

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

    // save it to image for viewing
    if(plot_car){
      vector<Vector3d> path3d = ToVector3dPath(path,grid.cs[1]);
      saveMapWithPath(dmap, "path2.ppm", path3d, geom, 3);
    }else{
      vector<Vector2d> path2d = ToVector2dPath(path,grid.cs[1]);
      save(dmap, "path2.ppm", &path2d);
    }
    cout << "Map and path saved to path2.ppm" << endl;
  }else{

    if(plot_car){
      vector<Vector3d> path3d; path3d.push_back(start); path3d.push_back(goal);
      saveMapWithPath(dmap, "path2.ppm", path3d, geom, 3);
    }else{

      Vector2d start2d = start.tail<2>(); Vector2d goal2d = goal.tail<2>();
      vector<Vector2d> path2d; path2d.push_back(start2d); path2d.push_back(goal2d);
      save(dmap, "path2.ppm", &path2d);
    }
    cout << "Map, start and goal (no path available ) saved to path2.ppm" << endl;
  }

  //Display the primitive at start
  vector<vector<Vector2d>> prims(0);
  bool gotprims = connectivity.GetPrims(start,prims);
  if(gotprims){
    saveMapWithPrims(dmap, "prim2.ppm",prims,3);
    cout << "Map, with primitives from start saved to prim2.ppm" << endl;
  }

  cmap.reset();
  return 0;
}
