#include <string.h>
#include "gridsearch.h"
#include "cargrid.h"
#include "carcost.h"
#include "carconnectivity.h"
#include "utils.h"
#include "params.h"
#include <fstream>

using namespace dsl;
using namespace std;
using namespace Eigen;

using CarPath = GridPath<Vector3d, Matrix3d, SE2Path>;

vector<Vector2d> ToVector2dPath(const CarPath &path) {
  vector<Vector2d> path2d;
  for (auto&& connection : path.connections)
    for (auto&& g : connection)      
      path2d.push_back(g.topRightCorner<2,1>());
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

  // get start and goal
  Vector3d goal, start;
  params.GetVector3d("start", start);
  params.GetVector3d("goal", goal);

  // get map
  string mapName;
  params.GetString("map", mapName);

  // get map
  string cmapName;
  bool cmapValid = params.GetString("cmap", cmapName);
  
  // occupancy cell size
  Vector3d ocs;
  params.GetVector3d("ocs", ocs);

  // grid cell size (normally larger than ocs)
  Vector3d gcs;
  params.GetVector3d("gcs", gcs);

  Vector4d whxy;
  bool useGeom = params.GetVector4d("geom", whxy);
  
  // load an occupancy map from ppm file
  dsl::Map<bool, 2> omap = load(mapName.c_str(), ocs.tail<2>());

  // a map that we'll use for display
  dsl::Map<bool, 2> dmap = omap;

  // dimensions are determined from occupancy map
  Vector3d xlb(-M_PI + gcs[0]/2, omap.xlb[0], omap.xlb[1]);
  Vector3d xub(M_PI + gcs[0]/2, omap.xub[0], omap.xub[1]);

  // configuration-space map
  dsl::Map<bool, 3> *cmap = 0;
  
  //  CarGrid::MakeMap(omap, cmap);
  // for non-point geometry comment this out 
  if (!cmapValid) {
    cmap = new dsl::Map<bool, 3>(xlb, xub, ocs);
    std::cout << "Making cmap... " << std::endl;

    if (useGeom) {
      CarGeom geom(whxy(0), whxy(1), whxy(2), whxy(3));
      CarGrid::MakeMap(geom, omap, *cmap);
    } else {
      CarGrid::MakeMap(omap, *cmap);      
    }
    
    cmapName = mapName;
    replaceExt(cmapName, string("cmap"));
    dsl::Map<bool,3>::Save(*cmap, cmapName.c_str());
    std::cout << "Saved cmap " << cmapName << " with xlb=" << cmap->xlb.transpose() << " xub=" << cmap->xub.transpose() << " gs=" << cmap->gs.transpose() << std::endl;
    
  } else {
    cmap = dsl::Map<bool,3>::Load(cmapName.c_str());
    std::cout << "Loaded map with xlb=" << cmap->xlb.transpose() << " xub=" << cmap->xub.transpose() << " gs=" << cmap->gs.transpose() << std::endl;
  }
  
  CarGrid grid(*cmap, gcs);

  // create cost and set custom angular cost mixing factor ac
  CarCost cost;
  params.GetDouble("ac", cost.ac);    

  // load car connectivity and set custom parameters
  CarConnectivity connectivity(grid);
  double dt = .25;
  double vx = 4;
  double kmax = 0.57;
  int kseg = 4;
  bool onlyfwd = false;
  params.GetDouble("dt", dt);  
  params.GetDouble("vx", vx);
  params.GetDouble("kmax", kmax);
  params.GetInt("kseg", kseg);
  params.GetBool("onlyfwd", onlyfwd);
  connectivity.SetPrimitives(dt, vx, kmax, kseg, onlyfwd);

  cout << "Creating a graph..." << endl;
  // create planner
  struct timeval timer;
  timer_start(&timer);

  bool initExpand = false;
  params.GetBool("initExpand", initExpand);

  
  GridSearch<Vector3d, Matrix3d, SE2Path> search(grid, connectivity, cost, initExpand);
  CarPath path;

  long time = timer_us(&timer);
  printf("graph construction time= %ld  us\n", time);

  search.SetStart(start);
  search.AddGoal(goal);

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
  vector<Vector2d> path2d = ToVector2dPath(path);
  save(dmap, "path1.ppm", &path2d);

  cout << "Map and path saved to path1.ppm" << endl;

  delete cmap;
  return 0;
}
