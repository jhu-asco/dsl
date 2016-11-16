#include <string.h>
#include "gridsearch.h"
#include "cargrid.h"
#include "carcost.h"
#include "carconnectivity.h"
#include "utils.h"
#include "gridutils.h"
#include "params.h"
#include <fstream>

using namespace dsl;
using namespace std;
using namespace Eigen;

using CarPath = GridPath<Vector3d, Matrix3d, SE2Path>;

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

  // get map
  string mapName;
  params.GetString("map", mapName);

  // occupancy cell size
  Vector3d ocs;
  params.GetVector3d("ocs", ocs);  

  // grid cell size (normally larger than ocs)
  Vector3d gcs;
  params.GetVector3d("gcs", gcs);  

  // load an occupancy map from ppm file
  dsl::Map<bool, 2> omap = load(mapName.c_str(), ocs.tail<2>());

  // a map that we'll use for display
  dsl::Map<bool, 2> dmap = omap;

  // dimensions are determined from occupancy map
  Vector3d xlb(-M_PI + gcs[0]/2, omap.xlb[0], omap.xlb[1]);
  Vector3d xub(M_PI + gcs[0]/2, omap.xub[0], omap.xub[1]);

  // configuration-space map
  dsl::Map<bool, 3> cmap(xlb, xub, ocs);
  
  //CarGrid::MakeMap(omap, cmap);
  /* for non-point geometry comment this out */
  double l=3;
  double b=1.5;
  double ox = -0.75;
  double oy = 0;
  CarGeom geom(l, b, ox, oy);
  CarGrid::MakeMap(geom, omap, cmap);  
  
  CarGrid grid(cmap, gcs);

  shared_ptr<dsl::Map<bool,2> > smap;
  cmap.Slice(-M_PI+0.05, 0);  // grid.Slice(cmap, -M_PI+0.05, smap);
  
  // save it to image for viewing
  save(dmap, "path1.ppm");
  save(*smap, "slice0.ppm");

  //  save(dmap, "path1.ppm");


  return 0;
}
