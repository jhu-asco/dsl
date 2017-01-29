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


int main(int argc, char** argv)
{
  if (argc!=2) {
    cout << "Usage: $./gridutils params.cfg" << endl;
    cout << "\t\t where params.cfg is a parameter file" << endl;
    cout << "This program will create several slices of the cmap and put them in folder named slices"<<endl;
    cout << "\t\t Make sure to create the folder called slices from where you are runing the executable"<<endl;
    return 0;
  }
  assert(argc == 2);

  Params params(argv[1]);


  // get map
  string map_img_filename;
  if(!params.GetString("map", map_img_filename)){
    cout<<"There is no string parameter named map. Exiting."<<endl;
    return -1;
  }

  // get map
  string cmap_filename;
  bool cmapValid = params.GetString("cmap", cmap_filename);
  
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
  dsl::Map<bool, 2>::Ptr omap = LoadPpm(map_img_filename, ocs.tail<2>());

  // dimensions are determined from occupancy map
  Vector3d xlb(-M_PI + gcs[0]/2, omap->xlb()[0], omap->xlb()[1]);
  Vector3d xub(M_PI + gcs[0]/2, omap->xub()[0], omap->xub()[1]);

  // configuration-space map
  shared_ptr<dsl::Map<bool, 3> > cmap;

  struct timeval timer;
  if (!cmapValid) {
    cmap.reset(new dsl::Map<bool, 3>(xlb, xub, ocs));
    std::cout << "Making cmap... " << std::endl;
    timer_start(&timer);
    if (useGeom) {
      int nthreads; params.GetInt("nthreads", nthreads);
      cmap = MakeCmap(*omap, ocs(0), geom, nthreads);
    } else {
      cmap = MakeCmap(*omap, ocs(0));
    }
    long time = timer_us(&timer);
    printf("cmap construction time= %ld  us\n", time);

    cmap_filename = ReplaceExtension(map_img_filename, string("cmap"));
    cmap->Save(cmap_filename);
    std::cout << "Saved cmap " << cmap_filename << " with xlb=" << cmap->xlb().transpose() << " xub=" << cmap->xub().transpose() << " gs=" << cmap->gs().transpose() << std::endl;

  } else {
    cmap = dsl::Map<bool,3>::Load(cmap_filename);
    std::cout << "Loaded map with xlb=" << cmap->xlb().transpose() << " xub=" << cmap->xub().transpose() << " gs=" << cmap->gs().transpose() << std::endl;
  }

  // saving all slices of the configuration map
  string folder = "slices";
  SavePpm(*cmap,folder);
  return 0;
}
