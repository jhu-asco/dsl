#include "grid.h"
#include <Eigen/Dense>


using namespace dsl;
using namespace std;
using namespace Eigen;

using SlicePtr = GridCore<Vector3d,double>::GridCoreSlicePtr ;

int main(int argc, char** argv){

  Vector3d xlb = Vector3d::Constant(-0.5);
  Vector3d xub = Vector3d::Constant( 2.5);
  Vector3d cs  = Vector3d::Ones();

  GridCore<Vector3d,double> grid(xlb,xub,cs);
  for(int id = 0; id < grid.nc ; id++){
    grid.cells[id] = id;
  }

  SlicePtr pslice = grid.Slice(1,0);

  Eigen::Map<MatrixXd > slicemat(pslice->cells.data(), pslice->gs[0], pslice->gs[1]);

  cout<<"slicemat:\n"<<slicemat<<endl;
  return 0;
}
