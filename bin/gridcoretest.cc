#include "grid.h"
#include <Eigen/Dense>
#include <chrono>


using namespace dsl;
using namespace std;
using namespace Eigen;

int main(int argc, char** argv){

  Vector3d xlb = Vector3d::Constant(-0.5);
  Vector3d xub = Vector3d::Constant( 2.5);
  Vector3d cs  = Vector3d::Ones();

  GridCore<Vector3d,double> grid(xlb,xub,cs);
  for(int id = 0; id < grid.nc ; id++){
    grid.cells[id] = id;
  }

  cout<<endl;
  cout<<"contents of the cell of grid"<<endl;
  for(int id = 0; id < grid.nc ; id++){
    cout<<grid.cells[id]<<", ";
  }
  cout<<endl<<endl;

  //Get slice along dimension: dim and index: idx along that dimension
  int dim = 0;
  int idx = 1;

  GridCore<Vector3d,double>::SlicePtr pslice = grid.GetSlice(1,dim);
  cout<<"Slicing along dimension:"<<dim<<" at index:"<<idx<<" from grid"<<endl;
  cout<<"contents of the cell of slice"<<endl;
  for(int id = 0; id < pslice->nc ; id++){
    cout<<pslice->cells[id]<<", ";
  }
  cout<<endl<<endl;

  GridCore<Vector3d,double>::Slice::StackPtr pstack = pslice->GetStack(dim,0,10,3);
  cout<<"Stacked up slices along dim:"<<dim<<" to create stack"<<endl;
  cout<<"contents of the cell of stack"<<endl;
  for(int id = 0; id < pstack->nc ; id++){
    cout<<pstack->cells[id]<<", ";
  }
  cout<<endl<<endl;

  return 0;
}
