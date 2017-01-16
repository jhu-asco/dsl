#include "grid.h"
#include <Eigen/Dense>
#include <chrono>


using namespace dsl;
using namespace std;
using namespace Eigen;

int main(int argc, char** argv){

  GridCore<Vector3d,Vector3i> grid(Vector3d(-0.5,-0.5,-0.5), Vector3d(99.5,99.5,9.5), Vector3d::Ones().eval());
  for(int id=0; id < grid.nc(); id++){
    Vector3i gidx;
    grid.Index(id, &gidx);
    grid.set_cells(id, gidx);
  }

  cout<<"grid.gs():"<<grid.gs().transpose()<<endl;
  cout<<"grid nc:"<<grid.nc()<<endl;
  cout<<"total memory to save index inside each cell:"<<
      sizeof(std::vector<Vector3i>) + (sizeof(Vector3i) * grid.nc())<<endl;


  chrono::time_point<chrono::system_clock> t_start, t_end;
  chrono::duration<double> elapsed;

  cout<<"iterating over all cells using Index()"<<endl;
  t_start = chrono::system_clock::now();
  for(int id=0; id < grid.nc(); id++){
    Vector3i gidx;
    grid.Index(id, &gidx);
  }
  t_end = chrono::system_clock::now();
  elapsed = t_end - t_start;
  cout<<"dt:"<<elapsed.count()<<endl;

  cout<<"iterating over all cells using saved index saved in each cell"<<endl;
  t_start = chrono::system_clock::now();
  for(int id=0; id < grid.nc(); id++){
    Vector3i idx;
    idx = grid.cells()[id];
  }
  t_end = chrono::system_clock::now();
  elapsed = t_end - t_start;
  cout<<"dt:"<<elapsed.count()<<endl;

  cout<<"iterating over all cells using variable for loop"<<endl;
  t_start = chrono::system_clock::now();
  const int n =3; /*Insert N here: how many loops do you need?*/;
  Vector4i idxe = Vector4i::Zero();
  Vector4i max; max<<grid.gs(),0; //0 to indicate iteration over
  int dim(0);
  int id(0);
  while (idxe[n]==0) {
//    Vector3i idx = idxe.head<3>();
    id++; idxe[0]++;
    while(idxe[dim]==max[dim]) {
      idxe[dim]=0;
      idxe[++dim]++;
      if(idxe[dim]!=max[dim])
        dim=0;
    }
  }
  t_end = chrono::system_clock::now();
  elapsed = t_end - t_start;
  cout<<"dt:"<<elapsed.count()<<endl;


  return 0;
}
