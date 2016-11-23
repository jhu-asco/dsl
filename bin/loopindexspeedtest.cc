#include "grid.h"
#include <Eigen/Dense>
#include <chrono>


using namespace dsl;
using namespace std;
using namespace Eigen;

using Vector3b = Matrix<bool,3,1>;

int main(int argc, char** argv){

  //make a 1000x1000x100 size grid with fast index disabled
  Vector3d xlb(-0.5,-0.5,-0.5);
  Vector3d xub(999.5,999.5,99.5);
  Vector3d cs(1.0,1.0,1.0);
  Vector3b wd(false,false,false);
  bool fastindex = false;
  GridCore<Vector3d,bool> grid(xlb,xub,cs, wd,fastindex);

  cout<<"grid.gs:"<<grid.gs.transpose()<<endl;
  cout<<"grid nc:"<<grid.nc<<endl;

  chrono::time_point<chrono::system_clock> t_start, t_end;
  chrono::duration<double> elapsed;

  //dummy variables
  bool cellval;
  Vector3i gidx_mod;

  cout<<"iterating over all cells using Index()"<<endl;
  t_start = chrono::system_clock::now();
  for(int id=0; id < grid.nc; id++){
    Vector3i gidx; grid.Index(gidx,id);
    gidx_mod  = 2* gidx; //do something with gidx
    cellval = grid.cells[id]; //do something with id
  }
  t_end = chrono::system_clock::now();
  elapsed = t_end - t_start;
  cout<<"dt:"<<elapsed.count()<<endl;

  cout<<"iterating over all cells using grid.LoopOver function"<<endl;
  t_start = chrono::system_clock::now();
  auto fun = [&](int id, const Vector3i& gidx){
    gidx_mod  = 2* gidx; //do something with gidx
    cellval = grid.cells[id]; //do something with id
  };
  grid.LoopOver(fun);
  t_end = chrono::system_clock::now();
  elapsed = t_end - t_start;
  cout<<"dt:"<<elapsed.count()<<endl;

  cout<<"iterating over all cells using nested for loops"<<endl;
  t_start = chrono::system_clock::now();
  for(int k=0; k < grid.gs[2]; k++){
    for(int j=0; j < grid.gs[1]; j++){
      for(int i=0; i < grid.gs[0]; i++){
        Vector3i gidx(i,j,k);
        int id = grid.Id(gidx);
        gidx_mod  = 2* gidx; //do something with gidx
        cellval = grid.cells[id]; //do something with id
      }
    }
  }
  t_end = chrono::system_clock::now();
  elapsed = t_end - t_start;
  cout<<"dt:"<<elapsed.count()<<endl;

  return 0;
}
