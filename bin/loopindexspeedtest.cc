#include "grid.h"
#include <Eigen/Dense>
#include <chrono>


using namespace dsl;
using namespace std;
using namespace Eigen;

using Vector2b = Matrix<bool,2,1>;
using Vector3b = Matrix<bool,3,1>;
using Vector4b = Matrix<bool,4,1>;

int main(int argc, char** argv){

  {
    //make a 10000x10000 size grid
    Vector2d xlb(-0.5,-0.5);
    Vector2d xub(9999.5,9999.5);
    Vector2d cs(1.0,1.0);
    Vector2b wd(false,false);
    Grid<Vector2d,bool> grid(xlb,xub,cs, wd);

    cout<<"grid.gs:"<<grid.gs().transpose()<<endl;
    cout<<"grid nc:"<<grid.nc()<<endl;

    chrono::time_point<chrono::system_clock> t_start, t_end;
    chrono::duration<double> elapsed;

    //dummy variables
    bool cellval;
    Vector2i gidx_mod;

    cout<<"iterating over all cells using Index()"<<endl;
    t_start = chrono::system_clock::now();
    for(int id=0; id < grid.nc(); id++){
      Vector2i gidx;
      grid.Index(id, &gidx);
      gidx_mod  = 2* gidx; //do something with gidx
      cellval = grid.cells()[id]; //do something with id
    }
    t_end = chrono::system_clock::now();
    elapsed = t_end - t_start;
    cout<<"dt:"<<elapsed.count()<<endl;

    cout<<"iterating over all cells using grid.LoopOver function"<<endl;
    t_start = chrono::system_clock::now();
    auto fun = [&](int id, const Vector2i& gidx){
      gidx_mod  = 2* gidx; //do something with gidx
      cellval = grid.cells()[id]; //do something with id
    };
    grid.LoopOver(fun);
    t_end = chrono::system_clock::now();
    elapsed = t_end - t_start;
    cout<<"dt:"<<elapsed.count()<<endl;

    cout<<"iterating over all cells using nested for loops"<<endl;
    t_start = chrono::system_clock::now();
    int id=0; Vector2i gidx;
      for(int j=0; j < grid.gs()[1]; j++){
        for(int i=0; i < grid.gs()[0]; i++){
          gidx[0] = i; gidx[1] = j;
          gidx_mod  = 2* gidx; //do something with gidx
          cellval = grid.cells()[id]; //do something with id
          id++;
        }
      }
    t_end = chrono::system_clock::now();
    elapsed = t_end - t_start;
    cout<<"dt:"<<elapsed.count()<<endl;
  }

  cout<<endl<<endl;
  {
    //make a 1000x1000x100 size grid
    Vector3d xlb(-0.5,-0.5,-0.5);
    Vector3d xub(999.5,999.5,99.5);
    Vector3d cs(1.0,1.0,1.0);
    Vector3b wd(false,false,false);
    Grid<Vector3d,bool> grid(xlb,xub,cs, wd);

    cout<<"grid.gs:"<<grid.gs().transpose()<<endl;
    cout<<"grid nc:"<<grid.nc()<<endl;

    chrono::time_point<chrono::system_clock> t_start, t_end;
    chrono::duration<double> elapsed;

    //dummy variables
    bool cellval;
    Vector3i gidx_mod;

    cout<<"iterating over all cells using Index()"<<endl;
    t_start = chrono::system_clock::now();
    for(int id=0; id < grid.nc(); id++){
      Vector3i gidx;
      grid.Index(id, &gidx);
      gidx_mod  = 2* gidx; //do something with gidx
      cellval = grid.cells()[id]; //do something with id
    }
    t_end = chrono::system_clock::now();
    elapsed = t_end - t_start;
    cout<<"dt:"<<elapsed.count()<<endl;


    cout<<"iterating over all cells using grid.LoopOver function"<<endl;
    t_start = chrono::system_clock::now();
    auto fun = [&](int id, const Vector3i& gidx){
      gidx_mod  = 2* gidx; //do something with gidx
      cellval = grid.cells()[id]; //do something with id
    };
    grid.LoopOver(fun);
    t_end = chrono::system_clock::now();
    elapsed = t_end - t_start;
    cout<<"dt:"<<elapsed.count()<<endl;

    cout<<"iterating over all cells using nested for loops"<<endl;
    t_start = chrono::system_clock::now();
    int id=0; Vector3i gidx;
      for(int k=0; k < grid.gs()[2]; k++){
        for(int j=0; j < grid.gs()[1]; j++){
          for(int i=0; i < grid.gs()[0]; i++){
            gidx[0] = i; gidx[1] = j; gidx[2]=k;
            gidx_mod  = 2* gidx; //do something with gidx
            cellval = grid.cells()[id]; //do something with id
            id++;
          }
        }
      }
    t_end = chrono::system_clock::now();
    elapsed = t_end - t_start;
    cout<<"dt:"<<elapsed.count()<<endl;
  }

  cout<<endl<<endl;
  {
    //make a 1000x100x100x10 size grid
    Vector4d xlb(-0.5,-0.5,-0.5,-0.5);
    Vector4d xub(999.5,99.5,99.5,9.5);
    Vector4d cs(1.0,1.0,1.0,1.0);
    Vector4b wd(false,false,false,false);
    Grid<Vector4d,bool> grid(xlb,xub,cs, wd);

    cout<<"grid.gs:"<<grid.gs().transpose()<<endl;
    cout<<"grid nc:"<<grid.nc()<<endl;

    chrono::time_point<chrono::system_clock> t_start, t_end;
    chrono::duration<double> elapsed;

    //dummy variables
    bool cellval;
    Vector4i gidx_mod;

    cout<<"iterating over all cells using Index()"<<endl;
    t_start = chrono::system_clock::now();
    for(int id=0; id < grid.nc(); id++){
      Vector4i gidx;
      grid.Index(id, &gidx);
      gidx_mod  = 2* gidx; //do something with gidx
      cellval = grid.cells()[id]; //do something with id
    }
    t_end = chrono::system_clock::now();
    elapsed = t_end - t_start;
    cout<<"dt:"<<elapsed.count()<<endl;

    cout<<"iterating over all cells using grid.LoopOver function"<<endl;
    t_start = chrono::system_clock::now();
    auto fun = [&](int id, const Vector4i& gidx){
      gidx_mod  = 2* gidx; //do something with gidx
      cellval = grid.cells()[id]; //do something with id
    };
    grid.LoopOver(fun);
    t_end = chrono::system_clock::now();
    elapsed = t_end - t_start;
    cout<<"dt:"<<elapsed.count()<<endl;

    cout<<"iterating over all cells using nested for loops"<<endl;
    t_start = chrono::system_clock::now();
    int id=0; Vector4i gidx;
    for(int l=0; l < grid.gs()[3]; l++){
      for(int k=0; k < grid.gs()[2]; k++){
        for(int j=0; j < grid.gs()[1]; j++){
          for(int i=0; i < grid.gs()[0]; i++){
            //gidx[0] = i; gidx[1]=j; gidx[2]=k; gidx[3]=l;
            gidx << i,j,k,l;
            gidx_mod  = 2* gidx; //do something with gidx
            cellval = grid.cells()[id]; //do something with id
            id++;
          }
        }
      }
    }
    t_end = chrono::system_clock::now();
    elapsed = t_end - t_start;
    cout<<"dt:"<<elapsed.count()<<endl;
  }
  return 0;
}
