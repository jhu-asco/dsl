/*
 * basic.cpp
 *
 *  Created on: Jan 3, 2017
 *      Author: subhransu
 */
#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <vector>
#include "grid.h"
#include <memory>

using namespace Eigen;
using namespace std;

using Vector3b = Eigen::Matrix<bool,3,1>;

struct HT{
  using Ptr = std::unique_ptr<HT>;
  using Cref = const HT&;
  HT():h(0),t(0){}

  HT(double h, double t): h(h), t(t){}

  double h;
  double t;

  bool SerializeToString(std::string* str) const{
    int double_size = sizeof(double);
    str->resize(2*double_size);
    str->replace(0,double_size,(char*)&h, double_size);
    str->replace(double_size,double_size,(char*)&t, double_size);
    return true;
  }

  bool ParseFromString(const std::string& str){
    int n_bytes = sizeof(double);
    std::memcpy(&h,str.c_str(), n_bytes);
    std::memcpy(&t,str.c_str()+n_bytes, n_bytes);
    return true;
  }
};


using namespace dsl;
using namespace std;
int main(int argc, char** argv){

  string filename("grid_data.bin");
  Vector3d xlb(0, 1, 2);
  Vector3d xub(10, 11, 12);
  Vector3i gs(2,1,2); //total of 4 cells
  Vector3b wd(true,false,false);

  //**********************************************************
  //*********************bool*********************************
  //**********************************************************
  cout<<"\nsaving and loading Grid<Vector3d, bool >"<<endl;
  Grid<Vector3d, bool, false>grid(xlb,xub,gs,wd);
  for(int id = 0; id < grid.nc(); id++)
    grid.Set(id,true);

  grid.Save(filename);

  auto grid_loaded = Grid<Vector3d,bool,false >::Load(filename);

  cout<<"Printing out the loaded grid for bools"<<endl;
  cout<<"xlb:"<< grid_loaded->xlb().transpose() << endl;
  cout<<"xub:"<< grid_loaded->xub().transpose() << endl;
  cout<<"ds:"<< grid_loaded->ds().transpose() << endl;
  cout<<"cs:"<< grid_loaded->cs().transpose() << endl;
  cout<<"gs:"<< grid_loaded->gs().transpose() << endl;
  cout<<"cgs:"<< grid_loaded->cgs().transpose() << endl;
  cout<<"wd:"<< grid_loaded->wd().transpose() << endl;
  for(auto&& c: grid_loaded->cells())
    cout<< c<<endl;


  //**********************************************************
  //*********************HT***********************************
  //**********************************************************
  cout<<"\nsaving and loading Grid<Vector3d,HT>"<<endl;
  Grid<Vector3d, HT, false > grid_ht(xlb,xub,gs,wd);
  for(int id = 0; id < grid.nc(); id++)
    grid_ht.Set(id,HT(10,20));

  grid_ht.Save(filename);

  auto grid_ht_loaded = Grid<Vector3d,HT, false >::Load(filename);

  cout<<"Printing out the loaded grid_ht"<<endl;
  cout<<"xlb:"<< grid_ht_loaded->xlb().transpose() << endl;
  cout<<"xub:"<< grid_ht_loaded->xub().transpose() << endl;
  cout<<"ds:"<< grid_ht_loaded->ds().transpose() << endl;
  cout<<"cs:"<< grid_ht_loaded->cs().transpose() << endl;
  cout<<"gs:"<< grid_ht_loaded->gs().transpose() << endl;
  cout<<"cgs:"<< grid_ht_loaded->cgs().transpose() << endl;
  cout<<"wd:"<< grid_ht_loaded->wd().transpose() << endl;
  for(auto&& c: grid_ht_loaded->cells()){
    cout<<"c.h:"<<c.h<<"  c.t:"<<c.t<<endl;
  }

  //**********************************************************
  //*********************HT::Ptr******************************
  //**********************************************************
  cout<<"\nsaving and loading Grid<Vector3d,HT::Ptr >"<<endl;
  Grid<Vector3d,HT, true> grid_pht(xlb,xub,gs,wd);

  //only the 0th and the 2nd cell allocated
  grid_pht.Set(0, HT(10,20));
  grid_pht.Set(2, HT(20,30));

  grid_pht.Save(filename);

  auto grid_pht_loaded = Grid<Vector3d,HT, true >::Load(filename);

  cout<<"Printing out the loaded grid_pht"<<endl;
  cout<<"xlb:"<< grid_pht_loaded->xlb().transpose() << endl;
  cout<<"xub:"<< grid_pht_loaded->xub().transpose() << endl;
  cout<<"ds:"<< grid_pht_loaded->ds().transpose() << endl;
  cout<<"cs:"<< grid_pht_loaded->cs().transpose() << endl;
  cout<<"gs:"<< grid_pht_loaded->gs().transpose() << endl;
  cout<<"cgs:"<< grid_pht_loaded->cgs().transpose() << endl;
  cout<<"wd:"<< grid_pht_loaded->wd().transpose() << endl;
  for(int id = 0; id< grid_pht_loaded->nc(); id++){
    HT::Cref c = grid_pht_loaded->Get(id);
    if(&c)
      cout<<"at id:"<<id<<"  c.h:"<<c.h<<"  c.t:"<<c.t<<endl;
  }


  return 0;
}


