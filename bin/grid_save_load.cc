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
  using Ptr = std::shared_ptr<HT>;

  HT():h(0),t(0){}

  HT(double h, double t): h(h), t(t){}

  double h;
  double t;
  bool SerializeToOstream(ostream* output) const{

    int double_size = sizeof(double);
    string str_h(double_size,'o'), str_t(double_size,'o');
    str_h.replace(0,double_size,(char*)&h, double_size);
    str_t.replace(0,double_size,(char*)&t, double_size);

    *output << str_h+str_t;
    return true;
  }

  bool ParseFromIstream(istream* input){
    int n_bytes = sizeof(double);
    std::stringstream ss;
    ss << input->rdbuf();
    std::memcpy(&h,ss.str().c_str(), n_bytes);
    std::memcpy(&t,ss.str().c_str()+n_bytes, n_bytes);

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
  cout<<"\nsaving and loading GridCore<Vector3d, bool >"<<endl;
  GridCore<Vector3d,bool >grid(xlb,xub,gs,wd);
  for(int id = 0; id < grid.nc(); id++)
    grid.set_cells(id,true);

  grid.Save(filename);

  auto grid_loaded = GridCore<Vector3d,bool >::Load(filename);

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
  cout<<"\nsaving and loading GridCore<Vector3d,HT>"<<endl;
  GridCore<Vector3d,HT > grid_ht(xlb,xub,gs,wd);
  for(int id = 0; id < grid.nc(); id++)
    grid_ht.set_cells(id,HT(10,20));

  grid_ht.Save(filename);

  auto grid_ht_loaded = GridCore<Vector3d,HT >::Load(filename);

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
  cout<<"\nsaving and loading GridCore<Vector3d,HT::Ptr >"<<endl;
  GridCore<Vector3d,HT::Ptr > grid_pht(xlb,xub,gs,wd);

  //only the 0th and the 2nd cell allocated
  HT::Ptr ptr(new HT(10,20));
  grid_pht.set_cells(0, ptr);
  ptr.reset(new HT(20,30));
  grid_pht.set_cells(2, ptr);

  grid_pht.Save(filename);

  auto grid_pht_loaded = GridCore<Vector3d,HT::Ptr >::Load(filename);

  cout<<"Printing out the loaded grid_pht"<<endl;
  cout<<"xlb:"<< grid_pht_loaded->xlb().transpose() << endl;
  cout<<"xub:"<< grid_pht_loaded->xub().transpose() << endl;
  cout<<"ds:"<< grid_pht_loaded->ds().transpose() << endl;
  cout<<"cs:"<< grid_pht_loaded->cs().transpose() << endl;
  cout<<"gs:"<< grid_pht_loaded->gs().transpose() << endl;
  cout<<"cgs:"<< grid_pht_loaded->cgs().transpose() << endl;
  cout<<"wd:"<< grid_pht_loaded->wd().transpose() << endl;
  for(int id = 0; id< grid_pht_loaded->nc(); id++){
    auto c = grid_pht_loaded->cells()[id];
    if(c)
      cout<<"at id:"<<id<<"  c.h:"<<c->h<<"  c.t:"<<c->t<<endl;
  }


  return 0;
}

