#include "cargrid.h"
#include "utils.h"

#include <iostream>

using namespace dsl;
using namespace Eigen;
using namespace std;
typedef Transform<double,2,Affine> Transform2d;
typedef Matrix<double,2,4> Matrix2x4d;
typedef Matrix<int,2,4> Matrix2x4i;

void getRotdVertsInPixWrtOrg(Matrix2x4d& verts2d_rotd_pix, double l,double b, double ox, double oy, double sx, double sy, double theta)
{
  Vector2d org_m(-ox,-oy);
  Matrix2x4d verts_wrt_center;
  verts_wrt_center.col(0) << -l/2, -b/2;//rb(right back)
  verts_wrt_center.col(1) <<  l/2, -b/2;//rf(right front)
  verts_wrt_center.col(2) <<  l/2,  b/2;//lf(left front)
  verts_wrt_center.col(3) << -l/2,  b/2;//lb(left back)

  //rotate vertices about origin and convert vertices and origin to pixel coordinates
  Transform2d tfm2d= Rotation2Dd(theta)*Translation2d(org_m);
  verts2d_rotd_pix = Scaling(1/sx,1/sy)* (tfm2d*verts_wrt_center);
}

template< typename T, int m, int n>
int inPoly(Matrix<T,m,n> verts, Matrix<T,m,1> pt)
{
  int i, j, c = 0;
    for (i = 0, j = n-1; i < n; j = i++) {
      if ( ((verts(1,i)>pt(1)) != (verts(1,j)>pt(1))) &&
       (pt(0) < (verts(0,j)-verts(0,i)) * (pt(1)-verts(1,i)) / (verts(1,j)-verts(1,i)) + verts(0,i)) )
         c = !c;
    }
    return c;
}

void fillQuad(double* data,int w, int h, Matrix2x4d  verts, double val)
{

  for(int r=0;r<h;r++)
  {
    for(int c=0;c<w;c++)
    {
      int idx_2d = c + r*w;
      if( inPoly(verts,Vector2d(c,r)))
        data[idx_2d] = val;
      else
        data[idx_2d] = 0;
    }
  }
}

template<typename T>
void myDilate(T* data_dil, T* data, int w, int h, T* data_k, int w_k, int h_k, int ox_k, int oy_k)
{
  vector<double> prod(h_k*w_k);

  //Visit each pixel in the main image
  for(int r=0; r<h; r++)
  {
    int r_rel = r - oy_k;
    for(int c=0; c <w; c++)
    {
      int c_rel = c - ox_k;
      int id=c+r*w;
      //for each pixel in the main image lay the kernel on the original image
      //  such that the origin of kernel(only 1s and 0s) coincides with pixel in question
      //  Then for pixel in question visit all the elements in the kernel. Each element of kernel image
      //  has a corresponding pixel in main image(except at the borders). Multiply those pair of values
      //  together. Take the max of all those values and that becomes the pixel value of the dilated image.
      //  What this does is checks if the any element of kernel is overlaid over an obstacle.
      for(int r_k=0; r_k < h_k; r_k++)
      {
        for(int c_k=0; c_k < w_k; c_k++)
        {
          int id_k=c_k+r_k*w_k;
          int id_roi = c_k+c_rel + (r_k+r_rel)*w;
          if(c_k + c_rel>=0 && c_k + c_rel<w && r_k+r_rel>=0 && r_k+r_rel<h  )
            prod[id_k] = data_k[id_k] * data[id_roi];
          else
            prod[id_k] = 0;
        }
      }
      data_dil[id] = *(max_element<typename vector<T>::iterator> (prod.begin(),prod.end()));

    }
  }
}


void
CarGrid::getDilatedMap(double* data_dil, double* data, double theta)
{
  double sx=cs(1), sy=cs(2);
  Matrix2x4d verts2d_rotd_pix;
  getRotdVertsInPixWrtOrg(verts2d_rotd_pix, l_, b_, ox_, oy_, sx, sy,theta);

  //round of the pixel values of the vertices above such that the rectange formed by the rounded off
  //  vertices surrounds the rotated rectange
  Vector2i org2i_rotd_pix; org2i_rotd_pix.setZero();//because it's wrt org itself and rounding doesn't matter
  Matrix2x4i verts2i_rotd_pix;
  for(int i=0;i<2;i++)
    for(int j=0;j<4;j++)
      verts2i_rotd_pix(i,j) = verts2d_rotd_pix(i,j)>0?ceil(verts2d_rotd_pix(i,j)):floor(verts2d_rotd_pix(i,j));

  //Size of kernel is given by the horizontal rectange that bounds the rounded off rectange above
  Vector2i verts2i_rotd_pix_min = verts2i_rotd_pix.rowwise().minCoeff();
  Vector2i verts2i_rotd_pix_max = verts2i_rotd_pix.rowwise().maxCoeff();
  Vector2i size2i_k = verts2i_rotd_pix_max - verts2i_rotd_pix_min;

  //Shift everything such that verts2i_rotd_pix_min is the [0,0] pixel of the kernel
  Matrix2x4i verts2i_rotd_pospix = verts2i_rotd_pix.colwise() - verts2i_rotd_pix_min;
  Vector2i   org2i_rotd_pospix   = org2i_rotd_pix - verts2i_rotd_pix_min;

  //create dilation kernel by filling the inside of the rotated rectanges with zero
  int w_k=size2i_k(0);
  int h_k=size2i_k(1);
  double data_k[w_k*h_k];
  fillQuad(data_k,size2i_k(0), size2i_k(1),verts2i_rotd_pospix.cast<double>(),1.0);

  //Dilate
  myDilate(data_dil, data, gs[1], gs[2], data_k, size2i_k(0), size2i_k(1),
           org2i_rotd_pospix(0), org2i_rotd_pospix(1));
}


CarGrid::CarGrid(double l,double b, double ox, double oy,
                 int width, int height, double* map, double sx, double sy, double sa,
                 double costScale,double maxCost)
                 :Grid<3, Matrix3d>(Vector3d(-M_PI/* + sa/2*/, 0, 0),
                                    Vector3d(M_PI/* + sa/2*/, sx*width, sy*height),
                                    Vector3i((int)round(2*M_PI/sa), width, height))
                 ,maxCost(maxCost),l_(l),b_(b),ox_(ox),oy_(oy)
{

  const int &angRes = gs[0];

  for (int k = 0; k < angRes; ++k)
  {
    //create a dilated map for a particular angle
    double theta = xlb(0) + (k + 0.5)*sa;
    double map_dil[width*height];
    getDilatedMap(map_dil,map,theta);

    for (int c = 0; c < width; ++c)
    {
      for (int r = 0; r < height; ++r)
      {
        int idx_2d = r*width + c;// since data is in row major format
        int idx_3d = r*angRes*width + c*angRes + k; //1,2 and 3 dim are a,x and y respectively

        double cost = costScale * map_dil[idx_2d];//Cell cost based on angle and geometry of car

        // add this as a cell only if cost is less than a given max cost
        // this is useful if maxCost defines map cells that are untreversable, so
        // they shouldn't be added to the list of cells
        if (cost < maxCost)
        {
          cells[idx_3d] = new SE2Cell(xlb + Vector3d((k + 0.5)*sa, (c + 0.5)*sx, (r + 0.5)*sy),
                                     Vector3d(sa/2, sx/2, sy/2), cost);
          se2_q2g(cells[idx_3d]->data, cells[idx_3d]->c);
        }
      }
    }
  }

}



CarGrid::CarGrid(int width, int height, double *map,
                 double sx, double sy, double sa, double costScale,
                 double maxCost) :
  Grid<3, Matrix3d>(Vector3d(-M_PI, 0, 0), Vector3d(M_PI, sx*width, sy*height),
                    Vector3i((int)round(2*M_PI/sa), width, height)),
                    maxCost(maxCost),l_(0),b_(0),ox_(0),oy_(0) {

  const int &angRes = gs[0];
  /*
  for (int i = 0; i < width; ++i) {
    for (int j = 0; j < height; ++j) {
      int mid = j*width + i;
      double cost = map[mid]*costScale; // cell cost = height/occupany/traversability
      for (int k = 0; k < angRes; ++k) {
        int id = k*width*height + j*width + i;
        // add this as a cell only if cost is less than a given max cost
        // this is useful if maxCost defines map cells that are untreversable, so
        // they shouldn't be added to the list of cells
        if (cost < maxCost) {
          cells[id] = new Cell<3>(xlb + Vector3d((k + 0.5)*sa, (i + 0.5)*sx, (j + 0.5)*sy),
                                  Vector3d(sa/2, sx/2, sy/2), cost);
        }
      }
    }
  }
  */

  for (int c = 0; c < width; ++c)
  {
    for (int r = 0; r < height; ++r)
    {
      int mid = r*width + c;
      //The pixel value in
      double cost = map[mid]*costScale; // cell cost = height/occupany/traversability
      for (int k = 0; k < angRes; ++k)
      {
        int id = r*angRes*width + c*angRes + k;

        // add this as a cell only if cost is less than a given max cost
        // this is useful if maxCost defines map cells that are untreversable, so
        // they shouldn't be added to the list of cells
        if (cost < maxCost)
        {
          cells[id] = new SE2Cell(xlb + Vector3d((k + 0.5)*sa, (c + 0.5)*sx, (r + 0.5)*sy),
                                  Vector3d(sa/2, sx/2, sy/2), cost);
          se2_q2g(cells[id]->data, cells[id]->c);


        }
      }
    }
  }

}

