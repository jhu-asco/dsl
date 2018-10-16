#ifndef GCOP_PARAMS_H
#define GCOP_PARAMS_H

#include <map>
#include <cstring>
#include <iostream>
#include <Eigen/Dense>

/*
 *  Provides storage and retrieval for arbitrary 
 *  number of parameters of different types (currently only primitives and strings)
 *  TODO: include matrix parameter
 *  Parameters can be loaded/saved from/to a file, or can be inserted at run-time 
 * 
 *
 *  Author: Marin Kobilarov (mkobilar@@robotics.usc.edu)
 */

namespace dsl {

#define DSL_PARAMS_MBL 256
  
  typedef Eigen::Matrix<double, 6, 6> Matrix6d;
  typedef Eigen::Matrix<double, 6, 1> Vector6d;
  typedef Eigen::Matrix<double, 5, 5> Matrix5d;
  typedef Eigen::Matrix<double, 5, 1> Vector5d;


class Params
{
 public:

  Params();
  Params(FILE *file);
  Params(const char *fileName);
  Params(std::iostream &io);

  virtual ~Params();

  void load(const char *fileName);
  void load(FILE *file);
  void load(std::iostream &io);

  void save(const char *fileName) const;
  void save(FILE* file) const;
  void save(std::iostream &io) const;

  bool exists(const char *name) const { return valueMap.find(std::string(name)) != valueMap.end(); }

  void setInt(const char* name, int v);
  bool getInt(const char* name, int& v) const;

  void setFloat(const char* name, float v);
  bool getFloat(const char* name, float& v) const;

  void setDouble(const char* name, double v);
  bool getDouble(const char* name, double& v) const;

  void setVectorXd(const char* name, const Eigen::VectorXd& v);
  bool getVectorXd(const char* name, Eigen::VectorXd& v) const;

  void setVector2d(const char* name, const Eigen::Vector2d& v);
  bool getVector2d(const char* name, Eigen::Vector2d& v) const;

  void setVector3d(const char* name, const Eigen::Vector3d& v);
  bool getVector3d(const char* name, Eigen::Vector3d& v) const;

  void setVector4d(const char* name, const Eigen::Vector4d& v);
  bool getVector4d(const char* name, Eigen::Vector4d& v) const;

  void setVector5d(const char* name, const Vector5d& v);
  bool getVector5d(const char* name, Vector5d& v) const;

  void setVector6d(const char* name, const Vector6d& v);
  bool getVector6d(const char* name, Vector6d& v) const;

  void setMatrix3d(const char* name, const Eigen::Matrix3d& m);
  bool getMatrix3d(const char* name, Eigen::Matrix3d& m) const;

  void setMatrix6d(const char* name, const Matrix6d& m);
  bool getMatrix6d(const char* name, Matrix6d& m) const;

  void setDoubleVec(const char* name, const std::vector< double >& v);
  bool getDoubleVec(const char* name, std::vector< double >& v) const;

  void setFloatArray(const char* name, int n, const float* v);
  bool getFloatArray(const char* name, int n, float* v) const;

  void setDoubleArray(const char* name, int n, const double* v);
  bool getDoubleArray(const char* name, int n, double* v) const;

  void setString(const char* name, const std::string& v);
  bool getString(const char* name, std::string& v) const;

  void setBool(const char* name, bool v);
  bool getBool(const char* name, bool& v) const;

  void setChar(const char* name, char v);
  char getChar(const char* name, char& v) const;

  void print(FILE* file = stdout) const;
  void print(std::iostream& io) const;

 protected:
   struct removeDelimiter {
    bool operator()(char c)
    {
      return (c =='\r' || c =='\t' || c == ' ' || c == '\n');
    }
  };

  void parse(char* line);

  std::map<std::string, std::string> valueMap;
  char buf[DSL_PARAMS_MBL];
};

}

#endif
