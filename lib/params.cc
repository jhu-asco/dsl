#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include "params.h"


namespace dsl {

using namespace std;
using namespace Eigen;

static void replace(std::string& str, const std::string& a, const std::string& b)
{
  size_t pos = 0;
  while((pos = str.find(a, pos)) != std::string::npos) {
    str.replace(pos, a.length(), b);
    pos += b.length();
  }
}

template <typename T>
static void process(const string& str, T& v) {
  if (str.find(string("pi")) != std::string::npos) {
    v = M_PI;
  } else {
    v = atof(str.c_str());
  }
}

static void tokenize(const string& str,
                     vector< string >& tokens,
                     const string& delimiters = " ") {
  // Skip delimiters at beginning.
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // find first "non-delimiter".
  string::size_type pos     = str.find_first_of(delimiters, lastPos);
  
  while (string::npos != pos || string::npos != lastPos)
    {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of(delimiters, pos);
      // find next "non-delimiter"
      pos = str.find_first_of(delimiters, lastPos);
    }
}




Params::Params()
{
}

Params::Params(const char *fileName)
{
  load(fileName);
}

Params::Params(FILE *file)
{
  load(file);
}

Params::Params(iostream &io)
{
  load(io);
}


Params::~Params()
{
  valueMap.clear();
}

void Params::load(const char *fileName)
{
  FILE *file = fopen(fileName, "r");
  if (!file) {
    cout << "[W] Params::Load: failed to load file " << fileName << endl;
    return;
  }
  load(file);
  fclose(file);
}

void Params::parse(char* line) {
  line[strcspn(line, "\n\r\t")] = 0;    
  if (line[0] == '#' || strlen(line) <= 2)
    return;
  char *name = strtok(line, "=");
  assert(name);
  char* value = strtok(NULL, "=");
  if (!value) {
    cerr << "Error:\tParams::Load:\tnull value for param <" << name << "> !" << endl;
    return;
  }
  
  //value = Trim(string(value)).c_str();
  //  name = Trim(string(name)).c_str();
  
  if (value[0] == '\"' || value[0] == ' ')
    value++;
  if (value[strlen(value) - 1] == '\"' || value[strlen(value) - 1] == ' ')
    value[strlen(value)-1] = 0;

  string svalue = string(value);
  svalue.erase(std::remove_if(svalue.begin(), svalue.end(), removeDelimiter()),
               svalue.end());

  string sname = string(name);
  sname.erase(std::remove_if(sname.begin(), sname.end(), removeDelimiter()),
              sname.end());

  cout << sname << " = " << svalue << endl;
  valueMap[sname] = svalue;
}

void Params::load(FILE *file)
{
  char line[256];
  while(fgets(line, 256, file)) {
    parse(line);
  }
}

void Params::load(iostream &io)
{
  char line[256];
  while(io.getline(line, 256)) {
    parse(line);
  }
}


void Params::save(const char *fileName) const
{
  FILE *file = fopen(fileName, "w+");
  if (!file) {
    cout << "[W] Params::Load: failed to save to file " << fileName << endl;
    return;
  }
  save(file);
  fclose(file);
}

void Params::save(FILE *file) const
{
  print(file);
}

void Params::save(iostream &io) const
{
  print(io);
}

void Params::setInt(const char* name, int value) {
  memset(buf, 0, DSL_PARAMS_MBL);
  sprintf(buf, "%d", value);
  valueMap[name] = string(buf);
}

void Params::setFloat(const char* name, float value) {
  memset(buf, 0, DSL_PARAMS_MBL);
  sprintf(buf, "%f", value);
  valueMap[name] = string(buf);
}

void Params::setDouble(const char* name, double value) {
  memset(buf, 0, DSL_PARAMS_MBL);
  sprintf(buf, "%lf", value);
  valueMap[name] = string(buf);
}

void Params::setVector2d(const char* name, const Vector2d& v) {
  vector<double>::const_iterator it;
  stringstream s;
  unsigned int i = 0;
  for (int i = 0; i < v.size(); ++i) {
    s << v[i];
    if (i < v.size()-1)
      s << ",";
  }
  valueMap[name] = s.str();
}

bool Params::getVector2d(const char* name, Vector2d& v) const {
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    // cerr << "Error:\tParams::get:\tparameter <" << name << "> not found!" <<
    // endl;
    return false;
  }

  vector<string> tokens;
  tokenize(i->second, tokens, ", ");
  vector<string>::iterator it;
  int vi = 0;
  for (it = tokens.begin(); it != tokens.end(); ++it) {
    string str = *it;
    replace(str, string("pi"), string("3.141592"));  
    v[vi] = atof(str.c_str());
    ++vi;
  }
  return true;
}

void Params::setVector3d(const char* name, const Vector3d& v) {
  stringstream s;
  unsigned int i = 0;
  for (int i = 0; i < v.size(); ++i) {
    s << v[i];
    if (i < v.size()-1)
      s << ",";
  }
  valueMap[name] = s.str();
}

bool Params::getVector3d(const char* name, Vector3d& v) const {
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    // cerr << "Error:\tParams::get:\tparameter <" << name << "> not found!" <<
    // endl;
    return false;
  }

  vector<string> tokens;
  tokenize(i->second, tokens, ", ");
  vector<string>::iterator it;
  int vi = 0;
  for (it = tokens.begin(); it != tokens.end(); ++it) {
    string str = *it;
    process<double>(str, v[vi]);
    ++vi;
  }
  return true;
}

void Params::setVector4d(const char* name, const Vector4d& v) {
  stringstream s;
  unsigned int i = 0;
  for (int i = 0; i < v.size(); ++i) {
    s << v[i];
    if (i < v.size()-1)
      s << ",";
  }
  valueMap[name] = s.str();
}

bool Params::getVector4d(const char* name, Vector4d& v) const {
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    // cerr << "Error:\tParams::get:\tparameter <" << name << "> not found!" <<
    // endl;
    return false;
  }

  vector<string> tokens;
  tokenize(i->second, tokens, ", ");
  vector<string>::iterator it;
  int vi = 0;
  for (it = tokens.begin(); it != tokens.end(); ++it) {
    string str = *it;
    process<double>(str, v[vi]);
    ++vi;
  }
  return true;
}

void Params::setVector5d(const char* name, const Vector5d& v) {
  stringstream s;
  unsigned int i = 0;
  for (int i = 0; i < v.size(); ++i) {
    s << v[i];
    if (i < v.size()-1)
      s << ",";
  }
  valueMap[name] = s.str();
}

bool Params::getVector5d(const char* name, Vector5d& v) const {
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    // cerr << "Error:\tParams::get:\tparameter <" << name << "> not found!" <<
    // endl;
    return false;
  }

  vector<string> tokens;
  tokenize(i->second, tokens, ", ");
  vector<string>::iterator it;
  int vi = 0;
  for (it = tokens.begin(); it != tokens.end(); ++it) {
    string str = *it;
    process<double>(str, v[vi]);
    ++vi;
  }
  return true;
}

void Params::setVector6d(const char* name, const Vector6d& v) {
  stringstream s;
  unsigned int i = 0;
  for (int i = 0; i < v.size(); ++i) {
    s << v[i];
    if (i < v.size()-1)
      s << ",";
  }
  valueMap[name] = s.str();
}

bool Params::getVector6d(const char* name, Vector6d& v) const {
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    // cerr << "Error:\tParams::get:\tparameter <" << name << "> not found!" <<
    // endl;
    return false;
  }

  vector<string> tokens;
  tokenize(i->second, tokens, ", ");

  if (tokens.size() != 6) {
    cout << "[E] Params::getMatrix6d: expecting 6 doubles, got: "
         << tokens.size() << endl;
    return false;
  }

  vector<string>::iterator it;
  int vi = 0;
  for (it = tokens.begin(); it != tokens.end(); ++it) {
    string str = *it;
    process<double>(str, v[vi]);
    ++vi;
  }
  return true;
}

void Params::setMatrix6d(const char* name, const Matrix6d& m) {
  stringstream s;
  for (int i = 0; i < m.rows(); ++i) {
    for (int j = 0; j < m.cols(); ++j) {
      s << m(i,j);
      if (i < m.rows() - 1 && j < m.cols() - 1)
        s << ",";
    }
  }
  valueMap[name] = s.str();
}

bool Params::getMatrix6d(const char* name, Matrix6d& m) const {
  std::map<string, string>::const_iterator si = valueMap.find(name);
  if (si == valueMap.end()) {
    // cerr << "Error:\tParams::get:\tparameter <" << name << "> not found!" <<
    // endl;
    return false;
  }

  vector<string> tokens;
  tokenize(si->second, tokens, ", ");

  if (tokens.size() != 36) {
    cout << "[E] Params::getMatrix6d: expecting 36 doubles, got: "
         << tokens.size() << endl;
    return false;
  }

  vector<string>::iterator it;
  int i = 0, j = 0;
  for (it = tokens.begin(); it != tokens.end(); ++it) {
    string str = *it;
    process<double>(str, m(i,j));
    j++;
    if (j==6) {
      i++;
      j = 0;
    }
  }
  return true;
}

void Params::setVectorXd(const char* name, const VectorXd& v) {
  stringstream s;
  unsigned int i = 0;
  for (int i = 0; i < v.size(); ++i) {
    s << v[i];
    if (i < v.size()-1)
      s << ",";
  }
  valueMap[name] = s.str();
}

bool Params::getVectorXd(const char* name, VectorXd& v) const {
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    // cerr << "Error:\tParams::get:\tparameter <" << name << "> not found!" <<
    // endl;
    return false;
  }

  vector<string> tokens;
  tokenize(i->second, tokens, ", ");
  vector<string>::iterator it;
  int vi = 0;
  for (it = tokens.begin(); it != tokens.end(); ++it) {
    if (vi >= v.size()) {
      cout << "[E] Params::getVectorXd: mismatched size in element " << name
           << " with expected v.size=" << v.size() << endl;
      return false;
    }
    string str = *it;
    process<double>(str, v[vi]);
    ++vi;
  }
  return true;
}

void Params::setDoubleVec(const char* name, const vector< double >& v) {
  vector<double>::const_iterator it;
  stringstream s;
  unsigned int i = 0;
  for (it = v.begin(); it != v.end(); ++it) {
    s << *it;
    ++i;
    if (i < v.size())
      s << ",";
  }
  valueMap[name] = s.str();
}

void Params::setDoubleArray(const char* name, int n, const double* v) {
  stringstream s;
  for (int i = 0; i < n; ++i) {
    s << v[i];
    if (i < n - 1)
      s << ",";
  }
  valueMap[name] = s.str();
}

void Params::setFloatArray(const char* name, int n, const float* v) {
  stringstream s;
  for (int i = 0; i < n; ++i) {
    s << v[i];
    if (i < n - 1)
      s << ",";
  }
  valueMap[name] = s.str();
}

void Params::setBool(const char* name, bool value) {
  memset(buf, 0, DSL_PARAMS_MBL);
  sprintf(buf, "%ud", value);
  valueMap[name] = string(buf);
}

void Params::setChar(const char* name, char value) {
  memset(buf, 0, DSL_PARAMS_MBL);
  sprintf(buf, "%c", value);
  valueMap[name] = string(buf);
}

void Params::setString(const char* name, const string& v) {
  valueMap[name] = v;
}

bool Params::getInt(const char* name, int& v) const {
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    //    cerr << "Error:\tParams::get:\tparameter <" << name << "> not found!"
    //    << endl;
    return false;
  }
  v = atoi(i->second.c_str());
  return true;
}

bool Params::getFloat(const char* name, float& v) const {
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    //    cerr << "Error:\tParams::get:\tparameter <" << name << "> not found!"
    //    << endl;
    return false;
  }
  string str = i->second;
  process<float>(str, v);
  return true;
}

bool Params::getDouble(const char* name, double& v) const {
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    // cerr << "Error:\tParams::get:\tparameter <" << name << "> not found!" <<
    // endl;
    return false;
  }
  string str = i->second;
  process<double>(str, v);
  return true;
}

bool Params::getDoubleVec(const char* name, vector< double >& vs) const {
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    // cerr << "Error:\tParams::get:\tparameter <" << name << "> not found!" <<
    // endl;
    return false;
  }

  vector<string> tokens;
  tokenize(i->second, tokens, ",");
  vector<string>::iterator it;
  for (it = tokens.begin(); it != tokens.end(); ++it) {
    string str = *it;
    double v;
    process<double>(str, v);
    vs.push_back(v);
  }
  return true;
}

bool Params::getDoubleArray(const char* name, int n, double* vs) const {
  std::map<string, string>::const_iterator mi = valueMap.find(name);
  if (mi == valueMap.end()) {
    // cerr << "Error:\tParams::get:\tparameter <" << name << "> not found!" <<
    // endl;
    return false;
  }

  vector<string> tokens;
  tokenize(mi->second, tokens, ",");
  vector<string>::iterator it;
  int i = 0;
  for (it = tokens.begin(); it != tokens.end(), i < n; ++it, ++i) {    
    string str = *it;
    assert(i < n);
    process<double>(str, vs[i]);
    ++i;
  }
  return true;
}

bool Params::getFloatArray(const char* name, int n, float* vs) const {
  std::map<string, string>::const_iterator mi = valueMap.find(name);
  if (mi == valueMap.end()) {
    // cerr << "Error:\tParams::get:\tparameter <" << name << "> not found!" <<
    // endl;
    return false;
  }

  vector<string> tokens;
  tokenize(mi->second, tokens, ",");
  vector<string>::iterator it;
  int i = 0;
  for (it = tokens.begin(); it != tokens.end(), i < n; ++it, ++i) {    
    string str = *it;
    assert(i < n);
    process<float>(str, vs[i]);
    ++i;
  }
  return true;
}

bool Params::getBool(const char* name, bool& v) const {
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    //    cerr << "Error:\tParams::get:\tparameter <" << name << "> not found!"
    //    << endl;
    return false;
  }
  v = atoi(i->second.c_str());
  return true;
}

char Params::getChar(const char* name, char& v) const {
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    //    cerr << "Error:\tParams::get:\tparameter <" << name << "> not found!"
    //    << endl;
    return false;
  }
  v = i->second[0];
  return true;
}

bool Params::getString(const char* name, string& v) const {
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    //    cerr << "Error:\tParams::get:\tparameter <" << name << "> not found!"
    //    << endl;
    return false;
  }
  v = i->second;
  return true;
}

void Params::print(FILE* file) const {
  fprintf(file, "# Params: count=%d\n\n", (int)valueMap.size());
  std::map<string, string>::const_iterator i;
  for (i = valueMap.begin(); i != valueMap.end(); ++i)
    fprintf(file, "%s=%s\n", i->first.c_str(), i->second.c_str());
}

void Params::print(iostream& io) const {
  io << "# Params: count=" << valueMap.size() << "\n" << endl;
  std::map<string, string>::const_iterator i;
  for (i = valueMap.begin(); i != valueMap.end(); ++i)
    io << i->first << " " << i->second << endl;
}

}
