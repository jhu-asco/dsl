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


/*
static void Trim(string& str)
{
  string::size_type pos = str.find_last_not_of(' ');
  if(pos != string::npos) {
    str.erase(pos + 1);
    pos = str.find_first_not_of(' ');
    if(pos != string::npos) str.erase(0, pos);
  }
  else str.erase(str.begin(), str.end());
}
*/

static void Tokenize(const string& str,
                     vector<string>& tokens,
                     const string& delimiters = " ")
{
  // Skip delimiters at beginning.
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  string::size_type pos     = str.find_first_of(delimiters, lastPos);
  
  while (string::npos != pos || string::npos != lastPos)
    {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of(delimiters, pos);
      // Find next "non-delimiter"
      pos = str.find_first_of(delimiters, lastPos);
    }
}




Params::Params()
{
}

Params::Params(const char *fileName)
{
  Load(fileName);
}

Params::Params(FILE *file)
{
  Load(file);
}

Params::Params(iostream &io)
{
  Load(io);
}


Params::~Params()
{
  valueMap.clear();
}

void Params::Load(const char *fileName)
{
  FILE *file = fopen(fileName, "r");
  if (!file) {
    cout << "[W] Params::Load: failed to load file " << fileName << endl;
    return;
  }
  Load(file);
  fclose(file);
}

void Params::Parse(char *line)
{
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
  svalue.erase( std::remove_if( svalue.begin(), svalue.end(), RemoveDelimiter()), svalue.end()); 

  string sname = string(name);
  sname.erase( std::remove_if( sname.begin(), sname.end(), RemoveDelimiter()), sname.end()); 

  cout << "  "<< sname << " = " << svalue << endl;
  valueMap[sname] = svalue;
}

void Params::Load(FILE *file)
{
  char line[256];
  cout<<"Params read:"<<endl;
  while(fgets(line, 256, file)) {
    Parse(line);
  }
  cout<<endl;
}

void Params::Load(iostream &io)
{
  char line[256];
  cout<<"Params read:"<<endl;
  while(io.getline(line, 256)) {
    Parse(line);        
  }
  cout<<endl;
}


void Params::Save(const char *fileName) const
{
  FILE *file = fopen(fileName, "w+");
  if (!file) {
    cout << "[W] Params::Load: failed to save to file " << fileName << endl;
    return;
  }
  Save(file);
  fclose(file);
}

void Params::Save(FILE *file) const
{
  Print(file);
}

void Params::Save(iostream &io) const
{
  Print(io);
}

void Params::SetInt(const char *name, int value)
{
  memset(buf, 0, DSL_PARAMS_MBL);
  sprintf(buf, "%d", value);
  valueMap[name] = string(buf);
}

void Params::SetFloat(const char *name, float value)
{
  memset(buf, 0, DSL_PARAMS_MBL);
  sprintf(buf, "%f", value);
  valueMap[name] = string(buf);
}

void Params::SetDouble(const char *name, double value)
{
  memset(buf, 0, DSL_PARAMS_MBL);
  sprintf(buf, "%lf", value);
  valueMap[name] = string(buf);
}


void Params::SetVector2d(const char *name, const Vector2d &v)
{
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

bool Params::GetVector2d(const char *name, Vector2d &v) const
{
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    // cerr << "Error:\tParams::Get:\tparameter <" << name << "> not found!" << endl;
    return false;
  }

  vector<string> tokens;
  Tokenize(i->second, tokens, ", ");
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

void Params::SetVector3d(const char *name, const Vector3d &v)
{
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

bool Params::GetVector3d(const char *name, Vector3d &v) const
{
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    // cerr << "Error:\tParams::Get:\tparameter <" << name << "> not found!" << endl;
    return false;
  }

  vector<string> tokens;
  Tokenize(i->second, tokens, ", ");
  vector<string>::iterator it;
  int vi = 0;
  for (it = tokens.begin(); it != tokens.end(); ++it) {
    string str = *it;
    process<double>(str, v[vi]);
    ++vi;
  }
  return true;
}

void Params::SetVector4d(const char *name, const Vector4d &v)
{
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


bool Params::GetVector4d(const char *name, Vector4d &v) const
{
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    // cerr << "Error:\tParams::Get:\tparameter <" << name << "> not found!" << endl;
    return false;
  }

  vector<string> tokens;
  Tokenize(i->second, tokens, ", ");
  vector<string>::iterator it;
  int vi = 0;
  for (it = tokens.begin(); it != tokens.end(); ++it) {
    string str = *it;
    process<double>(str, v[vi]);
    ++vi;
  }
  return true;
}

void Params::SetVector5d(const char *name, const Vector5d &v)
{
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


bool Params::GetVector5d(const char *name, Vector5d &v) const
{
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    // cerr << "Error:\tParams::Get:\tparameter <" << name << "> not found!" << endl;
    return false;
  }

  vector<string> tokens;
  Tokenize(i->second, tokens, ", ");
  vector<string>::iterator it;
  int vi = 0;
  for (it = tokens.begin(); it != tokens.end(); ++it) {
    string str = *it;
    process<double>(str, v[vi]);
    ++vi;
  }
  return true;
}


void Params::SetVector6d(const char *name, const Vector6d &v)
{
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


bool Params::GetVector6d(const char *name, Vector6d &v) const
{
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    // cerr << "Error:\tParams::Get:\tparameter <" << name << "> not found!" << endl;
    return false;
  }

  vector<string> tokens;
  Tokenize(i->second, tokens, ", ");

  if (tokens.size() != 6) {
    cout << "[E] Params::GetMatrix6d: expecting 6 doubles, got: " << tokens.size() << endl;
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


void Params::SetMatrix6d(const char *name, const Matrix6d &m)
{
  vector<double>::const_iterator it;
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


bool Params::GetMatrix6d(const char *name, Matrix6d &m) const
{
  std::map<string, string>::const_iterator si = valueMap.find(name);
  if (si == valueMap.end()) {
    // cerr << "Error:\tParams::Get:\tparameter <" << name << "> not found!" << endl;
    return false;
  }

  vector<string> tokens;
  Tokenize(si->second, tokens, ", ");

  if (tokens.size() != 36) {
    cout << "[E] Params::GetMatrix6d: expecting 36 doubles, got: " << tokens.size() << endl;
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


void Params::SetVectorXd(const char *name, const VectorXd &v)
{
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

bool Params::GetVectorXd(const char *name, VectorXd &v) const
{
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    // cerr << "Error:\tParams::Get:\tparameter <" << name << "> not found!" << endl;
    return false;
  }

  vector<string> tokens;
  Tokenize(i->second, tokens, ", ");
  vector<string>::iterator it;
  int vi = 0;
  for (it = tokens.begin(); it != tokens.end(); ++it) {
    if (vi >= v.size()) {
      cout << "[E] Params::GetVectorXd: mismatched size in element " << name << " with expected v.size=" << v.size() << endl;
      return false;
    }
    string str = *it;
    process<double>(str, v[vi]);
    ++vi;
  }
  return true;
}



void Params::SetDoubleVec(const char *name, const vector<double> &v)
{
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

void Params::SetDoubleArray(const char *name, int n, const double *v)
{
  stringstream s;
  for (int i = 0; i < n; ++i) {
    s << v[i];
    if (i < n - 1)
      s << ",";
  }
  valueMap[name] = s.str();
}

void Params::SetFloatArray(const char *name, int n, const float *v)
{
  stringstream s;
  for (int i = 0; i < n; ++i) {
    s << v[i];
    if (i < n - 1)
      s << ",";
  }
  valueMap[name] = s.str();
}



void Params::SetBool(const char *name, bool value)
{
  memset(buf, 0, DSL_PARAMS_MBL);
  sprintf(buf, "%ud", value);
  valueMap[name] = string(buf);
}

void Params::SetChar(const char *name, char value)
{
  memset(buf, 0, DSL_PARAMS_MBL);
  sprintf(buf, "%c", value);
  valueMap[name] = string(buf);
}


void Params::SetString(const char *name, const string &v)
{
  valueMap[name] = v;
}

bool Params::GetInt(const char *name, int &v) const
{
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    //    cerr << "Error:\tParams::Get:\tparameter <" << name << "> not found!" << endl;
    return false;
  }
  v = atoi(i->second.c_str());
  return true;
}


bool Params::GetFloat(const char *name, float &v) const
{
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    //    cerr << "Error:\tParams::Get:\tparameter <" << name << "> not found!" << endl;
    return false;
  }
  string str = i->second;
  process<float>(str, v);
  return true;
}


bool Params::GetDouble(const char *name, double &v) const
{
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    // cerr << "Error:\tParams::Get:\tparameter <" << name << "> not found!" << endl;
    return false;
  }
  string str = i->second;
  process<double>(str, v);
  return true;
}



bool Params::GetDoubleVec(const char *name, vector<double> &vs) const
{
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    // cerr << "Error:\tParams::Get:\tparameter <" << name << "> not found!" << endl;
    return false;
  }

  vector<string> tokens;
  Tokenize(i->second, tokens, ",");
  vector<string>::iterator it;
  for (it = tokens.begin(); it != tokens.end(); ++it) {
    string str = *it;
    double v;
    process<double>(str, v);
    vs.push_back(v);
  }
  return true;
}


bool Params::GetDoubleArray(const char *name, int n, double *vs) const
{
  std::map<string, string>::const_iterator mi = valueMap.find(name);
  if (mi == valueMap.end()) {
    // cerr << "Error:\tParams::Get:\tparameter <" << name << "> not found!" << endl;
    return false;
  }

  vector<string> tokens;
  Tokenize(mi->second, tokens, ",");
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

bool Params::GetFloatArray(const char *name, int n, float *vs) const
{
  std::map<string, string>::const_iterator mi = valueMap.find(name);
  if (mi == valueMap.end()) {
    // cerr << "Error:\tParams::Get:\tparameter <" << name << "> not found!" << endl;
    return false;
  }

  vector<string> tokens;
  Tokenize(mi->second, tokens, ",");
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


bool Params::GetBool(const char *name, bool &v) const
{
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    //    cerr << "Error:\tParams::Get:\tparameter <" << name << "> not found!" << endl;
    return false;
  }
  v = atoi(i->second.c_str());
  return true;
}

char Params::GetChar(const char *name, char &v) const
{
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    //    cerr << "Error:\tParams::Get:\tparameter <" << name << "> not found!" << endl;
    return false;
  }
  v = i->second[0];
  return true;
}


bool Params::GetString(const char *name, string &v) const
{
  std::map<string, string>::const_iterator i = valueMap.find(name);
  if (i == valueMap.end()) {
    //    cerr << "Error:\tParams::Get:\tparameter <" << name << "> not found!" << endl;
    return false;
  }
  v = i->second;
  return true;
}


void Params::Print(FILE *file) const
{
  fprintf(file, "# Params: count=%d\n\n", (int)valueMap.size());
  std::map<string, string>::const_iterator i;
  for (i = valueMap.begin(); i != valueMap.end(); ++i)
    fprintf(file, "%s=%s\n", i->first.c_str(), i->second.c_str());
}

void Params::Print(iostream &io) const
{
  io << "# Params: count=" << valueMap.size() << "\n" << endl;
  std::map<string, string>::const_iterator i;
  for (i = valueMap.begin(); i != valueMap.end(); ++i)
    io << i->first << " " << i->second << endl;
}

}
