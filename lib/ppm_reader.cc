#include "ppm_reader.h"
#include <exception>
#include <iostream>

namespace dsl{
using namespace std;

void ImageRGB::ChangeBitDepth(BitDepth bitdepth_new){
  switch(bitdepth){
    case BitDepth::BD8:
      if(bitdepth_new == BitDepth::BD8){
        return;
      }else{
        for(int id = 0; id < w*h; id++){
          rdata[id] = rdata[id]<<8;
          gdata[id] = gdata[id]<<8;
          bdata[id] = bdata[id]<<8;
        }
      }
      break;
    case BitDepth::BD16:
      if(bitdepth_new == BitDepth::BD16){
        return;
      }else{
        for(int id = 0; id < w*h; id++){
          rdata[id] = rdata[id]>>8;
          gdata[id] = gdata[id]>>8;
          bdata[id] = bdata[id]>>8;
        }
      }
      break;
  }
  bitdepth = bitdepth_new;
}

void ImageRGB::resize(int size){
  rdata.resize(size);
  gdata.resize(size);
  bdata.resize(size);
}

void ImageRGB::set_to_white(int id, uint16_t val){
  set_color(id, val, WHITE);
}

void ImageRGB::set_to_red(int id, uint16_t val){
  set_color(id, val, RED);
}

void ImageRGB::set_to_blue(int id, uint16_t val){
  set_color(id, val, BLUE);
}

void ImageRGB::set_to_green(int id, uint16_t val){
  set_color(id, val, GREEN);
}

void ImageRGB::set_to_yellow(int id, uint16_t val){
  set_color(id, val, YELLOW);
}

void ImageRGB::set_to_cyan(int id, uint16_t val){
  set_color(id, val, CYAN);
}

void ImageRGB::set_to_magenta(int id, uint16_t val){
  set_color(id, val, MAGENTA);
}

void ImageRGB::set_color(int id, uint16_t val, Color color){
switch (color) {
  case RED:
    rdata[id] = val;
    gdata[id] = 0;
    bdata[id] = 0;
    break;
  case GREEN:
    rdata[id] = 0;
    gdata[id] = val;
    bdata[id] = 0;
    break;
  case BLUE:
    rdata[id] = 0;
    gdata[id] = 0;
    bdata[id] = val;
    break;
  case WHITE:
    rdata[id] = val;
    gdata[id] = val;
    bdata[id] = val;
    break;
  case YELLOW:
    rdata[id] = val;
    gdata[id] = val;
    bdata[id] = 0;
    break;
  case CYAN:
    rdata[id] = 0;
    gdata[id] = val;
    bdata[id] = val;
    break;
  case MAGENTA:
    rdata[id] = val;
    gdata[id] = 0;
    bdata[id] = val;
    break;
  default:
    rdata[id] = val;
    gdata[id] = val;
    bdata[id] = val;
    break;
}

}

bool LoadPpm(ImageRGB &img, const string &filename){
  if(filename.compare(filename.size() - 4, 4, ".ppm")){
    cout<<"File doesn't have .ppm extension"<<endl;
    return false;
  }

  ifstream f(filename.c_str(), ios::binary);
  if (f.fail()){
    cout << "Could not open file: " << filename << endl;
    return false;
  }

  // get type of file
  RemoveComment(f);
  string header;
  f >> header;

  // get w
  RemoveComment(f);
  f >> img.w;

  // get h
  RemoveComment(f);
  f >> img.h;

  // get bitdepth
  RemoveComment(f);
  long bitdepth = 0;
  f >> bitdepth;
  img.bitdepth = bitdepth;

  // error checking
  if (header != "P6"){
    cout << "Unsupported magic number" << endl;
    f.close();
    return false;
  }

  if (img.w < 1){
    cout << "Unsupported width: " << img.w << endl;
    f.close();
    return false;
  }

  if (img.h < 1){
    cout << "Unsupported height: " << img.h << endl;
    f.close();
    return false;
  }

  if (bitdepth!=0xff && bitdepth != 0xffff){
    cout << "Unsupported bitdepth: " << bitdepth << endl;
    f.close();
    return false;
  }

  // load image data
  img.rdata.resize(img.w * img.h);
  img.gdata.resize(img.w * img.h);
  img.bdata.resize(img.w * img.h);

  f.get();

  if(bitdepth == 0xff){
    vector<uint8_t> data(img.w * img.h * 3);
    f.read((char*)&data[0], data.size());
    for(int i=0; i< img.w * img.h ; i++){
      int id = i*3;
      img.rdata[i] = data[id];
      img.gdata[i] = data[id+1];
      img.bdata[i] = data[id+2];
    }
  }

  if(bitdepth == 0xffff){
    vector<uint8_t> data(img.w * img.h * 3 *2); //3 channel and 2 byte per channel
    f.read((char*)&data[0], data.size());
    uint16_t temp;
    for(int i=0; i< img.w * img.h ; i++){
      int id = i*6;
      //Most Significant Byte fist according to Netpbm standard
      temp = data[id];
      img.rdata[i] = (temp<<8) | data[id+1];
      temp = data[id+2];
      img.gdata[i] = (temp<<8) | data[id+3];
      temp = data[id+4];
      img.bdata[i] = (temp<<8) | data[id+5];
    }

  }

  // close file
  f.close();

  return true;
}

bool SavePpm(ImageRGB& img,  const string& filename) {

  if(filename.compare(filename.size() - 4, 4, ".ppm")){
    cout<<"File doesn't have .ppm extension"<<endl;
    return false;
  }

  std::fstream fs(filename, std::fstream::out);
  if(!fs.is_open()){
    cout<<"couldn't open the image file to write"<<endl;
    return false;
  }

  if(img.bitdepth == 0xff){
    fs << "P6" << std::endl << img.w << " " << img.h << std::endl << "255" << std::endl;
    uint8_t data[img.h*img.w*3];
    for(int i=0; i< img.w * img.h ; i++){
      int id = i*3;
      data[id] = (uint8_t)img.rdata[i];
      data[id + 1] = (uint8_t)img.gdata[i];
      data[id + 2] = (uint8_t)img.bdata[i];
    }
    fs.write((char*)data, img.h*img.w*3);
  }else if(img.bitdepth == 0xffff ){
    fs << "P6" << std::endl << img.w << " " << img.h << std::endl << "65535" << std::endl;
    uint8_t data[img.h*img.w*3*2];
    for(int i=0; i< img.w * img.h ; i++){
      int id = i*6;
      //Most Significant Byte fist according to Netpbm standard
      data[id]   = (uint8_t)(img.rdata[i]>>8);
      data[id+1] = (uint8_t)img.rdata[i];
      data[id+2] = (uint8_t)(img.gdata[i]>>8);
      data[id+3] = (uint8_t)img.gdata[i];
      data[id+4] = (uint8_t)(img.bdata[i]>>8);
      data[id+5] = (uint8_t)img.bdata[i];
    }
    fs.write((char*)data, img.h*img.w*6);
  }else{
    cout<<"The bitdepth of the image is not supported"<<endl;
    return false;
  }

  fs.close();
  return true;
}


void RemoveComment(ifstream &f){
  char linebuf[1024];
  char ppp;
  while (ppp = f.peek(), ppp == '\n' || ppp == '\r')
    f.get();
  if (ppp == '#')
    f.getline(linebuf, 1023);
}

}

