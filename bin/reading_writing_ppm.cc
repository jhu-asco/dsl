#include "ppm_reader.h"
#include <iostream>
#include <assert.h>

using namespace dsl;
using namespace std;

int main(int argc, char** argv)
{
  if (argc!=2) {
    cout << "Usage: $./reading_writing_ppm image.ppm" << endl;
    cout << "\t\t where the image.ppm has a P6 signature and a bitdepth of 1 or 255 or 65535" << endl;
    return 0;
  }
  assert(argc == 2);

  cout<<"loaded the input image file to img1"<<endl;
  ImageRGB img1;
  LoadPpm(img1,argv[1]);

  cout<<"converting img1 to bitdepth of 255 and saving it to bitdepth_0xff.ppm"<<endl;
  img1.ChangeBitDepth(ImageRGB::BD8);
  SavePpm(img1,"bitdepth_0xff.ppm");

  cout<<"Loading the file bitdepth_255.ppm as img2"<<endl;
  ImageRGB img2;
  LoadPpm(img2,"bitdepth_0xff.ppm");

  cout<<"converting img2 to bitdepth of 65535 and saving it to bitdepth_0xffff.ppm"<<endl;
  img2.ChangeBitDepth(ImageRGB::BD16);
  SavePpm(img2, "saved_0xffff.ppm");
}
