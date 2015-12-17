#ifndef DSL_UTILS_H
#define DSL_UTILS_H

#include <sys/time.h>

namespace dsl {
  
  void save_map(const char* map, int width, int height, const char* filename);
  
  char* load_map(int* width, int* height, const char* filename);
  
  void timer_start(struct timeval *time);
  
  long timer_us(struct timeval *time) ;
  
}

#endif
