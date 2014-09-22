include(LibFindMacros)

# Dependencies
#libfind_package(itpp)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(OPENCV_PKGCONF itpp)

# Include dir
find_path(OPENCV_INCLUDE_DIR
  NAMES cv.h
  PATHS ${OPENCV_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(OPENCV_LIBRARY
  NAMES cv
  PATHS ${OPENCV_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(OPENCV_PROCESS_INCLUDES OPENCV_INCLUDE_DIR)
set(OPENCV_PROCESS_LIBS OPENCV_LIBRARY)
libfind_process(OPENCV)