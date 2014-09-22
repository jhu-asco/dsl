include(LibFindMacros)

# Dependencies
#libfind_package(itpp)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(ITPP_PKGCONF itpp)

# Include dir
find_path(ITPP_INCLUDE_DIR
  NAMES itbase.h
  PATHS ${ITPP_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(ITPP_LIBRARY
  NAMES itpp
  PATHS ${ITPP_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(ITPP_PROCESS_INCLUDES ITPP_INCLUDE_DIR)
set(ITPP_PROCESS_LIBS ITPP_LIBRARY)
libfind_process(ITPP)