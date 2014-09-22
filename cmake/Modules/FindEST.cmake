include(LibFindMacros)

# Dependencies
#libfind_package(est)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(EST_PKGCONF est)

# Include dir
find_path(EST_INCLUDE_DIR
  NAMES search.h
  PATHS ${EST_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(EST_LIBRARY
  NAMES est
  PATHS ${EST_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(EST_PROCESS_INCLUDES EST_INCLUDE_DIR)
set(EST_PROCESS_LIBS EST_LIBRARY)
libfind_process(EST)