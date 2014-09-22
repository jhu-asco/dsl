include(LibFindMacros)

# Dependencies
#libfind_package(dsl)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(SNOPT_PKGCONF snopt)

# Include dir
find_path(SNOPT_INCLUDE_DIR
  NAMES snopt/snopt.hh
  PATHS ${SNOPT_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(SNOPT_LIBRARY
  NAMES snopt
  PATHS ${SNOPT_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(SNOPT_PROCESS_INCLUDES SNOPT_INCLUDE_DIR)
set(SNOPT_PROCESS_LIBS SNOPT_LIBRARY)
libfind_process(SNOPT)