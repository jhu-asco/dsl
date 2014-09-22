find_package(PkgConfig)
pkg_check_modules(PC_LIBITPP QUIETLY itpp)
set(LIBITPP_DEFINITIONS ${PC_LIBITPP_CFLAGS_OTHER})

find_path(LIBITPP_INCLUDE_DIR itpp/itbase.h
          HINTS ${PC_LIBITPP_INCLUDEDIR} ${PC_LIBITPP_INCLUDE_DIRS}
          PATH_SUFFIXES libitpp )

find_library(LIBITPP_LIBRARY NAMES itpp libitpp
             HINTS ${PC_LIBITPP_LIBDIR} ${PC_LIBITPP_LIBRARY_DIRS} )

set(LIBITPP_LIBRARIES ${LIBITPP_LIBRARY} )
set(LIBITPP_INCLUDE_DIRS ${LIBITPP_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBXML2_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(LibITPP DEFAULT_MSG
                                  LIBITPP_LIBRARY LIBITPP_INCLUDE_DIR)

mark_as_advanced(LIBITPP_INCLUDE_DIR LIBITPP_LIBRARY )
