# Defines
#
#  ITPP_FOUND - the system has IT++
#  ITPP_INCLUDE_DIR - the include directory
#  ITPP_LINK_DIR - the link directory
#  ITPP_LIBRARIES - the libraries needed to use IT++

# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.

if(ITPP_INCLUDE_DIR AND ITPP_LIBRARIES)
   # in cache already
   set(Itpp_FIND_QUIETLY TRUE)
endif(ITPP_INCLUDE_DIR AND ITPP_LIBRARIES)

if(MSVC)
   SET(ITPP_INCLUDE_DIR "C:/Program Files/itpp")
   SET(ITPP_LINK_DIR "C:/Program Files/itpp/lib" $ENV{LIB})
   SET(ITPP_LIBRARIES mkl_c_dll) # itpp)
ELSE (MSVC)
   # use pkg-config to get the directories and then use these values
   # in the FIND_PATH() and FIND_LIBRARY() calls
   INCLUDE(UsePkgConfig)
   PKGCONFIG(itpp ITPP_INCLUDE_DIR ITPP_LINK_DIR _ItppLinkFlags _ItppCflags)
   SET(ITPP_LIBRARIES) #itpp)
   # we should somehow distinguish between release and debug
ENDIF (MSVC)

if(ITPP_INCLUDE_DIR AND ITPP_LIBRARIES)
   set(ITPP_FOUND TRUE)
else(ITPP_INCLUDE_DIR AND ITPP_LIBRARIES)
   set(ITPP_FOUND FALSE)
endif(ITPP_INCLUDE_DIR AND ITPP_LIBRARIES)

if (ITPP_FOUND)
   if (NOT Itpp_FIND_QUIETLY)
      message(STATUS "Found IT++: ${ITPP_LIBRARIES}")
   endif (NOT Itpp_FIND_QUIETLY)
else (ITPP_FOUND)
   if (Itpp_FIND_REQUIRED)
      message(SEND_ERROR "Could NOT find IT++")
   endif (Itpp_FIND_REQUIRED)
endif (ITPP_FOUND)

mark_as_advanced(ITPP_INCLUDE_DIR ITPP_LIBRARIES)