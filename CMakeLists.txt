cmake_minimum_required (VERSION 2.8) 
project (dsl)
#set(CMAKE_BUILD_TYPE Debug)
include_directories(lib)

SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-std=c++11 -pthread")
ADD_DEFINITIONS(-g -O -Wall -Werror -W -Wno-unused-parameter -Wno-unused-variable)

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Dynamic D* Lite Graph Search Library")
SET(CPACK_PACKAGE_VENDOR "Marin Kobilarov")
#SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/ReadMe.txt")
#SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "1")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")
#SET(CPACK_STRIP_FILES "bin/MyExecutable")
#SET(CPACK_SOURCE_STRIP_FILES "")
SET(CPACK_PACKAGE_EXECUTABLES test2d "2D grid test search")

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package( Protobuf REQUIRED )
include_directories(${PROTOBUF_INCLUDE_DIRS})
set(LIBS ${LIBS} ${PROTOBUF_LIBRARIES})

include_directories(third_party)

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
  )
endif(DOXYGEN_FOUND)

include_directories( ${CMAKE_BINARY_DIR}/lib )

add_subdirectory (lib) 
add_subdirectory (bin) 

#add_subdirectory (docs)

INSTALL(FILES dsl.pc
        DESTINATION lib/pkgconfig)

set(CPACK_SOURCE_IGNORE_FILES "/.git/;/build/;/CVS/;/\\\\.svn/;\\\\.swp$;~$;\\\\.\\\\#;/\\\\#")

INCLUDE(CPack)
