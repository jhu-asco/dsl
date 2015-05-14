cmake_minimum_required (VERSION 2.8) 
project (dsl)

include_directories(lib)

ADD_DEFINITIONS(-g -O2 -Wall -W)

add_subdirectory (lib) 
add_subdirectory (bin) 
#add_subdirectory (docs) 

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Dynamic D* Lite Graph Search Library")
SET(CPACK_PACKAGE_VENDOR "Marin Kobilarov")
#SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/ReadMe.txt")
#SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "0")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")
#SET(CPACK_STRIP_FILES "bin/MyExecutable")
#SET(CPACK_SOURCE_STRIP_FILES "")
SET(CPACK_PACKAGE_EXECUTABLES test "Test search")


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

INSTALL(FILES dsl.pc
        DESTINATION lib/pkgconfig)

set(CPACK_SOURCE_IGNORE_FILES "/build/;/CVS/;/\\\\.svn/;\\\\.swp$;~$;\\\\.\\\\#;/\\\\#")

INCLUDE(CPack)
