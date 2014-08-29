link_directories(${DSL_BINARY_DIR}/lib) 

add_executable(test test.cc)
add_executable(test3D test3D.cc)

target_link_libraries(test dsl)
target_link_libraries(test3D dsl)

INSTALL(TARGETS test
        RUNTIME
        DESTINATION bin)

INSTALL(FILES test.cc map.ppm 
        DESTINATION bin)
