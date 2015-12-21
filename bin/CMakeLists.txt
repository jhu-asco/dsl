link_directories(${DSL_BINARY_DIR}/lib) 


add_executable(test2d test2d.cc)
#add_executable(test test.cc)
add_executable(test3d test3d.cc)
add_executable(cartest cartest.cc)

target_link_libraries(test2d dsl)
target_link_libraries(test3d dsl)
target_link_libraries(cartest dsl)

INSTALL(TARGETS test2d
        RUNTIME
        DESTINATION bin)

INSTALL(FILES test.cc map.ppm 
        DESTINATION bin)
