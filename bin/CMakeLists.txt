link_directories(${DSL_BINARY_DIR}/lib) 

add_executable(test test.cc)

target_link_libraries(test dsl)

INSTALL(TARGETS test
        RUNTIME
        DESTINATION bin)

INSTALL(FILES test.cc map.ppm 
        DESTINATION bin)