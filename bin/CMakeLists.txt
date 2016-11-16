link_directories(${DSL_BINARY_DIR}/lib) 

add_executable(test2d test2d.cc)
target_link_libraries(test2d dsl)

add_executable(cartest cartest.cc)
target_link_libraries(cartest dsl)

add_executable(cartesttwistprim cartesttwistprim.cc)
target_link_libraries(cartesttwistprim dsl)

add_executable(test_cargrid test_cargrid.cc)
target_link_libraries(test_cargrid dsl)

# Testing the function of gridcore class
add_executable(gridcoretest gridcoretest.cc)
target_link_libraries(gridcoretest dsl)

#add_executable(test test.cc)

#add_executable(test3d test3d.cc)
#target_link_libraries(test3d dsl)


#add_executable(carwithgeomtest carwithgeomtest.cc)
#target_link_libraries(carwithgeomtest dsl ${OpenCV_LIBRARIES})


INSTALL(TARGETS test2d
        RUNTIME
        DESTINATION bin)

INSTALL(FILES map.ppm 
        DESTINATION bin)
