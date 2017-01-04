link_directories(${DSL_BINARY_DIR}/lib) 

add_executable(test2d test2d.cc)
target_link_libraries(test2d dsl)

add_executable(cartest cartest.cc)
target_link_libraries(cartest dsl)

add_executable(terraincartest terraincartest.cc)
target_link_libraries(terraincartest dsl)

add_executable(test_cargrid test_cargrid.cc)
target_link_libraries(test_cargrid dsl)

# Testing some function of gridcore class
add_executable(gridcoretest gridcoretest.cc)
target_link_libraries(gridcoretest dsl)

# Testing id to index conversion time for different methods
add_executable(idtoindexspeedtest idtoindexspeedtest.cc)
target_link_libraries(idtoindexspeedtest dsl)

# Testing speed for different methods to loop over all elements 
add_executable(loopindexspeedtest loopindexspeedtest.cc)
target_link_libraries(loopindexspeedtest dsl)

# Testing the function of gridutils class
add_executable(gridutilstest gridutilstest.cc)
target_link_libraries(gridutilstest dsl)

# Testing reading writing of .ppm files
add_executable(reading_writing_ppm reading_writing_ppm.cc)
target_link_libraries(reading_writing_ppm dsl)

# Testing yaml file parsing
#add_executable(yamlreadtest yamlreadtest.cc)
#target_link_libraries(yamlreadtest dsl)

#add_executable(test test.cc)

#add_executable(test3d test3d.cc)
#target_link_libraries(test3d dsl)


#add_executable(carwithgeomtest carwithgeomtest.cc)
#target_link_libraries(carwithgeomtest dsl ${OpenCV_LIBRARIES})


#INSTALL(TARGETS test2d
#        RUNTIME
#        DESTINATION bin)
#
#INSTALL(FILES map.ppm 
#        DESTINATION bin)
