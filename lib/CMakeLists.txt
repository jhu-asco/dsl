add_library(dsl SHARED fibheap.c 
                       search.cc
                       utils.cc
                       utilsimg.cc
                       gridutils.cc
                       grid2d.cc
                       grid3d.cc
                       cargrid.cc
                       grid2dconnectivity.cc
                       grid3dconnectivity.cc
                       carconnectivity.cc
                       carconntwist.cc
                       carprimcfg.cc
                       carcost.cc
                       params.cc
)

SET_TARGET_PROPERTIES(dsl
                      PROPERTIES
                      VERSION 1.1.0
                      SOVERSION 1)

INSTALL(TARGETS dsl
        DESTINATION lib)

INSTALL(FILES fibheap.h
              libiberty.h
              ansidecl.h
              vertex.h
              edge.h
              graph.h
              cost.h
              search.h
              utils.h
              utilsimg.h
              gridutils.h
              cell.h
              grid.h
              gridcost.h
              gridconnectivity.h
              gridpath.h
              lineconnectivity.h
              carconnectivity.h
              carconntwist.h
              carprimcfg.h
              carcost.h
              cargrid.h
              carcost.h
              grid2d.h
              grid2dconnectivity.h
              grid3d.h
              grid3dconnectivity.h
              gridsearch.h
              spline.h
              params.h
              DESTINATION include/dsl)              
