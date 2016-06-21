add_library(dsl SHARED fibheap.c 
                       search.cc
                       utils.cc
                       utilsimg.cc
                       grid2d.cc
                       grid3d.cc
                       cargrid.cc
                       grid2dconnectivity.cc
                       grid3dconnectivity.cc
                       carprimcfg.cc
                       carconnectivity.cc
                       carconnwithpath.cc
                       carconnwithprim.cc
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
              cell.h
              grid.h
              gridcost.h
              gridconnectivity.h
              gridpath.h
              lineconnectivity.h
              carprimcfg.h
              carconnectivity.h
              carconnwithpath.h
              carconnwithprim.h              
              carcost.h
              cargrid.h
              traversabilitycost.h
              grid2d.h
              grid2dconnectivity.h
              grid3d.h
              grid3dconnectivity.h
              gridsearch.h
              spline.h
              params.h
              map.h
              DESTINATION include/dsl)              
