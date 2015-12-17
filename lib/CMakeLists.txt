add_library(dsl SHARED fibheap.c 
                       search.cc
                       utils.cc
                       grid2d.cc
                       cargrid.cc
                       grid2dconnectivity.cc
                       carconnectivity.cc
                       gridcost3d.cc
                       gridtravcost.cc
                       gridsearch3d.cc
                       carcost.cc
)


SET_TARGET_PROPERTIES(dsl
                      PROPERTIES
                      VERSION 1.0.0
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
              cell.h
              grid.h
              gridcost.h
              gridconnectivity.h
              lineconnectivity.h
              carconnectivity.h
              carcost.h
              cargrid.h
              gridcost3d.h
              carcost.h
              grid2d.h
              grid2dconnectivity.h
              gridtravcost.h
              gridsearch.h
              gridsearch3d.h
              spline.h
              DESTINATION include/dsl)              
