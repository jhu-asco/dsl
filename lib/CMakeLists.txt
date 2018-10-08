add_library(dsl SHARED fibheap.c
                       search.cc
                       utils.cc
                       utilsimg.cc
                       grid2d.cc
                       grid3d.cc
                       cargrid.cc
                       grid2dconnectivity.cc
                       grid3dconnectivity.cc
                       carconnectivity.cc
                       carconnectivity2.cc
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
              lpastar.h
              dstar.h
              utils.h
              utilsimg.h
              cell.h
              grid.h
              gridcost.h
              gridconnectivity.h
              gridpath.h
              lineconnectivity.h
              carconnectivity.h
              carconnectivity2.h
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
