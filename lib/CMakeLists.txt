add_library(dsl SHARED fibheap.c 
                       vertex.cc
                       edge.cc
                       graph.cc
                       cost.cc
                       search.cc
                       gridcost.cc
                       gridcost3d.cc
                       gridcost3dbezier.cc
                       gridcost3dspline.cc
                  #     gridcost3dprm.cc
                       gridsearch.cc
                       gridsearch3d.cc
                       gridsearch3dvel.cc
                       gridsearch3dvelb.cc
                       gridsearch3dvels.cc)


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
              gridcost.h
              gridcost3d.h
              gridcost3dbezier.h
              gridcost3dspline.h
           #   gridcost3dprm.h
              gridsearch.h
              gridsearch3d.h
              gridsearch3dvel.h
              gridsearch3dvelb.h
              gridsearch3dvels.h
              spline.h
              DESTINATION include/dsl)              
