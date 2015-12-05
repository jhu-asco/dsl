add_library(dsl SHARED fibheap.c 
#                       vertex.cc
#                       edge.cc
#                       graph.cc
#                       cost.cc
                       search.cc
                       gridcost.cc
                       gridcost3d.cc
                       gridtravcost.cc
                       gridsearch.cc
                       gridsearch3d.cc
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
              gridcost.h
              cell2d.h
              gridcost3d.h
              gridtravcost.h
              gridsearch.h
              gridsearch3d.h
              spline.h
              DESTINATION include/dsl)              
