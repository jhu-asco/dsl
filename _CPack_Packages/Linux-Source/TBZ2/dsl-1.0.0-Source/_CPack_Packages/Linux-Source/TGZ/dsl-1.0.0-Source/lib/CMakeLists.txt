add_library(dsl SHARED fibheap.c 
                       vertex.cc
                       edge.cc
                       graph.cc
                       cost.cc
                       search.cc
                       gridcost.cc
                       gridsearch.cc)


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
              gridsearch.h
              DESTINATION include/dsl)              