set (SOURCES
  main.cpp
  Parameters.cpp
  Glauber.cpp
  Util.cpp
  )

if (test)
    add_executable (unittest_Util.e Util.cpp)
    install(TARGETS unittest_Util.e DESTINATION ${CMAKE_HOME_DIRECTORY})
else (test)
    add_library(3dMCGlb_lib SHARED ${SOURCES})
    add_executable(3dMCGlb.e ${SOURCES})
    target_link_libraries (3dMCGlb_lib ${LIBS})
    target_link_libraries (3dMCGlb.e ${LIBS})
    install(TARGETS 3dMCGlb.e DESTINATION ${CMAKE_HOME_DIRECTORY})
endif (test)
