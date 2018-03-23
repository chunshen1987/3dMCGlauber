set (SOURCES
  main.cpp
  Parameters.cpp
  Glauber.cpp
  )

add_library(3dMCGlb_lib SHARED ${SOURCES})
add_executable(3dMCGlb.e ${SOURCES})
target_link_libraries (3dMCGlb_lib ${LIBS})
target_link_libraries (3dMCGlb.e ${LIBS})

install(TARGETS 3dMCGlb.e DESTINATION ${CMAKE_HOME_DIRECTORY})
