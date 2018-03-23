set (SOURCES
  main.cpp
  Parameters.cpp
  ParametersMap.cpp
  Glauber.cpp
  Util.cpp
  )

if (test)
    add_executable (unittest_Util.e Util.cpp UnitTest_Util.cc)
    install(TARGETS unittest_Util.e DESTINATION ${CMAKE_HOME_DIRECTORY})
    add_executable (unittest_ParametersMap.e ParametersMap.cpp Util.cpp UnitTest_ParametersMap.cc)
    install(TARGETS unittest_ParametersMap.e DESTINATION ${CMAKE_HOME_DIRECTORY})
    add_executable (unittest_Parameters.e Parameters.cpp ParametersMap.cpp Util.cpp UnitTest_Parameters.cc)
    install(TARGETS unittest_Parameters.e DESTINATION ${CMAKE_HOME_DIRECTORY})
else (test)
    add_library(3dMCGlb_lib SHARED ${SOURCES})
    add_executable(3dMCGlb.e ${SOURCES})
    target_link_libraries (3dMCGlb_lib ${LIBS})
    target_link_libraries (3dMCGlb.e ${LIBS})
    install(TARGETS 3dMCGlb.e DESTINATION ${CMAKE_HOME_DIRECTORY})
endif (test)
