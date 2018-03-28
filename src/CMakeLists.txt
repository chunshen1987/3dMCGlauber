set (SOURCES
    main.cpp
    Parameters.cpp
    ParametersMap.cpp
    Glauber.cpp
    Random.cpp
    Nucleus.cpp
    Nucleon.cpp
    CollisionEvent.cpp
    Util.cpp
    )

if (test)
    add_executable (unittest_CollisionEvent.e CollisionEvent_UnitTest.cpp CollisionEvent.cpp Nucleon.cpp)
    install(TARGETS unittest_CollisionEvent.e DESTINATION ${CMAKE_HOME_DIRECTORY})
    add_executable (unittest_Quark.e Quark_UnitTest.cpp)
    install(TARGETS unittest_Quark.e DESTINATION ${CMAKE_HOME_DIRECTORY})
    add_executable (unittest_Nucleon.e Nucleon.cpp Nucleon_UnitTest.cpp)
    install(TARGETS unittest_Nucleon.e DESTINATION ${CMAKE_HOME_DIRECTORY})
    add_executable (unittest_Nucleus.e Nucleus.cpp Nucleon.cpp Random.cpp Nucleus_UnitTest.cpp)
    install(TARGETS unittest_Nucleus.e DESTINATION ${CMAKE_HOME_DIRECTORY})
    add_executable (unittest_Particle.e Particle_UnitTest.cpp)
    install(TARGETS unittest_Particle.e DESTINATION ${CMAKE_HOME_DIRECTORY})
    add_executable (unittest_Util.e Util.cpp Util_UnitTest.cpp)
    install(TARGETS unittest_Util.e DESTINATION ${CMAKE_HOME_DIRECTORY})
    add_executable (unittest_Random.e Random.cpp Random_UnitTest.cpp)
    install(TARGETS unittest_Random.e DESTINATION ${CMAKE_HOME_DIRECTORY})
    add_executable (unittest_ParametersMap.e ParametersMap.cpp Util.cpp ParametersMap_UnitTest.cpp)
    install(TARGETS unittest_ParametersMap.e DESTINATION ${CMAKE_HOME_DIRECTORY})
    add_executable (unittest_Parameters.e Parameters.cpp ParametersMap.cpp Util.cpp Parameters_UnitTest.cpp)
    install(TARGETS unittest_Parameters.e DESTINATION ${CMAKE_HOME_DIRECTORY})
else (test)
    add_library(3dMCGlb_lib SHARED ${SOURCES})
    target_link_libraries (3dMCGlb_lib ${LIBS})
    add_executable(3dMCGlb.e ${SOURCES})
    target_link_libraries (3dMCGlb.e ${LIBS})
    install(TARGETS 3dMCGlb.e DESTINATION ${CMAKE_HOME_DIRECTORY})
endif (test)
