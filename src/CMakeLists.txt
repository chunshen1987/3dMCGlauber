set (SOURCES
    Parameters.cpp
    ParametersMap.cpp
    EventGenerator.cpp
    Glauber.cpp
    MakeDensity.cpp
    Random.cpp
    Nucleus.cpp
    Nucleon.cpp
    CollisionEvent.cpp
    Util.cpp
    QCDString.cpp
    pretty_ostream.cpp
    ../eps09/eps09.cpp
    )

add_library(${libname} SHARED ${SOURCES})
target_link_libraries (${libname})
if (link_with_lib)
    add_executable(${exename} main.cpp)
    target_link_libraries (${exename} ${libname})
else (link_with_lib)
    add_executable(${exename} main.cpp ${SOURCES})
    target_link_libraries (${exename})
endif (link_with_lib)
add_executable(Metropolis.e SampleValenceQuarkMetropolis.cpp)
target_link_libraries (Metropolis.e ${libname} LHAPDF)
if (APPLE)
    set_target_properties(${libname} PROPERTIES COMPILE_FLAGS "-DAPPLE")
    set_target_properties (${exename} PROPERTIES COMPILE_FLAGS "-DAPPLE")
endif (APPLE)
set_target_properties(${libname}
    PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${CMAKE_HOME_DIRECTORY}/LHAPDF_Lib/lib/")
set_target_properties(${exename}
    PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${CMAKE_HOME_DIRECTORY}/LHAPDF_Lib/lib/")
set_target_properties(Metropolis.e
    PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${CMAKE_HOME_DIRECTORY}/LHAPDF_Lib/lib/")
install(TARGETS ${exename} Metropolis.e ${libname}
        RUNTIME DESTINATION ${CMAKE_HOME_DIRECTORY}
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static)

if (unittest)
    add_executable (unittest_CollisionEvent.e CollisionEvent_UnitTest.cpp)
    target_link_libraries (unittest_CollisionEvent.e ${libname} ${LIBS})
    add_executable (unittest_Quark.e Quark_UnitTest.cpp)
    add_executable (unittest_Nucleon.e Nucleon_UnitTest.cpp)
    target_link_libraries (unittest_Nucleon.e ${libname} ${LIBS})
    add_executable (unittest_Nucleus.e Nucleus_UnitTest.cpp)
    target_link_libraries (unittest_Nucleus.e ${libname} ${LIBS})
    add_executable (unittest_Particle.e Particle_UnitTest.cpp)
    add_executable (unittest_Util.e Util_UnitTest.cpp)
    target_link_libraries (unittest_Util.e ${libname} ${LIBS})
    add_executable (unittest_Random.e Random_UnitTest.cpp)
    target_link_libraries (unittest_Random.e ${libname} ${LIBS})
    add_executable (unittest_ParametersMap.e ParametersMap_UnitTest.cpp)
    target_link_libraries (unittest_ParametersMap.e ${libname} ${LIBS})
    add_executable (unittest_Parameters.e Parameters_UnitTest.cpp)
    target_link_libraries (unittest_Parameters.e ${libname} ${LIBS})
    install(TARGETS unittest_CollisionEvent.e unittest_Quark.e
                    unittest_Nucleon.e unittest_Nucleus.e unittest_Particle.e
                    unittest_Util.e unittest_Random.e unittest_ParametersMap.e
                    unittest_Parameters.e
                    ${libname}
            RUNTIME DESTINATION ${CMAKE_HOME_DIRECTORY}/test
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib/static)
endif (unittest)
