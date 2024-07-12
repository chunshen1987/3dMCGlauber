# Install script for directory: /home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/3dMCGlb.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/3dMCGlb.e")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/3dMCGlb.e"
         RPATH "/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/lib:/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/LHAPDF_Lib/lib/")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/3dMCGlb.e")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber" TYPE EXECUTABLE FILES "/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/build/src/3dMCGlb.e")
  if(EXISTS "$ENV{DESTDIR}/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/3dMCGlb.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/3dMCGlb.e")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/3dMCGlb.e"
         OLD_RPATH "/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/LHAPDF_Lib/lib:/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/build/src:"
         NEW_RPATH "/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/lib:/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/LHAPDF_Lib/lib/")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/3dMCGlb.e")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/Metropolis.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/Metropolis.e")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/Metropolis.e"
         RPATH "/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/lib:/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/LHAPDF_Lib/lib/")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/Metropolis.e")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber" TYPE EXECUTABLE FILES "/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/build/src/Metropolis.e")
  if(EXISTS "$ENV{DESTDIR}/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/Metropolis.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/Metropolis.e")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/Metropolis.e"
         OLD_RPATH "/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/LHAPDF_Lib/lib:/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/build/src:"
         NEW_RPATH "/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/lib:/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/LHAPDF_Lib/lib/")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/Metropolis.e")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/lib3dMCGlb.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/lib3dMCGlb.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/lib3dMCGlb.so"
         RPATH "/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/lib:/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/LHAPDF_Lib/lib/")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/build/src/lib3dMCGlb.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/lib3dMCGlb.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/lib3dMCGlb.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/lib3dMCGlb.so"
         OLD_RPATH ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
         NEW_RPATH "/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/lib:/home/wenbin/Downloads/Wenbin_working/Work/WSU_BNL_work/UPC/code/my_branch/iEBE-MUSIC/codes/wenbin_main/3dMCGlauber/LHAPDF_Lib/lib/")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/lib3dMCGlb.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

