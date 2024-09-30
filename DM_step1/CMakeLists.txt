#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.9 FATAL_ERROR)
project(muPos)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()
#----------------------------------------------------------------------------
# Find ROOT (required package)
#
find_package(ROOT REQUIRED)
# Include ROOT's CMake functions for dictionary generation
#include("${ROOT_DIR}/modules/RootNewMacros.cmake")  #ROOT5
include("${ROOT_DIR}/ROOTUseFile.cmake")   #ROOT6

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIRS}
                    )
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(muPos muPos.cc ${sources} ${headers})
target_link_libraries(muPos ${Geant4_LIBRARIES}  ${ROOT_LIBRARIES} )

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build muPos. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(muPos_SCRIPTS
   vis.mac  CryMu.mac run.sh run_cry.sh
  )

foreach(_script ${muPos_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()




