#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.13)
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
#include("${ROOT_DIR}/modules/RootNewMacros.cmake")  # ROOT5
include("${ROOT_DIR}/ROOTUseFile.cmake")   # ROOT6

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})
include(${ROOT_USE_FILE})

#----------------------------------------------------------------------------
# Find yaml-cpp (required package)
#
#set(YAML_CPP_INCLUDE_DIR /opt/homebrew/opt/yaml-cpp/include)
#set(YAML_CPP_LIBRARY /opt/homebrew/opt/yaml-cpp/lib/libyaml-cpp.dylib)
find_package(yaml-cpp REQUIRED)

#----------------------------------------------------------------------------
# General compilation and linkage options
#
option(ENABLE_SANITIZER "compile and link with sanitizing flags" OFF)
option(ENABLE_GPROF "compile and link with gprof flags" OFF)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE DEBUG)
endif()
if(CMAKE_BUILD_TYPE STREQUAL DEBUG)
    add_compile_options(-O0 -fno-omit-frame-pointer)
    if(ENABLE_SANITIZER)
        set(SANITIZER_OPTIONS -fno-sanitize-recover=all)
        list(APPEND SANITIZER_OPTIONS -fsanitize=address)
        list(APPEND SANITIZER_OPTIONS -fsanitize=leak)
        list(APPEND SANITIZER_OPTIONS -fsanitize=undefined)
        list(APPEND SANITIZER_OPTIONS -fsanitize=float-divide-by-zero)
        add_compile_options(${SANITIZER_OPTIONS})
        add_link_options(${SANITIZER_OPTIONS})
    endif()
else()
    add_compile_options(-O3)
endif()
if(ENABLE_GPROF)
    set(GPROF_OPTIONS -pg)
    add_compile_options(${GPROF_OPTIONS})
    add_link_options(${GPROF_OPTIONS})
endif()
add_compile_options(-g -Wall -Wshadow -Wextra -Werror)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(SYSTEM ${Geant4_INCLUDE_DIR} ${ROOT_INCLUDE_DIRS} ${YAML_CPP_INCLUDE_DIR})
file(GLOB sources CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(muPos muPos.cc ${sources} ${headers})
target_link_libraries(muPos ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} ${YAML_CPP_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build muPos. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
file(GLOB scripts CONFIGURE_DEPENDS *.sh *.mac)
foreach(script ${scripts})
    get_filename_component(script ${script} NAME)
    configure_file(${PROJECT_SOURCE_DIR}/${script} ${PROJECT_BINARY_DIR}/${script} COPYONLY)
endforeach()
