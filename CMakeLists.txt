cmake_minimum_required (VERSION 2.6)
project (SBNfit)
# The version number.
set (SBNfit_VERSION_MAJOR 2)
set (SBNfit_VERSION_MINOR 0)

if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
          set (CMAKE_CXX_FLAGS "--std=gnu++11 ${CMAKE_CXX_FLAGS}")
    endif ()
    else ()
	    set (CMAKE_CXX_STANDARD 11)
endif ()


list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package(ROOT REQUIRED COMPONENTS Core Cint RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread)
include(${ROOT_USE_FILE})


# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
#ROOT_GENERATE_DICTIONARY(G__Event MCEventWeight.h LINKDEF MCEventWeightLinkDef.h)

#add_subdirectory(inc)
add_subdirectory(tinyxml)
add_subdirectory(src)
add_subdirectory(examples)

