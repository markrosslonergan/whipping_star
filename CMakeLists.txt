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

#set (CMAKE_CXX_FLAGS "-L/grid/fermiapp/products/larsoft/root/v5_34_32/Linux64bit+2.6-2.12-e9-nu-prof/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -pthread -std=c++1y -Wno-deprecated-declarations -m64 -I/grid/fermiapp/products/larsoft/root/v5_34_32/Linux64bit+2.6-2.12-e9-nu-prof/include ")

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package(ROOT REQUIRED)# COMPONENTS Core Cint RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread)
include(${ROOT_USE_FILE})

message("Root lib: ${ROOT_LIBRARIES}")
#ls /home/mark/programs/root/binary/root/lib/*.so | awk '{printf "%s+",$0} END {print ""}'

## So these next two lines on duncan
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
#and this line for nevis
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lgsl -lgslcblas -lm -I/usr/include")

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -L/usr/local/lib -lglobes  -I/usr/local/include")


# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inc/)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/tinyxml/)
#ROOT_GENERATE_DICTIONARY(G__Event MCEventWeight.h LINKDEF MCEventWeightLinkDef.h)

#add_subdirectory(inc)
add_subdirectory(tinyxml)
add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(miniboone)

