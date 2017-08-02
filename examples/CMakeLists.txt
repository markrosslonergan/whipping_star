# "demo.cxx" and "demo_b.cxx". The extensions are automatically found.


add_executable (exampledemo examples.cxx)
add_executable (quickplot quickplot.cxx)
add_executable (mcdemo testmc.cxx)
add_executable (unit unit_test.cxx)
add_executable (gentest gentest.cxx)
add_executable (probtest probtest.cxx)

# Link the executable to the Hello library. Since the Hello library has
# public include directories we will use those link directories when building
# helloDemo
include_directories ( ${CMAKE_SOURCE_DIR}/inc)
target_link_libraries (exampledemo SBNfitlib)
target_link_libraries (exampledemo tinyxmllib)

#ROOT_GENERATE_DICTIONARY(H__vec mymap.h LINKDEF mymapLinkDef.h)
target_link_libraries (mcdemo SBNfitlib)
target_link_libraries (gentest SBNfitlib)
target_link_libraries (probtest SBNfitlib)

target_link_libraries (unit SBNfitlib)
target_link_libraries (unit tinyxmllib)

target_link_libraries (quickplot SBNfitlib)
target_link_libraries (quickplot tinyxmllib)

