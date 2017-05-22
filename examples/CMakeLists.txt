# "demo.cxx" and "demo_b.cxx". The extensions are automatically found.
add_executable (exampledemo examples.cxx)
add_executable (mcdemo testmc.cxx)

# Link the executable to the Hello library. Since the Hello library has
# public include directories we will use those link directories when building
# helloDemo
target_include_directories (exampledemo PUBLIC ${PROJECT_BINARY_DIR}/inc)
target_link_libraries (exampledemo LINK_PUBLIC SBNfitlib)

#ROOT_GENERATE_DICTIONARY(H__vec mymap.h LINKDEF mymapLinkDef.h)

target_include_directories (mcdemo PUBLIC ${PROJECT_BINARY_DIR}/inc)
target_link_libraries (mcdemo LINK_PUBLIC SBNfitlib)
