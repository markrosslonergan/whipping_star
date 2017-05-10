# "demo.cxx" and "demo_b.cxx". The extensions are automatically found.
add_executable (exampledemo examples.cxx)

# Link the executable to the Hello library. Since the Hello library has
# public include directories we will use those link directories when building
# helloDemo
target_include_directories (exampledemo PUBLIC ${PROJECT_BINARY_DIR}/inc)
target_link_libraries (exampledemo LINK_PUBLIC SBNfitlib)
