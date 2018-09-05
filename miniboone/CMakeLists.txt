# "demo.cxx" and "demo_b.cxx". The extensions are automatically found.
add_executable (minirepo miniboone.cxx)
add_executable (minirepoL minibooneL.cxx)

# Link the executable to the Hello library. Since the Hello library has
# public include directories we will use those link directories when building
# helloDemo
include_directories ( ${CMAKE_SOURCE_DIR}/inc)
target_link_libraries (minirepo SBNfitlib)
target_link_libraries (minirepo tinyxmllib)
target_link_libraries (minirepoL SBNfitlib)
target_link_libraries (minirepoL tinyxmllib)

