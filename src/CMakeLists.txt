# The extension is already found. Any number of sources could be listed here.

add_library (SBNfitlib SBNconfig.cxx SBNcovar.cxx SBNspec.cxx SBNfit.cxx SBNchi.cxx prob.cxx SBNosc.cxx SBNfit3pN.cxx)
# Make sure the compiler can find include files for our Hello library
# when other libraries or executables link to Hello
target_include_directories (SBNfitlib PUBLIC ${PROJECT_BINARY_DIR}/inc)
target_link_libraries(SBNfitlib ${ROOT_LIBRARIES})
target_link_libraries(SBNfitlib tinyxmllib)
