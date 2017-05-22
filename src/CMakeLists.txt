# The extension is already found. Any number of sources could be listed here.



add_library (SBNfitlib SBNconfig.cxx SBNcovar.cxx SBNspec.cxx SBNfit.cxx SBNchi.cxx prob.cxx SBNfit3pN.cxx SBNosc.cxx mdict_h.so )
# Make sure the compiler can find include files for our Hello library
# when other libraries or executables link to Hello


include_directories ( ${CMAKE_SOURCE_DIR}/inc/)
include_directories ( ${CMAKE_SOURCE_DIR}/src/)
include_directories ( ${CMAKE_SOURCE_DIR}/tinyxml/)

#ROOT_GENERATE_DICTIONARY(LINKDEF mdict.h )
#target_link_libraries(SBNfitlib mdict)
target_link_libraries(SBNfitlib ${ROOT_LIBRARIES})
target_link_libraries(SBNfitlib tinyxmllib)
