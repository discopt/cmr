# Find Doxygen package.
find_package(Doxygen)

if(DOXYGEN_FOUND)
  set(DOXYGEN_EXTRA_PACKAGES "amsmath" "amssymb")
  set(DOXYGEN_EXTRACT_ALL "yes")
  set(DOXYGEN_EXTRACT_STATIC "yes")
  set(DOXYGEN_USE_MATHJAX "yes")
  doxygen_add_docs(doc
    ${CMAKE_CURRENT_SOURCE_DIR}/README.md ${CMAKE_SOURCE_DIR}/include/cmr ${CMAKE_SOURCE_DIR}/src/cmr
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating doxygen documentation in ${CMAKE_CURRENT_BINARY_DIR}.")
else()
  message(WARNING "Could not find Doxygen. Generating a documentation is not supported.")
endif()
