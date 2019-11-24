# Find Doxygen package.
find_package(Doxygen)

if(DOXYGEN_FOUND)
   doxygen_add_docs(doc
      ${CMAKE_CURRENT_SOURCE_DIR}/README.md ${CMAKE_SOURCE_DIR}/include/tu ${CMAKE_SOURCE_DIR}/src/tu
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating doxygen documentation in ${CMAKE_CURRENT_BINARY_DIR}.")
else()
   message(WARNING "Could not find Doxygen. Generating a documentation is not supported.")
endif()
