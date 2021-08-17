# Find Doxygen package.
find_package(Doxygen)

if(DOXYGEN_FOUND)
  file(COPY
    ${CMAKE_SOURCE_DIR}/doc/doxygen-awesome-css
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
  set(DOXYGEN_EXTRA_PACKAGES "amsmath" "amssymb" "mathtools" "amsfont")
  set(DOXYGEN_EXTRACT_ALL "yes")
  set(DOXYGEN_EXTRACT_STATIC "yes")
  set(DOXYGEN_USE_MATHJAX "yes")
  set(DOXYGEN_INCLUDE_GRAPH "no")
  set(DOXYGEN_INCLUDED_BY_GRAPH "no")
  set(DOXYGEN_HTML_EXTRA_STYLESHEET "doxygen-awesome-css/doxygen-awesome.css")
  doxygen_add_docs(doc
    ${CMAKE_CURRENT_SOURCE_DIR}/index.md
    ${CMAKE_CURRENT_SOURCE_DIR}/file-formats.md
    ${CMAKE_CURRENT_SOURCE_DIR}/camion.md
    ${CMAKE_CURRENT_SOURCE_DIR}/tu.md
    ${CMAKE_CURRENT_SOURCE_DIR}/network.md
    ${CMAKE_CURRENT_SOURCE_DIR}/regular.md
    ${CMAKE_CURRENT_SOURCE_DIR}/graphic.md
    ${CMAKE_CURRENT_SOURCE_DIR}/series-parallel.md
    ${CMAKE_CURRENT_SOURCE_DIR}/balanced.md
    ${CMAKE_CURRENT_SOURCE_DIR}/perfect.md
    ${CMAKE_CURRENT_SOURCE_DIR}/consecutive-ones.md
    ${CMAKE_CURRENT_SOURCE_DIR}/conversions.md
    ${CMAKE_SOURCE_DIR}/include/cmr
    ${CMAKE_SOURCE_DIR}/src/cmr
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating doxygen documentation in ${CMAKE_CURRENT_BINARY_DIR}.")
else()
  message(WARNING "Could not find Doxygen. Generating a documentation is not supported.")
endif()
