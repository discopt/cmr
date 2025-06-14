# Find Doxygen package.
find_package(Doxygen)

if(DOXYGEN_FOUND)
  file(COPY
    ${CMAKE_SOURCE_DIR}/doc/doxygen-awesome-css
    ${CMAKE_SOURCE_DIR}/doc/header.html
    ${CMAKE_SOURCE_DIR}/doc/logo.png
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
  file(COPY
    ${CMAKE_SOURCE_DIR}/doc/header.html
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
  set(DOXYGEN_EXTRA_PACKAGES "amsmath" "amssymb" "mathtools" "amsfont")
  set(DOXYGEN_EXTRACT_ALL "yes")
  set(DOXYGEN_EXTRACT_STATIC "yes")
  set(DOXYGEN_USE_MATHJAX "yes")
  set(DOXYGEN_MATHJAX_RELPATH "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/")
  set(DOXYGEN_INCLUDE_GRAPH "no")
  set(DOXYGEN_INCLUDED_BY_GRAPH "no")
  set(DOXYGEN_HTML_HEADER "header.html")
  set(DOXYGEN_PROJECT_LOGO "logo.png")
  set(DOXYGEN_HTML_EXTRA_STYLESHEET "doxygen-awesome-css/doxygen-awesome.css")
  doxygen_add_docs(doc
    ${CMAKE_CURRENT_SOURCE_DIR}/index.md
    ${CMAKE_CURRENT_SOURCE_DIR}/balanced.md
    ${CMAKE_CURRENT_SOURCE_DIR}/build.md
    ${CMAKE_CURRENT_SOURCE_DIR}/camion.md
    ${CMAKE_CURRENT_SOURCE_DIR}/changes.md
    ${CMAKE_CURRENT_SOURCE_DIR}/consecutive-ones.md
    ${CMAKE_CURRENT_SOURCE_DIR}/ctu.md
    ${CMAKE_CURRENT_SOURCE_DIR}/k_ary.md
    ${CMAKE_CURRENT_SOURCE_DIR}/file-formats.md
    ${CMAKE_CURRENT_SOURCE_DIR}/graphic.md
    ${CMAKE_CURRENT_SOURCE_DIR}/generators.md
    ${CMAKE_CURRENT_SOURCE_DIR}/equimodular.md
    ${CMAKE_CURRENT_SOURCE_DIR}/matroids.md
    ${CMAKE_CURRENT_SOURCE_DIR}/max-flow-min-cut.md
    ${CMAKE_CURRENT_SOURCE_DIR}/named.md
    ${CMAKE_CURRENT_SOURCE_DIR}/network.md
    ${CMAKE_CURRENT_SOURCE_DIR}/perfect.md
    ${CMAKE_CURRENT_SOURCE_DIR}/regular.md
    ${CMAKE_CURRENT_SOURCE_DIR}/series-parallel.md
    ${CMAKE_CURRENT_SOURCE_DIR}/seymour.md
    ${CMAKE_CURRENT_SOURCE_DIR}/totally-balanced.md
    ${CMAKE_CURRENT_SOURCE_DIR}/tu.md
    ${CMAKE_CURRENT_SOURCE_DIR}/utilities.md
    ${CMAKE_SOURCE_DIR}/include/cmr
    ${CMAKE_SOURCE_DIR}/src/cmr
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating doxygen documentation in ${CMAKE_CURRENT_BINARY_DIR}.")
else()
  message(WARNING "Could not find Doxygen. Generating a documentation is not supported.")
endif()
