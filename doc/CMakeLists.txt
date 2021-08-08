# Find Doxygen package.
find_package(Doxygen)

if(DOXYGEN_FOUND)
  file(COPY
    ${CMAKE_SOURCE_DIR}/doc/doxygen-awesome-css
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
  set(DOXYGEN_EXTRA_PACKAGES "amsmath" "amssymb" "mathtools" "amsfont")
#  set(DOXYGEN_EXTRACT_ALL "yes")
  set(DOXYGEN_EXTRACT_STATIC "yes")
  set(DOXYGEN_USE_MATHJAX "yes")
  set(DOXYGEN_HTML_EXTRA_STYLESHEET "doxygen-awesome-css/doxygen-awesome.css")
  doxygen_add_docs(doc
    ${CMAKE_CURRENT_SOURCE_DIR}/INDEX.md
    ${CMAKE_CURRENT_SOURCE_DIR}/FILEFORMATS.md
    ${CMAKE_CURRENT_SOURCE_DIR}/TU.md
    ${CMAKE_CURRENT_SOURCE_DIR}/NETWORK.md
    ${CMAKE_CURRENT_SOURCE_DIR}/REGULAR.md
    ${CMAKE_CURRENT_SOURCE_DIR}/GRAPHIC.md
    ${CMAKE_CURRENT_SOURCE_DIR}/SP.md
    ${CMAKE_CURRENT_SOURCE_DIR}/BALANCED.md
    ${CMAKE_CURRENT_SOURCE_DIR}/PERFECT.md
    ${CMAKE_SOURCE_DIR}/include/cmr
    ${CMAKE_SOURCE_DIR}/src/cmr
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating doxygen documentation in ${CMAKE_CURRENT_BINARY_DIR}.")
else()
  message(WARNING "Could not find Doxygen. Generating a documentation is not supported.")
endif()
