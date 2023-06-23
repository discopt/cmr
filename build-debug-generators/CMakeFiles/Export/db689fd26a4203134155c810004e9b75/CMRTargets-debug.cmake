#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "CMR::cmr" for configuration "Debug"
set_property(TARGET CMR::cmr APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(CMR::cmr PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libcmr.so"
  IMPORTED_SONAME_DEBUG "libcmr.so"
  )

list(APPEND _cmake_import_check_targets CMR::cmr )
list(APPEND _cmake_import_check_files_for_CMR::cmr "${_IMPORT_PREFIX}/lib/libcmr.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
