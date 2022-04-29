if (NOT RStarTree_FOUND)
  find_path(RStarTree_INCLUDE_DIR NAMES "RStarTree/RStarTree.h" HINTS ENV RSTARTREE_INC PATH_SUFFIXES include)
  mark_as_advanced(RStarTree_INCLUDE_DIR)
  if (RStarTree_INCLUDE_DIR)
    string(REGEX REPLACE "^v?([0-9._-]+).*$" "\\1" RStarTree_VERSION "$ENV{RSTARTREE_VERSION}")
    string(REGEX REPLACE "[_-]" "." RStarTree_VERSION "${RStarTree_VERSION}")
    if (CETMODULES_CURRENT_PROJECT_NAME AND
        ${CETMODULES_CURRENT_PROJECT_NAME}_OLD_STYLE_CONFIG_VARS)
      include_directories("${RStarTree_INCLUDE_DIR}")
    endif()
  endif()
endif()

include(FindPackageHandleStandardArgs)

if (NOT "${RStarTree_VERSION}" STREQUAL "")
  set(_cet_frst_version_args VERSION_VAR RStarTree_VERSION)
else()
  unset(_cet_frst_version_args)
endif()

find_package_handle_standard_args(RStarTree
  REQUIRED_VARS RStarTree_INCLUDE_DIR
  ${cet_frst_version_args}
)

unset(_cet_frst_version_args)

if (RStarTree_FOUND AND NOT TARGET RStarTree::RStarTree)
  add_library(RStarTree::RStarTree INTERFACE IMPORTED)
  set_target_properties(RStarTree::RStarTree PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${RStarTree_INCLUDE_DIR}"
  )
endif()
