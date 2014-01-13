project(dsk)

cmake_minimum_required (VERSION 2.8.11.2)

################################################################################
# GATB-CORE
################################################################################

# we depend on gatb-core; here, we define where to find all the required material
add_subdirectory (../../thirdparty/gatb-core/gatb-core  "${CMAKE_CURRENT_BINARY_DIR}/ext/gatb-core")


################################################################################
# TOOL
################################################################################

# we get compilation definitions from the gatb-core part
add_definitions (${gatb-core-flags})

# we give the source directories: from project source and from GATB-CORE source
include_directories (src  ${gatb-core-includes} )

# we define the files to be compiled
file (GLOB_RECURSE  ProjectFiles  src/*)

# we define the artifact to be built: the project binary
add_executable        (${CMAKE_PROJECT_NAME}  ${ProjectFiles})

# we define which libraries to be linked with project binary
target_link_libraries (${CMAKE_PROJECT_NAME}  ${gatb-core-libraries})

