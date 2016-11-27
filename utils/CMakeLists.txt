################################################################################
#  EXAMPLES 
################################################################################

# We add the compilation options for the library
add_definitions (${gatb-core-flags} " -Wno-ambiguous-member-template")

include_directories (${gatb-core-includes})

#list (APPEND PROGRAMS dsk2ascii filter_by_previous_contigs)
list (APPEND PROGRAMS dsk2ascii)

FOREACH (program ${PROGRAMS})
  add_executable(${program} ${program}.cpp)
  target_link_libraries(${program} ${gatb-core-libraries})
ENDFOREACH (program)

################################################################################
#  INSTALLATION 
################################################################################
#INSTALL (TARGETS dsk2ascii filter_by_previous_contigs DESTINATION bin)
INSTALL (TARGETS dsk2ascii DESTINATION bin)
