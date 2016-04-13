#-----------------------------------------------------------------------
# gEMfitter
gem_find_source(gEMfitter_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
add_executable (gEMfitter ${gEMfitter_SOURCES})

target_link_libraries(gEMfitter ${EXTRA_LIBS})

if(NOT GEM_LINK_STATIC)
    set_target_properties(gEMfitter PROPERTIES INSTALL_RPATH "${GEMEXE_INSTALL_RPATH}"
                                               INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

#-----------------------------------------------------------------------
# Installation
install(TARGETS gEMfitter
        RUNTIME DESTINATION bin        COMPONENT Runtime
        LIBRARY DESTINATION libraries  COMPONENT Runtime
        ARCHIVE DESTINATION libraries  COMPONENT Development)

#-----------------------------------------------------------------------
# Testing
set(TestName          gEMfitter)
add_test(${TestName}  ./gEMfitter)
set_tests_properties(${TestName}
                     PROPERTIES    PASS_REGULAR_EXPRESSION  "PARSE ERROR:")
