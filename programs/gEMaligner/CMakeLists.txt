#-----------------------------------------------------------------------
# gEMaligner
gem_find_source(gEMaligner_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
add_executable (gEMaligner ${gEMaligner_SOURCES})

target_link_libraries(gEMaligner ${EXTRA_LIBS})

if(NOT GEM_LINK_STATIC)
    set_target_properties(gEMaligner PROPERTIES INSTALL_RPATH "${GEMEXE_INSTALL_RPATH}"
                                                INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

#-----------------------------------------------------------------------
# Installation
install(TARGETS gEMaligner
        RUNTIME DESTINATION bin        COMPONENT Runtime
        LIBRARY DESTINATION libraries  COMPONENT Runtime
        ARCHIVE DESTINATION libraries  COMPONENT Development)

#-----------------------------------------------------------------------
# Testing
set(TestName          gEMaligner)
add_test(${TestName}  ./gEMaligner)
set_tests_properties(${TestName}
                     PROPERTIES    PASS_REGULAR_EXPRESSION  "PARSE ERROR:")
