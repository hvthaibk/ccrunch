#-----------------------------------------------------------------------
# gEMpicker
gem_find_source(gEMpicker_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
add_executable (gEMpicker ${gEMpicker_SOURCES})

target_link_libraries(gEMpicker ${EXTRA_LIBS})

if(NOT GEM_LINK_STATIC)
    set_target_properties(gEMpicker PROPERTIES INSTALL_RPATH "${GEMEXE_INSTALL_RPATH}"
                                               INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

#-----------------------------------------------------------------------
# Installation
install(TARGETS gEMpicker
        RUNTIME DESTINATION bin        COMPONENT Runtime
        LIBRARY DESTINATION libraries  COMPONENT Runtime
        ARCHIVE DESTINATION libraries  COMPONENT Development)

#-----------------------------------------------------------------------
# Testing
set(TestName          gEMpicker)
add_test(${TestName}  ./gEMpicker)
set_tests_properties(${TestName}
                     PROPERTIES    PASS_REGULAR_EXPRESSION  "PARSE ERROR:")
