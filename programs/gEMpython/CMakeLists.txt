gem_find_source(GemPython_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

if(NOT GEM_LINK_STATIC)
    add_library(GemPython MODULE ${GemPython_SOURCES})

    target_link_libraries(GemPython ${BOOST_LIBRARY} ${PYTHON_LIBRARIES})

    set_target_properties(GemPython PROPERTIES INSTALL_RPATH "${GEMLIB_INSTALL_RPATH}"
                                               INSTALL_RPATH_USE_LINK_PATH TRUE)
    install(TARGETS       GemPython DESTINATION ${PROJECT_LIB_DIR})
endif()
