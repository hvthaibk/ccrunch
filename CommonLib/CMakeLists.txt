gem_find_source(GemCommon_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

if(GEM_LINK_STATIC)
    add_library(GemCommon STATIC ${GemCommon_SOURCES})
else()
    add_library(GemCommon SHARED ${GemCommon_SOURCES})
endif()

target_link_libraries(GemCommon ${BOOST_LIBRARY} ${FFTW_LIBRARY})

if(NOT GEM_LINK_STATIC)
    set_target_properties(GemCommon PROPERTIES INSTALL_RPATH "${GEMLIB_INSTALL_RPATH}"
                                               INSTALL_RPATH_USE_LINK_PATH TRUE)
    install(TARGETS       GemCommon DESTINATION ${PROJECT_LIB_DIR})
endif()
