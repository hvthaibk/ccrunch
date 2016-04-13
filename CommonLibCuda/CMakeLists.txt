gem_find_source(GemCuda_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cu")

if(GEM_LINK_STATIC)
    cuda_add_library(GemCuda STATIC ${GemCuda_SOURCES})
else()
    cuda_add_library(GemCuda SHARED ${GemCuda_SOURCES})
endif()

target_link_libraries(GemCuda GemCommon)

cuda_add_cufft_to_target (GemCuda)

if(NOT GEM_LINK_STATIC)
    set_target_properties(GemCuda PROPERTIES INSTALL_RPATH "${GEMLIB_INSTALL_RPATH}"
                                                   INSTALL_RPATH_USE_LINK_PATH TRUE)
    install(TARGETS       GemCuda DESTINATION ${PROJECT_LIB_DIR})
endif()
