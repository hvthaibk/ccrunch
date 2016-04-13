gem_find_source(GemData_SOURCES     "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
gem_find_source(GemDataCuda_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/CUDA/*.cpp")

if(NOT CUDA_FOUND)
    foreach(file_path ${GemDataCuda_SOURCES})
        list(REMOVE_ITEM GemData_SOURCES ${file_path})
    endforeach()
endif()

if(GEM_LINK_STATIC)
    add_library(GemData STATIC ${GemData_SOURCES})
else()
    add_library(GemData SHARED ${GemData_SOURCES})
endif()

if(CUDA_FOUND)
    cuda_add_cublas_to_target(GemData)
endif()

target_link_libraries(GemData ${EXTRA_LIBS})

if(NOT GEM_LINK_STATIC)
    set_target_properties(GemData PROPERTIES INSTALL_RPATH "${GEMLIB_INSTALL_RPATH}"
                                                   INSTALL_RPATH_USE_LINK_PATH TRUE)
    install(TARGETS       GemData DESTINATION ${PROJECT_LIB_DIR})
endif()
