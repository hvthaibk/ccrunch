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

if(EXISTS "${PROJECT_SOURCE_DIR}/build64/gEMpicker")
    include(GetPrerequisites)
    get_prerequisites("${PROJECT_SOURCE_DIR}/build64/gEMpicker"
                      DEPENDENCIES 1 1 "" "${3RDPARTY_LIB_DIR}")
    foreach(DEPENDENCY ${DEPENDENCIES})
        get_filename_component(DEPENDENCY_NAME    "${DEPENDENCY}" NAME)
        get_filename_component(DEPENDENCY_ACTUAL  "${DEPENDENCY}" REALPATH)
        install(FILES       "${DEPENDENCY_ACTUAL}"
                DESTINATION "${CMAKE_INSTALL_PREFIX}/shared"  COMPONENT Runtime)
    endforeach()
endif(EXISTS "${PROJECT_SOURCE_DIR}/build64/gEMpicker")

#-----------------------------------------------------------------------
# Packaging
include(InstallRequiredSystemLibraries)
set(CPACK_GENERATOR "DEB;STGZ;ZIP")
set(CPACK_PACKAGE_VENDOR "INRIA")
set(CPACK_PACKAGE_NAME gEMpicker)
set(CPACK_PACKAGE_VERSION_MAJOR "${ccrunch_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${ccrunch_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${ccrunch_VERSION_PATCH}")
set(CPACK_PACKAGE_CONTACT "Thai V. Hoang <hvthaibk@gmail.com>")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "gEMpicker if a software developed at INRIA.")
#set(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6 (>= 2.11.0), libgcc1 (>= 1:4.4.3), libsdl1.2debian (>=1.2.14), libsdl-mixer1.2 (>=1.2.8), libsdl-net1.2 (>=1.2.7), bash")
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/License.txt")
set(CPACK_COMPONENTS_ALL_IN_ONE_PACKAGE 1)
#set(CPACK_SET_DESTDIR "ON")

set(BUILD_ARCH "${CMAKE_SYSTEM_PROCESSOR}")
if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
    set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
    set(CPACK_RPM_PACKAGE_ARCHITECTURE "x86_64" )
    set(BUILD_ARCH "${CPACK_DEBIAN_PACKAGE_ARCHITECTURE}")
endif(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
include(CPack)

#-----------------------------------------------------------------------
# Testing
set(TestName          gEMpicker)
add_test(${TestName}  ./gEMpicker)
set_tests_properties(${TestName}
                     PROPERTIES    PASS_REGULAR_EXPRESSION  "PARSE ERROR:")
