#-----------------------------------------------------------------------
# gEMcrystal
project(gEMcrystal)
find_package(Qt4)
set(main_SOURCES cData3Qt.cpp functions.cpp loader.cpp main.cpp main_view.cpp main_window.cpp mask.cpp micrograph.cpp output_text_edit.cpp particle_label.cpp picker.cpp point.cpp result.cpp result_extra.cpp scene.cpp)
set(main_HEADERS cData3Qt.h loader.h main_window.h main_view.h micrograph.h output_text_edit.h particle_label.h picker.h point.h result.h scene.h)
QT4_WRAP_CPP(main_HEADERS_MOC ${main_HEADERS})
include(${QT_USE_FILE})
add_definitions(${QT_DEFINITIONS})
add_executable(gEMcrystal  ${main_SOURCES}
                    ${main_HEADERS_MOC})
target_link_libraries(gEMcrystal ${EXTRA_LIBS})
set_target_properties(gEMcrystal PROPERTIES INSTALL_RPATH "${GEMEXE_INSTALL_RPATH}"
                                            INSTALL_RPATH_USE_LINK_PATH TRUE)
#-----------------------------------------------------------------------
# Installation
install(TARGETS gEMcrystal
        RUNTIME DESTINATION bin        COMPONENT Runtime
        LIBRARY DESTINATION libraries  COMPONENT Runtime
        ARCHIVE DESTINATION libraries  COMPONENT Development)

#-----------------------------------------------------------------------
# Testing
#set(TestName          gEMcrystal)
#add_test(${TestName}  ./gEMcrystal)
#set_tests_properties(${TestName}
#                     PROPERTIES    PASS_REGULAR_EXPRESSION  "PARSE ERROR:")

