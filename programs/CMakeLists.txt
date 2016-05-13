gem_find_dir(DirPrograms "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

list(REMOVE_ITEM DirPrograms "${CMAKE_CURRENT_SOURCE_DIR}")
if(GEM_LINK_STATIC)
    list(REMOVE_ITEM DirPrograms "${CMAKE_CURRENT_SOURCE_DIR}/gEMaligner")
    list(REMOVE_ITEM DirPrograms "${CMAKE_CURRENT_SOURCE_DIR}/gEMpython")
endif()
if(GEM_LINK_STATIC OR ((NOT GEM_LINK_STATIC) AND (NOT __GEM_USE_QT4__)))
    list(REMOVE_ITEM DirPrograms "${CMAKE_CURRENT_SOURCE_DIR}/gEMcrystal")
endif()

foreach(dir ${DirPrograms})
    add_subdirectory(${dir})
endforeach()
