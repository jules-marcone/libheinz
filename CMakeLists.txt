cmake_minimum_required(VERSION 3.6 FATAL_ERROR)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
include(PreventInSourceBuilds)

project(LibHeinz VERSION 0.1 LANGUAGES CXX)

## Options.

if(NOT DEFINED BUILD_SHARED_LIBS)
    option(BUILD_SHARED_LIBS "Build as shared library" ON)
endif()
option(PEDANTIC "Compile with pedantic warnings" ON)
option(WERROR "Treat warnings as errors" OFF)

## Compiler settings.

set(CMAKE_CXX_STANDARD 17)
add_compile_options(-O2)
if(WIN32)
    set(WINDOWS_EXPORT_ALL_SYMBOLS ON)
    message(STATUS "compiling for Windows")
else()
    add_compile_options(-Wall -Wno-unknown-pragmas)
    if(PEDANTIC)
        add_compile_options(-pedantic)
    endif()
endif()
if(WERROR)
    add_compile_options(-Werror)
endif()

## Install settings.

include(CTest)

add_subdirectory(inc)
add_subdirectory(test)

include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "${LibHeinz_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${LibHeinz_VERSION_MINOR}")
set(CPACK_SOURCE_GENERATOR "TGZ")
include(CPack)
include(CMakePackageConfigHelpers)
export(PACKAGE LibHeinz)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/LibHeinzConfig.cmake"
    INSTALL_DESTINATION "lib/cmake/example"
    NO_SET_AND_CHECK_MACRO
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
    )

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/LibHeinzConfigVersion.cmake"
    VERSION "${LibHeinz_VERSION_MAJOR}.${LibHeinz_VERSION_MINOR}"
    COMPATIBILITY AnyNewerVersion
    )

install(FILES
    ${PROJECT_BINARY_DIR}/LibHeinzConfig.cmake
    ${PROJECT_BINARY_DIR}/LibHeinzConfigVersion.cmake
    DESTINATION cmake COMPONENT dev)
