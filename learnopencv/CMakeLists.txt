cmake_minimum_required(VERSION 3.0)
project(learnOpencv)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")


add_subdirectory(src)
