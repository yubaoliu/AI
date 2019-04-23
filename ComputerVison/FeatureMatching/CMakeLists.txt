cmake_minimum_required(VERSION 3.0)
project( FeatureMatching)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")

file(GLOB_RECURSE srcs "src/*.cpp")
file(GLOB_RECURSE hdrs "src/*.h*")


find_package( OpenCV REQUIRED )
add_executable( FeatureMatching   ${srcs} ${hdrs}  )
target_link_libraries( FeatureMatching ${OpenCV_LIBS} )
