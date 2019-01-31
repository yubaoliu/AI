cmake_minimum_required(VERSION 3.0)

if(WIN32)
    add_definitions(-DFACELIB_EXPORTS -D_STLP_USE_NEWALLOC)
    # avoid error C2589 on Windows
    add_definitions(-DNOMINMAX)
endif()

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

file(GLOB_RECURSE srcs "*.cpp")
file(GLOB_RECURSE hdrs "*.h*")

add_executable(WarpTriangle  ${srcs} ${hdrs})
target_link_libraries(WarpTriangle ${OpenCV_LIBS})
