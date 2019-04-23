# 声明要求的 cmake 最低版本
cmake_minimum_required( VERSION 3.14 ) 

# 声明一个 cmake 工程
project( MonoCameraCalibration )

#PACKAGE
find_package( OpenCV REQUIRED )


#Customize Variables
set(OPENCV_MONO_CALIBRATION_SRC OpenCV_MonoCamera_Calibration.cpp)



# 设置编译模式
set( CMAKE_BUILD_TYPE "Debug" )


# 添加一个可执行程序
# 语法：add_executable( 程序名 源代码文件 ）
add_executable( ${PROJECT_NAME} ${OPENCV_MONO_CALIBRATION_SRC} )


# 将库文件链接到可执行程序上
target_link_libraries( ${PROJECT_NAME} ${OpenCV_LIBS} )
