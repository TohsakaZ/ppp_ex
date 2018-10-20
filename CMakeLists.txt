#Cmake 最低版本号要求
cmake_minimum_required(VERSION 3.0)

project(ppp)

#添加C++11编译器支持
add_compile_options(-std=c++11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#添加cpp文件并输出
file(GLOB src_files code/*.cpp)
message(${src_files})

add_executable(ppp ${src_files})

#指定可执行文件输出路径
set(EXECUTABLE_OUTPUT_PATH bin)