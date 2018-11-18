#Cmake 最低版本号要求
cmake_minimum_required(VERSION 3.0)

project(ppp)

#添加C++11编译器支持
add_compile_options(-std=c++11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_BUILD_TYPE "debug")

#生成可调试程序
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")

#添加cpp文件并输出
file(GLOB src_files code/*.cpp)
message(${src_files})

add_executable(ppp ${src_files})

#指定可执行文件输出路径
set(EXECUTABLE_OUTPUT_PATH bin)
