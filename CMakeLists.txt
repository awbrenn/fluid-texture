cmake_minimum_required(VERSION 2.8.4)
project(fluid_simulator)

set(SOURCE_FILES fluid_simulator.cpp cfd.h cfd.cpp cfdUtility.h)


if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    include_directories("/usr/local/include")
    find_library(OIIO "OpenImageIO" "/usr/local/lib")
    find_library(FOUNDATION "Foundation")
    find_library(GLUT "GLUT")
    find_library(OPENGL "OpenGL")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -Wall")
    include_directories("/usr/include")
    find_library(OIIO "OpenImageIO" "/usr/lib")
    find_library(GLUT "glut")
    find_library(GL "GL")
    find_library(GLU "GLU")
    find_library(GLEW "glew")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

add_executable(fluid_simulator ${SOURCE_FILES})

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    target_link_libraries(fluid_simulator ${OIIO} ${FOUNDATION} ${GLUT} ${OPENGL})
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    target_link_libraries(fluid_simulator ${OIIO} ${GLUT} ${GL} ${GLU} ${GLEW})
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
