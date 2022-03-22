cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(graphics_engine)
#SET(VTK_DIR <location of your vtk build directory>)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(graphics_engine graphics_engine)

if(VTK_LIBRARIES)
  target_link_libraries(graphics_engine ${VTK_LIBRARIES})
else()
  target_link_libraries(graphics_engine vtkHybrid)
endif()

